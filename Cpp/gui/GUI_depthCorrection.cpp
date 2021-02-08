#include <zebrafish/GUI.h>
#include <zebrafish/TiffReader.h>
#include <zebrafish/FileDialog.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/Quantile.h>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

#include <LBFGS.h>
#include <string>
#include <sstream>
#include <algorithm>
#include <set>


namespace zebrafish {

namespace {

double GetMean(const Eigen::MatrixXd &energy, int i, int j) {
// Find the mean of (2*l+1) consecutive energies, excluding the min and max

    const int l = 2;
    int t1 = std::max(0, j-l);
    int t2 = std::min(int(energy.cols()-1), j+l);

    double minE = 1.0, maxE = -1.0, sum = 0.0;
    for (int t=t1; t<=t2; t++) {
        if (energy(i, t) > maxE) maxE = energy(i, t);
        if (energy(i, t) < minE) minE = energy(i, t);
        sum += energy(i, t);
    }

    return (sum - minE - maxE) / double(t2-t1-1);
}


void SmoothCurve(Eigen::MatrixXd &energy_cache) {

    const int rows = energy_cache.rows();
    const int cols = energy_cache.cols();
    Eigen::MatrixXd tmpEnergy(rows, cols);

    for (int i=0; i<rows; i++)
        for (int j=0; j<cols; j++) {
            tmpEnergy(i, j) = GetMean(energy_cache, i, j);
        }

    energy_cache = tmpEnergy;
}


void CalcLaplacian(const Eigen::MatrixXd &energy_cache, Eigen::MatrixXd &secondDerivative) {

    const int rows = energy_cache.rows();
    const int cols = energy_cache.cols();

    for (int i=0; i<rows; i++)
        for (int j=0; j<cols; j++) {

            if (j==0 || j==cols-1) {
                secondDerivative(i, j) = 0.0;
                continue;  // no padding
            }
            secondDerivative(i, j) = energy_cache(i, j) * 2.0 - energy_cache(i, j-1) - energy_cache(i, j+1);
        }
}


bool ValidDerivative(const Eigen::MatrixXd &secondDerivative, int i, int j) {
// whether the second derivative is small near the minimum point?

    const int cols = secondDerivative.cols();
    int colStart = std::max(0, j-8);
    int colEnd = std::min(cols-1, j+8);
    int count = 0;

    for (int col=colStart; col<=colEnd; col++) {
        if (std::fabs(secondDerivative(i, col)) > 5e-3) count++;
    }
    if (count >= 2)
        return false;
    else
        return true;
}


double FindDepthFromMesh(const Eigen::MatrixXi &markerMeshArray, const Eigen::MatrixXd &depthArray, int index) {

    const int N = markerMeshArray.rows();
    double depth = 0;
    std::set<int> indexSet;

    for (int i=0; i<N; i++) {

        if (markerMeshArray(i, 0) == index || markerMeshArray(i, 1) == index || markerMeshArray(i, 2) == index) {
            // add two adjacent markers' depths
            if (markerMeshArray(i, 0) != index) indexSet.insert(markerMeshArray(i, 0));
            if (markerMeshArray(i, 1) != index) indexSet.insert(markerMeshArray(i, 1));
            if (markerMeshArray(i, 2) != index) indexSet.insert(markerMeshArray(i, 2));
        }
    }

    /*
    // DEBUG PURPOSE
    for (auto it=indexSet.begin(); it!=indexSet.end(); it++) 
        std::cerr << *it << " ";
    std::cerr << std::endl;
    std::cerr << "markerMeshArray size" << std::endl;
    std::cerr << markerMeshArray.size() << std::endl;
    std::cerr << "depth Array" << std::endl;
    std::cerr << depthArray.transpose() << std::endl;
    */

    for (auto it=indexSet.begin(); it!=indexSet.end(); it++) {
        depth += depthArray(*it);
    }

    if (depth == 0 || indexSet.empty()) {
        depth = 0;
        logger().warn("No adjacent marker found for a bad marker. Marker index {}", index);
        std::cerr << "No adjacent marker found for a bad marker. Marker index " << index << std::endl;
    } else {
        depth /= double(indexSet.size());
        // do not divide by zero
    }

    return depth;
}

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Special module about depth correction

bool GUI::MarkerDepthCorrection(int frameIdx, int depthNum, double depthGap, bool logEnergy) {
// try different z's and pick the one with minimal energy

    const int N = markerArray[frameIdx].num;
    const int M = depthNum * 2 + 1;
    if (N == 0) return true;

    /////////////////////////////////////////////////
    // prepare LBFGS
    LBFGSpp::LBFGSParam<double> param;
    param.epsilon = optimEpsilon;
    param.max_iterations = optimMaxIt;

    // prepare for parallel optimization
    Eigen::MatrixXd x_cache, y_cache, z_cache, r_cache, energy_cache;
    Eigen::VectorXd depthArray(M);
    x_cache.resize(N, M);
    y_cache.resize(N, M);
    z_cache.resize(N, M);
    r_cache.resize(N, M);
    energy_cache.resize(N, M);
    for (int i=-depthNum; i<=depthNum; i++) {
        depthArray[i+depthNum] = depthGap * i;
    }
    logger().info("Depth correction for markers in frame {}  #starting points = {}, searching gap = {}, searching num = {}", frameIdx, N, depthGap, depthNum);

    // Optimization
    tbb::parallel_for( tbb::blocked_range<int>(0, N),
        //////////////////////////////////////
        // lambda function for parallel_for //
        //////////////////////////////////////
        [this/*.markerArray[frameIdx], .bsplineArray[frameIdx]*/, &param, frameIdx, M, &depthArray, 
         &x_cache, &y_cache, &z_cache, &r_cache, &energy_cache]
        (const tbb::blocked_range<int> &r) {

        // NOTE: LBFGSSolver is NOT thread safe. This must be instantiated for every thread
        LBFGSpp::LBFGSSolver<double> solver(param);

        // NOTE: the "variable count" used by "Autodiff" will be stored in 
        //       thread-local memory, so this must be set for every thread
        DiffScalarBase::setVariableCount(3);

        for (int ii = r.begin(); ii != r.end(); ++ii) {    

            // iterate over all z's
            for (int jj=0; jj<M; jj++) {

                Eigen::VectorXd vec(3, 1);
                vec(0) = markerArray[frameIdx].loc(ii, 0);  // x
                vec(1) = markerArray[frameIdx].loc(ii, 1);  // y
                vec(2) = markerArray[frameIdx].loc(ii, 3);  // r
                double res;

                ///////////////////////////////////
                // lambda function for optimizer //
                ///////////////////////////////////
                auto func = [this/*.markerArray[frameIdx], .bsplineArray[frameIdx]*/, ii, jj, frameIdx, &depthArray]
                (const Eigen::VectorXd& x, Eigen::VectorXd& grad) {

                        DScalar ans;
                        double z_ = markerArray[frameIdx].loc(ii, 2) + depthArray[jj];  // add z correction

                        if (!cylinder::IsValid(bsplineArray[frameIdx], x(0), x(1), z_, x(2), cylinder::H)) {
                            grad.setZero();
                            return 1.0;
                        }
                        cylinder::EvaluateCylinder(bsplineArray[frameIdx], DScalar(0, x(0)), DScalar(1, x(1)), z_, DScalar(2, x(2)), cylinder::H, ans, reverseColor);
                        grad.resize(3, 1);
                        grad = ans.getGradient();
                        return ans.getValue();
                    };
                // NOTE: the template of "solver.minimize" does not accept a temprary variable (due to non-const argument)
                //       so we define a "func" and pass it in
                int it = solver.minimize(func, vec, res);
                ///////////////////////////////////

                x_cache(ii, jj) = vec(0);    // x
                y_cache(ii, jj) = vec(1);    // y
                z_cache(ii, jj) = markerArray[frameIdx].loc(ii, 2) + depthArray[jj];  // z
                r_cache(ii, jj) = vec(2);    // r
                energy_cache(ii, jj) = res;  // energy
            }
        }
    });  // end of tbb::parallel_for

    //////////////////////////////////////////////////////////////////////////////
    // find min for each marker

    bool res = true;
    int minColIdx, correctedCount = 0;
    double minEnergy;
    std::vector<bool> derivativeTable(N, true);
    bool badDerivative = false;

    // DEBUG PURPOSE
    if (logEnergy) {
        using namespace std;
        cout << " >>>>>>> depth correction energy_cache >>>>>>>" << endl;
        cout << "Frame = " << frameIdx << endl << endl;
        cout << energy_cache << endl << endl;
        /*
        cout << "x cache" << endl << x_cache << endl;
        cout << "y cache" << endl << y_cache << endl;
        cout << "z cache" << endl << z_cache << endl;
        cout << "r cache" << endl << r_cache << endl;
        */
    }
    // DEBUG PURPOSE

    // smoothen the energy cache twice
        /// Why twice? This is an empirical decision.
        /// (1) for most (almost all) cases, do it twice will not change anything
        /// (2) there are very rare cases where do it once will fail
    SmoothCurve(energy_cache);
    SmoothCurve(energy_cache);
    // 2nd derivative
    Eigen::MatrixXd secondDerivative(energy_cache.rows(), energy_cache.cols());
    CalcLaplacian(energy_cache, secondDerivative);

    for (int i=0; i<N; i++) {

        const double thresDist = optimMaxXYDisp;  // cannot move over "thresDist" pixels in xy plane
        minEnergy = 1.0;  // reset
        for (int j=0; j<M; j++) {

            // valid?
            Eigen::Vector2d xyDist;
            xyDist << x_cache(i, j)-markerArray[frameIdx].loc(i, 0), y_cache(i, j)-markerArray[frameIdx].loc(i, 1);
            if (xyDist.norm() > thresDist) continue;
            // does this depth lead to a smaller energy?
            if (energy_cache(i, j) < minEnergy) {
                minEnergy = energy_cache(i, j);
                minColIdx = j;
            }
        }

        if (minEnergy == 1.0) {
            // this should not happen
            char errorMsg[100];
            std::sprintf(errorMsg, "> [warning] Depth correction min energy 1.0 exception: Frame %d, Marker index %d at [%.2f, %.2f, %.2f].", frameIdx, i, markerArray[frameIdx].loc(i, 0), markerArray[frameIdx].loc(i, 1), markerArray[frameIdx].loc(i, 2));
            logger().warn(errorMsg);
            std::cerr << errorMsg << std::endl;
            res = false;
        } else {
            // minEnergy < 1.0

            if (minColIdx != depthNum) {
                // if not the original z
                /*
                logger().debug("frameIdx = {} | old x={} y={} z={} r={} e={} | new x={} y={} z={} r={} e={}", 
                        frameIdx, markerArray[frameIdx].loc(i, 0), markerArray[frameIdx].loc(i, 1), markerArray[frameIdx].loc(i, 2), markerArray[frameIdx].loc(i, 3), markerArray[frameIdx].energy(i), 
                        x_cache(i, minColIdx), y_cache(i, minColIdx), z_cache(i, minColIdx), r_cache(i, minColIdx), energy_cache(i, minColIdx));
                */
                // correctedCount++;
            }

            // min is at end point?
            if (depthNum > 0 && (minColIdx == 0 || minColIdx == 1 || minColIdx == M-1 || minColIdx == M-2)) {
                char errorMsg[200];
                std::sprintf(errorMsg, "> [warning] Reached depth correction search range limit. Consider using a larger search range. Frame %d, Marker index %d at [%.2f, %.2f, %.2f].", frameIdx, i, markerArray[frameIdx].loc(i, 0), markerArray[frameIdx].loc(i, 1), markerArray[frameIdx].loc(i, 2));
                logger().warn(errorMsg);
                std::cerr << errorMsg << std::endl;
                res = false;
            }
            // derivative looks good?
            if (depthNum > 0 && !ValidDerivative(secondDerivative, i, minColIdx)) {
                derivativeTable[i] = false;
                badDerivative = true;
                char warnMsg[200];
                std::sprintf(warnMsg, "> [note] This marker's second derivative in depth search is abnormal: Frame %d, Marker index %d at [%.2f, %.2f, %.2f].", frameIdx, i, markerArray[frameIdx].loc(i, 0), markerArray[frameIdx].loc(i, 1), markerArray[frameIdx].loc(i, 2));
                logger().info(warnMsg);
                std::cerr << warnMsg << std::endl;
            }

            markerArray[frameIdx].loc(i, 0) = x_cache(i, minColIdx);
            markerArray[frameIdx].loc(i, 1) = y_cache(i, minColIdx);
            markerArray[frameIdx].loc(i, 2) = z_cache(i, minColIdx);
            markerArray[frameIdx].loc(i, 3) = r_cache(i, minColIdx);
            markerArray[frameIdx].energy(i) = energy_cache(i, minColIdx);
        }
    }

    /////////////////////////////////////////////////////////////////////
    // Do not trust the depth of the markers whose derivative is abnormal
        /// Why frameIdx>0? For the first frame, the topology relation has not been inferred yet
        ///                 It is useless to run "FindDepthFromMesh"
    if (badDerivative && frameIdx > 0) {
        for (int i=0; i<N; i++) {

            if (derivativeTable[i]) continue;  // if this marker's derivative is OK

            // replace its depth by the mean of depths of adjacent markers
            double newDepth = FindDepthFromMesh(markerMeshArray, markerArray[frameIdx].loc.col(2), i);
            if (newDepth > 0) markerArray[frameIdx].loc(i, 2) = newDepth;
        }
        // optimize all marker with fixed depth
        // this will not be a recursion
        try {
            MarkerDepthCorrection(frameIdx, 0, 0, false);
        } catch (const std::exception &e) {
            logger().error("Fatal error: MarkerDepthCorrection for bad derivative markers");
            std::cerr << "Fatal error: MarkerDepthCorrection for bad derivative markers" << std::endl;
        }
    }

    return res;
}

}  // namespace zebrafish
