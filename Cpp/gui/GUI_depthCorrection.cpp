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


namespace zebrafish {

namespace {

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

                        if (!cylinder::IsValid(bsplineArray[frameIdx], x(0), x(1), z_, x(2), 3)) {
                            grad.setZero();
                            return 1.0;
                        }
                        cylinder::EvaluateCylinder(bsplineArray[frameIdx], DScalar(0, x(0)), DScalar(1, x(1)), z_, DScalar(2, x(2)), 3, ans, reverseColor);
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
            std::sprintf(errorMsg, "[warning] Depth correction exception: Frame %d, Marker index %d at [%.2f, %.2f, %.2f].", frameIdx, i, markerArray[frameIdx].loc(i, 0), markerArray[frameIdx].loc(i, 1), markerArray[frameIdx].loc(i, 2));
            logger().warn(errorMsg);
            std::cerr << errorMsg << std::endl;
            res = false;
        } else {
            if (minColIdx != depthNum) {
                // if not the original z
                /*
                logger().debug("frameIdx = {} | old x={} y={} z={} r={} e={} | new x={} y={} z={} r={} e={}", 
                        frameIdx, markerArray[frameIdx].loc(i, 0), markerArray[frameIdx].loc(i, 1), markerArray[frameIdx].loc(i, 2), markerArray[frameIdx].loc(i, 3), markerArray[frameIdx].energy(i), 
                        x_cache(i, minColIdx), y_cache(i, minColIdx), z_cache(i, minColIdx), r_cache(i, minColIdx), energy_cache(i, minColIdx));
                */
                // correctedCount++;
            }
            if (depthNum > 0 && (minColIdx == 0 || minColIdx == M-1)) {
                logger().warn("Reached depth correction search range limit. Consider using a larger search range.");
            }

            markerArray[frameIdx].loc(i, 0) = x_cache(i, minColIdx);
            markerArray[frameIdx].loc(i, 1) = y_cache(i, minColIdx);
            markerArray[frameIdx].loc(i, 2) = z_cache(i, minColIdx);
            markerArray[frameIdx].loc(i, 3) = r_cache(i, minColIdx);
            markerArray[frameIdx].energy(i) = energy_cache(i, minColIdx);
        }
    }

    return res;
}

}  // namespace zebrafish
