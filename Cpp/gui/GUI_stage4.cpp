#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

#include <LBFGS.h>

namespace zebrafish {

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 4: Optimization

void GUI::DrawStage4() {

    // Visualize optimized points
    static int optimPointSize = 7;
    if (showOptimizedPoints) {

        static double optimEnergyThres_cache = -1;
        if (optimEnergyThres != optimEnergyThres_cache) {
            // Update the points to visualize when "optimEnergyThres" has changed or 
            // "Optimization" has been re-computed
            UpdateOptimPointLoc();
            optimEnergyThres_cache = optimEnergyThres;
        }

        viewer.data().point_size = optimPointSize;
        Eigen::MatrixXd pointColor(1, 3);
        pointColor << 0.87, 0.33, 0.33;

        if (optimPointLoc.rows() > 0) {
            // show optimized points
            viewer.data().add_points(
                optimPointLoc,
                pointColor
            );
        }

        ////// DEBUG ONLY //////
        Eigen::MatrixXd tempLoc;
        tempLoc.resize(3, 3);
        tempLoc << 0, 0, 1, 
                   imgCols, imgRows, 1, 
                   imgCols-1, imgRows-1, 1;
        Eigen::MatrixXd pointColor_t(1, 3);
        pointColor << 0.33, 0.83, 0.33;
        viewer.data().add_points(tempLoc, pointColor_t);
        ////// DEBUG ONLY //////
    }

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::Text("Optimization");
    if (ImGui::Button("Start Optimization")) {

        Optimization();
    }

    ImGui::Text("Stage 4");
}


////////////////////////////////////////////////////////////////////////////////////////
// Optimization


void GUI::Optimization() {

    /////////////////////////////////////////////////
    // prepare LBFGS
    LBFGSpp::LBFGSParam<double> param;
    param.epsilon = 1e-4;
    param.max_iterations = 15;

    // prepare sampleInput & sampleOutput for Newtons
    const int N = pointRecord.num;
    pointRecord.optimization.resize(N, 6);  // x, y, z, r, energy, iteration
    logger().info("Optimization #starting points = {}", N);

    // Optimization
    logger().info(">>>>>>>>>> Before optimization >>>>>>>>>>");
    tbb::parallel_for( tbb::blocked_range<int>(0, N),
        //////////////////////////////////////
        // lambda function for parallel_for //
        //////////////////////////////////////
        [this/*.grid_search, .optimization, &bsplineSolver*/, &param]
        (const tbb::blocked_range<int> &r) {

        // NOTE: LBFGSSolver is NOT thread safe. This must be instantiated for every thread
        LBFGSpp::LBFGSSolver<double> solver(param);

        // NOTE: the "variable count" used by "Autodiff" will be stored in 
        //       thread-local memory, so this must be set for every thread
        DiffScalarBase::setVariableCount(3);

        for (int ii = r.begin(); ii != r.end(); ++ii) {    

                Eigen::VectorXd vec(3, 1);
                vec(0) = pointRecord.grid_search(ii, 0);  // x
                vec(1) = pointRecord.grid_search(ii, 1);  // y
                vec(2) = pointRecord.grid_search(ii, 3);  // r
                double res;

                ///////////////////////////////////
                // lambda function for optimizer //
                ///////////////////////////////////
                auto func = [this/*.sampleInput_Newton, .bsplineSolver*/, ii]
                (const Eigen::VectorXd& x, Eigen::VectorXd& grad) {

                        DScalar ans;

                        if (!cylinder::IsValid(bsplineSolver, x(0), x(1), pointRecord.grid_search(ii, 2), x(2), 3)) {
                            grad.setZero();
                            return 1.0;
                        }
                        cylinder::EvaluateCylinder(bsplineSolver, DScalar(0, x(0)), DScalar(1, x(1)), pointRecord.grid_search(ii, 2), DScalar(2, x(2)), 3, ans);
                        grad.resize(3, 1);
                        grad = ans.getGradient();
                        return ans.getValue();
                    };
                // NOTE: the template of "solver.minimize" does not accept a temprary variable (due to non-const argument)
                //       so we define a "func" and pass it in
                int it = solver.minimize(func, vec, res);
                ///////////////////////////////////

                pointRecord.optimization(ii, 0) = vec(0);  // x
                pointRecord.optimization(ii, 1) = vec(1);  // y
                pointRecord.optimization(ii, 2) = pointRecord.grid_search(ii, 2);  // z
                pointRecord.optimization(ii, 3) = vec(2);  // r
                pointRecord.optimization(ii, 4) = res;     // energy
                pointRecord.optimization(ii, 5) = it;      // iteration
        }
    });
    logger().info("<<<<<<<<<< After optimization <<<<<<<<<<");
}


void GUI::UpdateOptimPointLoc() {

    const int N = pointRecord.num;

    optimPointLoc.resize(N, 3);
    if (N == 0) return;

    optimPointLoc.col(0) = pointRecord.optimization.col(1).array() + 0.5;
    optimPointLoc.col(1) = (imgRows-0.5) - pointRecord.optimization.col(0).array();
    optimPointLoc.col(2) = pointRecord.optimization.col(2);

    logger().info("Optimization resultant points updated: total number = {}", N);
}

}  // namespace zebrafish
