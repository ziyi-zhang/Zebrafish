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

            UpdateOptimEnergyHist();
            // NOTE: This has been deprecated
            //       Do not use this cache
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

    if (ImGui::CollapsingHeader("Optimization", ImGuiTreeNodeFlags_DefaultOpen)) {
        
        if (ImGui::TreeNode("Advanced optim config")) {
            
            const float inputWidth = ImGui::GetWindowWidth() / 3.0;
            ImGui::PushItemWidth(inputWidth);
            ImGui::InputDouble("BFGS epsilon", &optimEpsilon);
            ImGui::InputDouble("BFGS max iter", &optimMaxIt);
            ImGui::PopItemWidth();
            
            ImGui::TreePop();
            ImGui::Separator();
        }

        if (ImGui::Button("Start Optimization")) {

            Optimization();
            UpdateOptimPointLoc();
            UpdateOptimEnergyHist();
        }

        ImGui::Checkbox("Show optimized locations", &showOptimizedPoints);

        ImGui::PushItemWidth(zebrafishWidth / 3.0);
        ImGui::SliderInt("Point Size", &optimPointSize, 1, 30);
        ImGui::PopItemWidth();

        // Histogram
        ImGui::Text("Histogram of energy");

        const float width = ImGui::GetWindowWidth() * 0.75f - 2;
        ImGui::PushItemWidth(width + 2);

        ImVec2 before = ImGui::GetCursorScreenPos();
        ImGui::PlotHistogram("", optimEnergyHist.hist.data(), optimEnergyHist.hist.size(), 0, NULL, 0, optimEnergyHist.hist.maxCoeff(), ImVec2(0, 80));
        ImVec2 after = ImGui::GetCursorScreenPos();
        after.y -= ImGui::GetStyle().ItemSpacing.y;

        float ratio = 0.0f;
        ImDrawList *drawList = ImGui::GetWindowDrawList();
        drawList->PushClipRectFullScreen();
        drawList->AddLine(
            ImVec2(before.x + width * ratio, before.y), 
            ImVec2(before.x + width * ratio, after.y), 
            IM_COL32(50, 205, 50, 255), 
            2.0f
        );
        drawList->PopClipRect();
        ImGui::PopItemWidth();

        ImGui::PushItemWidth(zebrafishWidth * 0.75);
        ImGui::PopItemWidth();
    }

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::Text("Stage 4: Optimization");
}


////////////////////////////////////////////////////////////////////////////////////////
// Optimization


void GUI::Optimization() {

    /////////////////////////////////////////////////
    // prepare LBFGS
    LBFGSpp::LBFGSParam<double> param;
    param.epsilon = optimEpsilon;
    param.max_iterations = optimMaxIt;

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
    if (N == 0) return;

    optimPointLoc.resize(N, 3);

    optimPointLoc.col(0) = pointRecord.optimization.col(1).array() + 0.5;
    optimPointLoc.col(1) = (imgRows-0.5) - pointRecord.optimization.col(0).array();
    optimPointLoc.col(2) = pointRecord.optimization.col(2);

    logger().info("[Visualization] Optimization resultant points updated: total number = {}", N);
}


void GUI::UpdateOptimEnergyHist() {

    Eigen::MatrixXd energyCol = pointRecord.optimization.col(4);
    double minValue = energyCol.minCoeff();
    double maxValue = energyCol.maxCoeff();
    const int N = pointRecord.num;
    assert(N > 0);

    const double epsilon = 0.0001;  // to make sure every number lies inside
    maxValue += epsilon;
    minValue -= epsilon;
    const double gap = (maxValue - minValue) / double(histBars);

    optimEnergyHist.hist = Eigen::MatrixXf::Zero(histBars, 1);
    optimEnergyHist.minValue = minValue;
    optimEnergyHist.maxValue = maxValue;

    for (int i=0; i<N; i++) {
        optimEnergyHist.hist( std::floor((energyCol(i) - minValue)/gap) )++;
    }
}

}  // namespace zebrafish
