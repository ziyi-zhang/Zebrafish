#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

#include <LBFGS.h>
#include <string>


namespace zebrafish {

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 8: Dispalcement

void GUI::DrawStage8() {

    // Visualize marker cluster points
    static int pointSize = 7;
    if (showMarkerPoints) {

        viewer.data().point_size = pointSize;

        if (!markerPointLocArray.empty()) {
            // show optimized markers
            if (!manualOverrideMarkerVis) {
                // show markers in the frame that is currently focused
                viewer.data().add_points(
                    markerPointLocArray[frameToShow],
                    markerPointColor
                );
            } else {
                // show markers in manually selected frames
                for (int i=0; i<markerPointStatusArray.rows(); i++) {
                    if (!markerPointStatusArray(i)) continue;
                    viewer.data().add_points(
                        markerPointLocArray[i],
                        markerPointColor
                    );
                }
            }
        }

        ////// DEBUG ONLY //////
        Eigen::MatrixXd tempLoc;
        tempLoc.resize(3, 3);
        tempLoc << 0, 0, 1, 
                   imgCols, imgRows, 1, 
                   imgCols-1, imgRows-1, 1;
        Eigen::MatrixXd debugPointColor(1, 3);
        debugPointColor << 0.33, 0.83, 0.33;
        viewer.data().add_points(tempLoc, debugPointColor);
        ////// DEBUG ONLY //////
    }

    // Visualize meshes
    DrawMarkerMesh();

    // visualize optical flow points
    // FIXME: DO NOT DO THIS
    if (showOpticalFlow) {

        viewer.data().point_size = pointSize;
        Eigen::MatrixXd pointColor(1, 3);
        pointColor << 0.91, 0.69, 0.01;

        if (!opticalFlowCorrection.empty() && frameToShow < opticalFlowCorrection.size()) {
            // show optical flow corrected location in the frame that is currently focused
            Eigen::MatrixXd temp(markerPointLocArray[frameToShow].rows(), markerPointLocArray[frameToShow].cols());
            temp.col(0) = markerPointLocArray[frameToShow].col(0) + opticalFlowCorrection[frameToShow].col(1);
            temp.col(1) = markerPointLocArray[frameToShow].col(1) - opticalFlowCorrection[frameToShow].col(0);
            temp.col(2) = markerPointLocArray[frameToShow].col(2) + opticalFlowCorrection[frameToShow].col(2);
            viewer.data().add_points(
                temp,
                pointColor
            );
        }
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Calculate Displacement", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::SliderInt("Depth correction search range", &depthCorrectionNum, 0, 12, "%d * gap");
        ImGui::SliderFloat("Depth correction gap", &depthCorrectionGap, 0, 2.0, "%.3f pixels");
        if (ImGui::Button("Find new locations")) {
            OptimizeAllFrames();
            logger().debug("   <button> Find new locations");
        }
    }
}


////////////////////////////////////////////////////////////////////////////////////////
// Optimization


void GUI::OptimizeAllFrames() {

    int currentFrame;

    for (currentFrame=0; currentFrame<currentLoadedFrames-1; currentFrame++) {

        OptimizeOneFrame(currentFrame);
        // depth correction for the frame that was just updated
        MarkerDepthCorrection(currentFrame + 1, depthCorrectionNum, depthCorrectionGap);
    }

    // update visualization variable
    UpdateMarkerPointLocArray();
}


void GUI::OptimizeOneFrame(int prevFrameIdx) {

    /////////////////////////////////////////////////
    // prepare LBFGS
    LBFGSpp::LBFGSParam<double> param;
    param.epsilon = optimEpsilon;
    param.max_iterations = optimMaxIt;

    // prepare for parallel optimization
    const int N = markerArray[prevFrameIdx].num;
    markerArray[prevFrameIdx+1].num = N;

    logger().info("Optimization from frame {} to {}", prevFrameIdx, prevFrameIdx+1);
    logger().info("Optimization #starting points = {}", N);

    // Optimization
    logger().info(">>>>>>>>>> Before optimization >>>>>>>>>>");
    tbb::parallel_for( tbb::blocked_range<int>(0, N),
        //////////////////////////////////////
        // lambda function for parallel_for //
        //////////////////////////////////////
        [this/*.markerArray[?], &bsplineArray[?], .opticalFlowCorrection*/, &param, prevFrameIdx]
        (const tbb::blocked_range<int> &r) {

        // NOTE: LBFGSSolver is NOT thread safe. This must be instantiated for every thread
        LBFGSpp::LBFGSSolver<double> solver(param);

        // NOTE: the "variable count" used by "Autodiff" will be stored in 
        //       thread-local memory, so this must be set for every thread
        DiffScalarBase::setVariableCount(3);

        for (int ii = r.begin(); ii != r.end(); ++ii) {    

                Eigen::VectorXd vec(3, 1);
                vec(0) = markerArray[prevFrameIdx].loc(ii, 0) + opticalFlowCorrection[prevFrameIdx](ii, 0);  // x
                vec(1) = markerArray[prevFrameIdx].loc(ii, 1) + opticalFlowCorrection[prevFrameIdx](ii, 1);  // y
                vec(2) = markerArray[prevFrameIdx].loc(ii, 3);  // r
                double res;

                ///////////////////////////////////
                // lambda function for optimizer //
                ///////////////////////////////////
                auto func = [this/*.markerArray[?], .bsplineArray[?], .opticalFlowCorrection*/, ii, prevFrameIdx]
                (const Eigen::VectorXd& x, Eigen::VectorXd& grad) {

                        DScalar ans;

                        if (!cylinder::IsValid(bsplineArray[prevFrameIdx+1], x(0), x(1), markerArray[prevFrameIdx].loc(ii, 2) + opticalFlowCorrection[prevFrameIdx](ii, 2), x(2), 3)) {
                            grad.setZero();
                            return 1.0;
                        }
                        cylinder::EvaluateCylinder(bsplineArray[prevFrameIdx+1], DScalar(0, x(0)), DScalar(1, x(1)), markerArray[prevFrameIdx].loc(ii, 2) + opticalFlowCorrection[prevFrameIdx](ii, 2), DScalar(2, x(2)), 3, ans);
                        grad.resize(3, 1);
                        grad = ans.getGradient();
                        return ans.getValue();
                    };
                // NOTE: the template of "solver.minimize" does not accept a temprary variable (due to non-const argument)
                //       so we define a "func" and pass it in
                int it = solver.minimize(func, vec, res);
                ///////////////////////////////////

                markerArray[prevFrameIdx+1].loc(ii, 0) = vec(0);  // x
                markerArray[prevFrameIdx+1].loc(ii, 1) = vec(1);  // y
                markerArray[prevFrameIdx+1].loc(ii, 2) = markerArray[prevFrameIdx].loc(ii, 2) + opticalFlowCorrection[prevFrameIdx](ii, 2);  // z
                markerArray[prevFrameIdx+1].loc(ii, 3) = vec(2);  // r
                markerArray[prevFrameIdx+1].energy(ii) = res;     // energy
        }
    });
    logger().info("<<<<<<<<<< After optimization <<<<<<<<<<");
}

}  // namespace zebrafish
