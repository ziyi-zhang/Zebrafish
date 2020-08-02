#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/TiffReader.h>
#include <zebrafish/Quantile.h>
#include <zebrafish/OpticalFlow.h>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

#include <string>


namespace zebrafish {

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 7: Optical Flow

void GUI::DrawStage7() {

    // Visualize marker cluster points
    static int pointSize = 7;
    if (showMarkerPoints) {

        viewer.data().point_size = pointSize;
        Eigen::MatrixXd pointColor(1, 3);
        pointColor << 0.99, 0.41, 0.01;

        if (!markerPointLocArray.empty()) {
            // show optimized markers
            if (!manualOverrideMarkerVis) {
                // show markers in the frame that is currently focused
                viewer.data().add_points(
                    markerPointLocArray[frameToShow],
                    pointColor
                );
            } else {
                // show markers in manually selected frames
                for (int i=0; i<markerPointStatusArray.rows(); i++) {
                    if (!markerPointStatusArray(i)) continue;
                    viewer.data().add_points(
                        markerPointLocArray[i],
                        pointColor
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

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Load frames", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::SliderInt("Desired frame number", &desiredFrames, 2, ttlFrames, "%d frames");
        if (ImGui::Button("Load subsequent frames")) {
            LoadSubsequentFrames();
        }
        if (ImGui::Button("Compute B-spline for all")) {
            ComputeBsplineForAllFrames();
        }
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Optical Flow", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::InputDouble("alpha", &opticalFlowAlpha);
        ImGui::SliderInt("iteration", &opticalFlowIter, 1, 100);
        if (ImGui::Button("Run Optical Flow")) {
            RunOpticalFlow();
        }
    }

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::Text("Stage 7: Optical Flow");
}


////////////////////////////////////////////////////////////////////////////////////////
// Load heler


void GUI::LoadSubsequentFrames() {

    // reserve space
    imgData.resize(desiredFrames);
    compressedImgTextureArray.resize(desiredFrames);
    bsplineArray.resize(desiredFrames);
        // initialize with the markers in the first frame
    markerArray.resize(desiredFrames, markerArray[0]);
        // initialize with zero matrices
    opticalFlowCorrection.resize(desiredFrames - 1);
    for (int i=0; i<desiredFrames - 1; i++) {
        opticalFlowCorrection[i] = Eigen::MatrixXd::Zero(markerArray[0].num, 3);
    }

    // only support loading one channel
    std::vector<bool> channelVec(channelPerSlice, false);
    channelVec[channelToLoad] = true;

    ReadTif(imagePath, layerPerImg, channelVec, desiredFrames, imgData, r0, c0, r1, c1);
    currentLoadedFrames = desiredFrames;

    // quantile trim
    for (int i=0; i<currentLoadedFrames; i++) {
        double normalizeQuantileRes = QuantileImage(imgData[i], normalizeQuantile);
        NormalizeImage(imgData[i], normalizeQuantileRes);
        ComputeCompressedTexture(imgData[i], i);
        logger().info("[Quantile Trim] 3D Image (index = {}) normalized with normalizeQuantile =  {:.4f}  thres =  {:.4f}", i, normalizeQuantile, normalizeQuantileRes);
    }
}


void GUI::ComputeBsplineForAllFrames() {

    const int bsplineDegree = 2;

    for (int i=1; i<currentLoadedFrames; i++)
        bsplineArray[i].SetResolution(resolutionX, resolutionY, resolutionZ);

    // Parallel B-spline computation
    logger().info(">>>>>>>>>> Before B-spline >>>>>>>>>>");
    logger().info("B-spline #frames = {}", currentLoadedFrames-1);
    tbb::parallel_for( tbb::blocked_range<int>(1, currentLoadedFrames),
        [this/*.bsplineArray[ii], .imgData[ii]*/, bsplineDegree](const tbb::blocked_range<int> &r) {

            for (int ii = r.begin(); ii != r.end(); ++ii) {
                bsplineArray[ii].CalcControlPts(imgData[ii], 0.7, 0.7, 0.7, bsplineDegree);
            }
        });
    logger().info("<<<<<<<<<< After B-spline <<<<<<<<<<");

    // DEBUG PURPOSE
    // used to test the correctness of parallel B-spline
    /*
    logger().debug("reference value = {}", markerArray[0].energy(0));
    double res;
    for (int i=0; i<currentLoadedFrames; i++) {
        cylinder::EvaluateCylinder(bsplineArray[i], markerArray[0].loc(0, 0), markerArray[0].loc(0, 1), markerArray[0].loc(0, 2), markerArray[0].loc(0, 3), 3.0, res);
        logger().debug("frame {} res = {}", i, res);
    }
    */
}


////////////////////////////////////////////////////////////////////////////////////////
// Optical Flow


void GUI::RunOpticalFlow() {

    FlowVelocity_t ux, uy, uz;
    int prevFrameIdx, i;
    const int N = markerArray[0].num;
    int x, y, z;

    /// note: "opticalFlowCorrection" vector has been resized before calling this function

    for (prevFrameIdx=0; prevFrameIdx<currentLoadedFrames - 1; prevFrameIdx++) {

        OpticalFlow::RunOpticalFlow(imgData[prevFrameIdx], imgData[prevFrameIdx+1], opticalFlowAlpha, opticalFlowIter, ux, uy, uz);

        // calculate for each marker the correction displacement
        for (i=0; i<N; i++) {

            // alias
            x = std::round( markerArray[prevFrameIdx].loc(i, 0) );
            y = std::round( markerArray[prevFrameIdx].loc(i, 1) );
            z = std::round( markerArray[prevFrameIdx].loc(i, 2) );

            // correction
            /// (1) Use "nearest cell" instead of interpolation
            /// (2) x, y, z should always be in bound
            opticalFlowCorrection[prevFrameIdx](i, 0) = ux[z](x, y);
            opticalFlowCorrection[prevFrameIdx](i, 1) = uy[z](x, y);
            opticalFlowCorrection[prevFrameIdx](i, 2) = uz[z](x, y);
        }

        std::cout << prevFrameIdx << "=========" << std::endl;
        std::cout << opticalFlowCorrection[prevFrameIdx] << std::endl;
    }
}

}  // namespace zebrafish