#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/TiffReader.h>

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

        if (markerPointLoc.rows() > 0) {
            // show optimized cluster points
            viewer.data().add_points(
                markerPointLoc,
                pointColor
            );
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

    if (ImGui::CollapsingHeader("Optical Flow", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::SliderInt("Desired frame number", &desiredFrames, 2, ttlFrames, "%d frames");
        if (ImGui::Button("Load subsequent frames")) {
            LoadSubsequentFrames();
        }
    }

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::Text("Stage 7: Optical Flow");
}


////////////////////////////////////////////////////////////////////////////////////////
// Optical Flow


void GUI::LoadSubsequentFrames() {

    // reserve space
    imgData.resize(desiredFrames);
    compressedImgTextureArray.resize(desiredFrames);

    // only support loading one channel
    std::vector<bool> channelVec(channelPerSlice, false);
    channelVec[channelToLoad] = true;

    ReadTif(imagePath, layerPerImg, channelVec, desiredFrames, imgData, r0, c0, r1, c1);
    currentLoadedFrames = desiredFrames;
}

}  // namespace zebrafish
