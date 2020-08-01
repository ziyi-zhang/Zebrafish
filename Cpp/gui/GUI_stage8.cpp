#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>

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

    if (ImGui::CollapsingHeader("Calculate Displacement", ImGuiTreeNodeFlags_DefaultOpen)) {

        if (ImGui::Button("Find new locations")) {
            OptimizeAllFrames();
        }
    }

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::Text("Stage 8: Displacement");
}


////////////////////////////////////////////////////////////////////////////////////////
// Optimization


void GUI::OptimizeAllFrames() {

    int currentFrame;

    for (currentFrame=0; currentFrame<currentLoadedFrames-1; currentFrame++) {

        
    }
}

}  // namespace zebrafish
