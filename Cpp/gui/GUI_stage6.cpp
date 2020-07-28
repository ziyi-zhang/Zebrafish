#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>


namespace zebrafish {

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 6: Iterative Closest Point

void GUI::DrawStage6() {

    if (stage5to6Flag) {

        FinalizeClusterLoc();
        UpdateMarkerPointLoc();
        stage5to6Flag = false;
    }

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
        Eigen::MatrixXd pointColor_t(1, 3);
        pointColor << 0.33, 0.83, 0.33;
        viewer.data().add_points(tempLoc, pointColor_t);
        ////// DEBUG ONLY //////
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Iterative Closest Point", ImGuiTreeNodeFlags_DefaultOpen)) {

        
        ImGui::Text("TBD");
    }

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::Text("Stage 6: ICP");
}


////////////////////////////////////////////////////////////////////////////////////////
// 

void GUI::UpdateMarkerPointLoc() {

    const int N = markerRecord.num;
    Eigen::MatrixXd tempLoc;

    markerPointLoc.resize(N, 3);
    
    tempLoc.resize(N, 3);
    tempLoc.col(0) = markerRecord.loc.col(0);  // x
    tempLoc.col(1) = markerRecord.loc.col(1);  // y
    tempLoc.col(2) = markerRecord.loc.col(2);  // z

    markerPointLoc.col(0) = tempLoc.col(1).array() + 0.5;
    markerPointLoc.col(1) = (imgRows-0.5) - tempLoc.col(0).array();
    markerPointLoc.col(2) = tempLoc.col(2);

    logger().info("   [Visualization] Marker clusters location updated: total number = {}", N);
}

}  // namespace zebrafish
