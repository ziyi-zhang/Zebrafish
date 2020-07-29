#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/FileDialog.h>
#include <zebrafish/ICP.h>

#include <string>


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
        Eigen::MatrixXd debugPointColor(1, 3);
        debugPointColor << 0.33, 0.83, 0.33;
        viewer.data().add_points(tempLoc, debugPointColor);
        ////// DEBUG ONLY //////
    }

    // Visualize reference points
    if (showReferencePoints) {

        viewer.data().point_size = pointSize;
        Eigen::MatrixXd pointColor(1, 3);
        pointColor << 0.41, 0.41, 0.99;

        if (refPointLoc.rows() > 0) {
            // show optimized cluster points
            viewer.data().add_points(
                refPointLoc,
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

    if (ImGui::CollapsingHeader("Iterative Closest Point", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::PushItemWidth(zebrafishWidth / 3.0);
        ImGui::Checkbox("Show background image", &showBackgroundImage);
        ImGui::Checkbox("Show detected centers", &showMarkerPoints);
        ImGui::Checkbox("Show pattern centers", &showReferencePoints);
        ImGui::SliderInt("Point Size", &pointSize, 1, 30);
        ImGui::PopItemWidth();

        if (ImGui::Button("Load pattern OFF")) {
            std::string filename = FileDialog::openFileName("./.*", {"*.off"});
            if (!filename.empty()) {
                patternFilename = filename;
                Eigen::MatrixXd tempF;
                if (igl::readOFF(patternFilename, refV, tempF)) {

                    PreprocessPatternLoc();
                } else {
                    logger().error("Error open OFF file {}", patternFilename);
                    std::cerr << "Error open OFF file" << std::endl;
                }
            } 
        }
        ImGui::SameLine();
        ImGui::Text("%s", patternFilename.c_str());

        if (ImGui::Button("Run ICP")) {
            
            SearchICP();
            UpdateRefPointLoc();
            // clear the background image
            showBackgroundImage = false;
        }
    }

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::Text("Stage 6: ICP");
}


////////////////////////////////////////////////////////////////////////////////////////
// ICP

void GUI::SearchICP() {

    double RMSerror;
    Eigen::MatrixXd markerLoc;
    markerLoc = markerRecord.loc.block(0, 0, markerRecord.num, 3);

    RMSerror = ICP::RunICP(markerLoc.transpose(), refV.transpose(), Rmat, Tmat);
    logger().info("RunICP error {}", RMSerror);
}


void GUI::PreprocessPatternLoc() {

    refV.col(0).array() -= refV.col(0).minCoeff();  // x
    refV.col(1).array() -= refV.col(1).minCoeff();  // y

    // FIXME: should we scale here?
}


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


void GUI::UpdateRefPointLoc() {

    refPointLoc = (Rmat.transpose() * (refV.transpose().colwise() - Tmat)).transpose();
    // refPointLoc = refV;
    // markerPointLoc = ((Rmat * markerPointLoc.transpose()).colwise() + Tmat).transpose();

    std::cout << Rmat << std::endl;
    std::cout << Tmat << std::endl;

    // std::cout << refPointLoc << std::endl;
}

}  // namespace zebrafish
