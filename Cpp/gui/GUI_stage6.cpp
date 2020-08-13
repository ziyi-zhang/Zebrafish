#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/FileDialog.h>
#include <zebrafish/ICP.h>

#include <string>
#include <map>


namespace zebrafish {

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 6: Iterative Closest Point

void GUI::DrawStage6() {

    if (stage5to6Flag) {

        FinalizeClusterLoc();
        UpdateMarkerPointLocArray();
        InitializeICPPattern();
        propertyListType = 2;
        rejectActive = false;
        stage5to6Flag = false;
    }

    // Visualize marker cluster points (first frame)
    static int pointSize = 10;
    if (showMarkerPoints) {

        viewer.data().point_size = pointSize;

        if (!markerPointLocArray.empty()) {
            // show optimized markers

            viewer.data().add_points(
                markerPointLocArray[0],
                markerPointColor
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

    // Visualize meshes
    DrawMarkerMesh();

    // Visualize reference points
    if (showReferencePoints) {

        viewer.data().point_size = pointSize;
        Eigen::MatrixXd pointColor(1, 3);
        pointColor << 0.0, 0.447, 0.741;

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

    if (showICPLines) {
        if (!markerPointLocArray.empty() && refPointLoc.rows()>0)
            DrawICPLines();
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Iterative Closest Point", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::PushItemWidth(zebrafishWidth / 3.0);
        ImGui::Checkbox("Show background image", &showBackgroundImage);
        ImGui::Checkbox("Show detected centers", &showMarkerPoints);
        ImGui::Checkbox("Show pattern centers", &showReferencePoints);
        ImGui::Checkbox("Show ICP correspondence", &showICPLines);
        ImGui::Checkbox("Show mesh", &showMarkerMesh);
        ImGui::SliderInt("Point Size", &pointSize, 1, 30);
        ImGui::PopItemWidth();

        ImGui::Separator();  /////////////////////////////////////////

        ImGui::SliderInt("Pattern rows", &ICP_patternRows, 5, ICP_patternRef * 3);
        ImGui::SliderInt("Pattern cols", &ICP_patternCols, 5, ICP_patternRef * 3);
        if (ImGui::Button("Generate triangular pattern")) {

            GenerateICPPattern();
            PreprocessPatternLoc();
            UpdateRefPointLoc();
            logger().debug("   <button> Generate triangular pattern");
        }

        if (ImGui::Button("Load pattern OFF")) {
            std::string filename = FileDialog::openFileName("./.*", {"*.off"});
            if (!filename.empty()) {
                patternFilename = filename;
                Eigen::MatrixXd tempF;
                if (igl::readOFF(patternFilename, refV, tempF)) {

                    PreprocessPatternLoc();
                    UpdateRefPointLoc();
                } else {
                    logger().error("Error open OFF file {}", patternFilename);
                    std::cerr << "Error open OFF file" << std::endl;
                }
            }
            logger().debug("   <button> Load pattern OFF");
        }
        ImGui::SameLine();
        ImGui::Text("%s", patternFilename.c_str());

        ImGui::Separator();  /////////////////////////////////////////

        if (ImGui::Button("Reset ICP transformation")) {
            ResetICP();
        }

        // align two point clouds manually
        if (ImGui::SliderFloat("X displacement", &ICP_xDisp, -0.5 * imgCols, imgRows, "%.2f pixels")) {
            UpdateRefPointManualAlignment();
            UpdateRefPointLoc();
        }
        if (ImGui::SliderFloat("Y displacement", &ICP_yDisp, -imgRows, 0.5 * imgRows, "%.2f pixels")) {
            UpdateRefPointManualAlignment();
            UpdateRefPointLoc();
        }
        if (ImGui::SliderFloat("Rotation", &ICP_angleRot, -igl::PI, igl::PI, "%.3f rads")) {
            UpdateRefPointManualAlignment();
            UpdateRefPointLoc();
        }
        if (ImGui::SliderFloat("Scale", &ICP_scale, 0.2f, 4.0f, "%.3f x")) {
            UpdateRefPointManualAlignment();
            UpdateRefPointLoc();
        }

        if (ImGui::Button("Run ICP")) {
            
            SearchICP();
            UpdateRefPointLoc();
            logger().debug("   <button> Run ICP");
        }
    }
}


////////////////////////////////////////////////////////////////////////////////////////
// ICP

void GUI::InitializeICPPattern() {

    ICP_patternRef = std::ceil(std::sqrt(markerArray[0].num) * 1.2);
    if (ICP_patternRef < 5) ICP_patternRef = 5;
    ICP_patternRows = ICP_patternRef;
    ICP_patternCols = ICP_patternRef;
}


void GUI::GenerateICPPattern() {

    refV.resize(ICP_patternRows * ICP_patternCols, 3);
    int count = 0;
    for (int row=0; row<ICP_patternRows; row++)
        for (int col=0; col<ICP_patternCols; col++) {

            refV(count, 0) = ICP_patternSpacing * double(col);
            if (row & 1)  // odd row
                refV(count, 0) += ICP_patternSpacing / 2.0;
            refV(count, 1) = std::sqrt(3) / 2.0 * ICP_patternSpacing * double(row);
            refV(count, 2) = 1.0;
            count++;
        }
}


void GUI::ResetICP() {

    ICP_Rmat = Eigen::MatrixXd::Identity(3, 3);
    ICP_Tmat = Eigen::MatrixXd::Zero(3, 1);

    // will not reset manaul alignment config

    UpdateRefPointLoc();
}


void GUI::SearchICP() {

    double RMSerror;
    Eigen::MatrixXd markerLoc;
    markerLoc = markerArray[0].loc.block(0, 0, markerArray[0].num, 3);
    ICP_matchIdx.resize(markerArray[0].num, 1);

    RMSerror = ICP::RunICP(markerPointLocArray[0].transpose(), refV_aligned.transpose(), ICP_Rmat, ICP_Tmat, ICP_matchIdx);
    UpdateMarkerMesh();
    logger().info("RunICP error {}", RMSerror);
}


void GUI::UpdateMarkerMesh() {

    int row, col, i, count;
    std::map<int, int> indexMap;
    const int N = ICP_matchIdx.size();
    Eigen::MatrixXi tempMarkerMeshArray(N*2, 3);

    for (i=0; i<N; i++) {
        indexMap.insert({ICP_matchIdx(i), i});  // pattern idx -> marker idx
    }

    count = 0;
    // iterate triangles pointing up
        //    it3
        //  it1  it2
    for (row=0; row<=ICP_patternRows-2; row++)
        for (col=0; col<=ICP_patternCols-2; col++) {

              int refIdx = row * ICP_patternCols + col;
            auto it1 = indexMap.find(refIdx);
            if (it1 == indexMap.end()) continue;
            auto it2 = indexMap.find(refIdx + 1);
            if (it2 == indexMap.end()) continue;
              refIdx += ICP_patternCols;
              if (row & 1) refIdx++;  // odd row
            auto it3 = indexMap.find(refIdx);
            if (it3 == indexMap.end()) continue;

            // we find a valid triangle
            tempMarkerMeshArray(count, 0) = it1->second;
            tempMarkerMeshArray(count, 1) = it2->second;
            tempMarkerMeshArray(count, 2) = it3->second;
            count++;
    }
    // iterate triangles pointing down
        //  it1  it2
        //    it3
    for (row=1; row<=ICP_patternRows-1; row++)
        for (col=0; col<=ICP_patternCols-2; col++) {

              int refIdx = row * ICP_patternCols + col;
            auto it1 = indexMap.find(refIdx);
            if (it1 == indexMap.end()) continue;
            auto it2 = indexMap.find(refIdx + 1);
            if (it2 == indexMap.end()) continue;
              refIdx -= ICP_patternCols;
              if (row & 1) refIdx++;  // odd row
            auto it3 = indexMap.find(refIdx);
            if (it3 == indexMap.end()) continue;

            // we find a valid triangle
            tempMarkerMeshArray(count, 0) = it1->second;
            tempMarkerMeshArray(count, 1) = it2->second;
            tempMarkerMeshArray(count, 2) = it3->second;
            count++;
    }

    // update mesh array
    markerMeshArray.resize(count, 3);
    markerMeshArray = tempMarkerMeshArray.block(0, 0, count, 3);
    // DEBUG PURPOSE
    // std::cout << "mesh array" << std::endl;
    // std::cout << markerMeshArray << std::endl;
    // std::cout << ICP_matchIdx.transpose() << std::endl;
}


void GUI::PreprocessPatternLoc() {

    assert(refV.size() > 0);

    refV.col(0).array() -= refV.col(0).minCoeff();  // x
    refV.col(1).array() -= refV.col(1).minCoeff();  // y
    refV.col(2) = Eigen::MatrixXd::Ones(refV.rows(), 1);  // force "z" to be 1

    // FIXME: should we scale here?

    // save a copy to "refV_aligned"
    refV_aligned = refV;
}


void GUI::UpdateRefPointManualAlignment() {

    assert(refV.size() > 0);

    static RMat_t manualAlignTransMat;
    manualAlignTransMat << std::cos(-ICP_angleRot), - std::sin(-ICP_angleRot), ICP_xDisp, 
                         std::sin(-ICP_angleRot),   std::cos(-ICP_angleRot), ICP_yDisp, 
                         0.0,          0.0,            1.0;

    refV_aligned = refV.array() * ICP_scale;
    refV_aligned = ( manualAlignTransMat * refV_aligned.transpose() ).transpose();
}


void GUI::UpdateRefPointLoc() {

    /// NOTE: ICP_Rmat is a unitary matrix
    ///       ICP_Rmat.transpose() * ICP_Rmat = Identity
    refPointLoc = (ICP_Rmat.transpose() * (refV_aligned.transpose().colwise() - ICP_Tmat)).transpose();

    // std::cout << ICP_Rmat << std::endl;
    // std::cout << ICP_Tmat << std::endl;
    // std::cout << refPointLoc << std::endl;
}


void GUI::DrawICPLines() {

    static Eigen::MatrixXd refTemp;
    static Eigen::MatrixXd lineColor(1, 3);
    lineColor << 0.77, 0.28, 0.24;
    viewer.data().line_width = 1.0f;

    const int N = markerArray[0].num;
    if (ICP_matchIdx.rows() != N) return;  // not ready yet
    refTemp.resize(N, 3);

    for (int i=0; i<N; i++) {
        refTemp(i, 0) = refPointLoc(ICP_matchIdx(i), 0);
        refTemp(i, 1) = refPointLoc(ICP_matchIdx(i), 1);
        refTemp(i, 2) = refPointLoc(ICP_matchIdx(i), 2);
    }

    viewer.data().add_edges(markerPointLocArray[0], refTemp, lineColor);
}

}  // namespace zebrafish
