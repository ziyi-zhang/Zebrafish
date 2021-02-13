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
        MarkerRecursiveDepthCorrection(0, depthCorrectionNum, depthCorrectionGap, false, true);
        UpdateMarkerPointLocArray();
        InitializeICPPattern();
        propertyListType = 2;
        rejectActive = false;
        stage5to6Flag = false;

        logger().debug("   <auto> Finalize cluster locations");
    }

    // Visualize marker cluster points (first frame)
    static int pointSize = 7;
    if (showMarkerPoints) {

        viewer.data().point_size = pointSize;

        if (!markerPointLocArray.empty()) {
            // show optimized markers

            viewer.data().add_points(
                markerPointLocArray[0],
                markerPointColor
            );
        }
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
    }

    if (showICPLines) {
        if (!markerPointLocArray.empty() && refPointLoc.rows()>0)
            DrawICPLines();
    }

    // Visualize manual marker drag code
    MarkerDragVisualization();

    DrawReferenceDots();

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Pattern", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::SliderInt("Pattern rows", &ICP.patternRows, 5, ICP.patternRef * 3);
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Generate a standard triangular mesh pattern with specified rows.\nThe reference pattern should be larger than the detected markers.");
        }
        ImGui::SliderInt("Pattern cols", &ICP.patternCols, 5, ICP.patternRef * 3);
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Generate a standard triangular mesh pattern with specified columns.\nThe reference pattern should be larger than the detected markers.");
        }
        if (ImGui::Button("Generate triangular pattern")) {

            GenerateICPPattern();
            PreprocessPatternLoc();
            UpdateRefPointLoc();
            logger().debug("   <button> Generate triangular pattern");
        }
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Generate a standard triangular mesh pattern with specified size.\nThe reference pattern should be larger than the detected markers.");
        }

        if (ImGui::TreeNode("Advanced pattern")) {

            if (ImGui::Button("Load pattern from file")) {
                std::string filename = FileDialog::openFileName("./.*", {"*.off"});
                if (!filename.empty()) {
                    patternFilename = filename;
                    Eigen::MatrixXd tempF;
                    if (igl::readOFF(patternFilename, ICP.refV, tempF)) {

                        PreprocessPatternLoc();
                        UpdateRefPointLoc();
                    } else {
                        logger().error("Error open OFF file {}", patternFilename);
                        std::cerr << "Error open OFF file" << std::endl;
                    }
                }
                logger().debug("   <button> Load pattern from file");
            }
            ImGui::SameLine();
            ImGui::Text("%s", patternFilename.c_str());

            ImGui::TreePop();
            ImGui::Separator();
        }
    }

    if (ImGui::CollapsingHeader("Iterative Closest Point", ImGuiTreeNodeFlags_DefaultOpen)) {

        if (ImGui::TreeNode("Advanced ICP")) {

            const float inputWidth = ImGui::GetWindowWidth() / 2.0;
            ImGui::PushItemWidth(inputWidth);

            if (ImGui::Button("Reset trans")) {
                ResetICPTransformation();
                showICPLines = false;
                showMarkerMesh = false;
            }
            if (showTooltip && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Reset ICP transformation matrix. This resets the state right before ""Run ICP""\nNote: This does not reset manual transformation matrix.");
            }
            ImGui::SameLine();
            if (ImGui::Button("Reset all")) {
                ResetICPTransformation();
                ICP.xDisp = 0.0f;
                ICP.yDisp = 0.0f;
                ICP.angleRot = 0.0f;
                ICP.scale = 1.0f;
                ICP.xStretch = 0.0f;
                ICP.yStretch = 0.0f;
                ICP.xdistort = 0.0f;
                ICP.ydistort = 0.0f;
                UpdateRefPointManualAlignment();
                UpdateRefPointLoc();
                showICPLines = false;
                showMarkerMesh = false;
            }
            if (showTooltip && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Reset both ICP and manual transformation matrix.");
            }

            // align two point clouds manually
            if (ImGui::SliderFloat("Scale", &ICP.scale, 0.7f, 1.5f, "%.4f x")) {
                UpdateRefPointManualAlignment();
                UpdateRefPointLoc();
            }
            if (ImGui::SliderFloat("Rotation", &ICP.angleRot, -igl::PI / 6.0, igl::PI / 6.0, "%.4f rads")) {
                UpdateRefPointManualAlignment();
                UpdateRefPointLoc();
            }
            if (ImGui::SliderFloat("X displacement", &ICP.xDisp, -0.2 * imgCols, 0.2 * imgCols, "%.3f pixels")) {
                UpdateRefPointManualAlignment();
                UpdateRefPointLoc();
            }
            if (ImGui::SliderFloat("Y displacement", &ICP.yDisp, -0.2 * imgRows, 0.2 * imgRows, "%.3f pixels")) {
                UpdateRefPointManualAlignment();
                UpdateRefPointLoc();
            }
            if (ImGui::SliderFloat("X distort", &ICP.xdistort, -0.1, 0.1, "%.3f")) {
                GenerateICPPattern();
                PreprocessPatternLoc();
                UpdateRefPointManualAlignment();
                UpdateRefPointLoc();
            }
            if (ImGui::SliderFloat("Y distort", &ICP.ydistort, -0.1, 0.1, "%.3f")) {
                GenerateICPPattern();
                PreprocessPatternLoc();
                UpdateRefPointManualAlignment();
                UpdateRefPointLoc();
            }
            if (ImGui::SliderFloat("X stretch", &ICP.xStretch, -0.5, 0.5, "%.2f")) {
                UpdateRefPointManualAlignment();
                UpdateRefPointLoc();
            }
            if (ImGui::SliderFloat("Y stretch", &ICP.yStretch, -0.5, 0.5, "%.2f")) {
                UpdateRefPointManualAlignment();
                UpdateRefPointLoc();
            }

            ImGui::PopItemWidth();
            
            ImGui::TreePop();
        }

        if (ImGui::Button("Run ICP")) {
            
            SearchICP();
            UpdateRefPointLoc();
            showICPLines = true;
            showMarkerMesh = true;
            logger().debug("   <button> Run ICP");
        }
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Align the reference pattern with the detected markers.\nAutomatically create a mesh with the markers.");
        }

        if (ImGui::TreeNode("Advanced ICP visualization")) {

            const float inputWidth = ImGui::GetWindowWidth() / 3.0;
            ImGui::PushItemWidth(inputWidth);
            ImGui::Checkbox("Show background image", &showBackgroundImage);
            ImGui::Checkbox("Show cluster centers", &showMarkerPoints);
            ImGui::Checkbox("Show pattern centers", &showReferencePoints);
            ImGui::Checkbox("Show ICP correspondence", &showICPLines);
            ImGui::Checkbox("Show mesh", &showMarkerMesh);
            ImGui::SliderInt("Point Size", &pointSize, 1, 30);
            ImGui::PopItemWidth();

            ImGui::TreePop();
            ImGui::Separator();
        }

        // this loads the code to render the GUI about mouse draging
        ImGui::Separator();  /////////////////////////////////////////
        RenderMarkerDragGUI();
    }
}


////////////////////////////////////////////////////////////////////////////////////////
// ICP

void GUI::InitializeICPPattern() {

    // need a larger reference pattern for ICP to work
    ICP.patternRef = std::ceil(std::sqrt(markerArray[0].num) * 1.3) + 1;
    if (ICP.patternRef < 5) ICP.patternRef = 5;
    ICP.patternRows = ICP.patternRef;
    ICP.patternCols = ICP.patternRef;
}


void GUI::GenerateICPPattern() {

    ICP.refV.resize(ICP.patternRows * ICP.patternCols, 3);
    ICP.refV_RC.resize(ICP.patternRows * ICP.patternCols, 2);

    // Estimate ICP.patternSpacing
    const int N = markerPointLocArray[0].rows();
    if (N > 2) {

        const Eigen::MatrixXd &loc = markerPointLocArray[0];
        double sumDist = 0.0, dist, minDist;
        for (int i=0; i<N; i++) {

            minDist = std::numeric_limits<double>::max();
            for (int j=0; j<N; j++) {

                if (i==j) continue;
                dist = (loc.row(i) - loc.row(j)).norm();
                if (dist < minDist) minDist = dist;
            }
            sumDist += minDist;
        }

        ICP.patternSpacing = sumDist / double(N) * 1.07;
        logger().debug("ICP pattern spacing estimated as {}", ICP.patternSpacing);
    }

    // Generate pattern
    int count = 0;
    double xdistort = ICP.xdistort * ICP.patternSpacing;
    double ydistort = ICP.ydistort * ICP.patternSpacing;
    for (int row=0; row<ICP.patternRows; row++)
        for (int col=0; col<ICP.patternCols; col++) {

            ICP.refV(count, 0) = (ICP.patternSpacing - xdistort * col) * double(col);
            if (row & 1)  // odd row
                ICP.refV(count, 0) += ICP.patternSpacing / 2.0;
            ICP.refV(count, 1) = std::sqrt(3) / 2.0 * (ICP.patternSpacing - ydistort * row) * double(row);
            ICP.refV(count, 2) = 1.0;
            // also store row & col for extrusion purpose
            ICP.refV_RC(count, 0) = row;
            ICP.refV_RC(count, 1) = col;
            count++;
        }
}


void GUI::ResetICPTransformation() {

    ICP.Rmat = Eigen::MatrixXd::Identity(3, 3);
    ICP.Tmat = Eigen::MatrixXd::Zero(3, 1);

    // will not reset manaul alignment config

    UpdateRefPointLoc();
}


void GUI::SearchICP() {

    double RMSerror;
    Eigen::MatrixXd markerLoc;
    markerLoc = markerArray[0].loc.block(0, 0, markerArray[0].num, 3);
    ICP.matchIdx.resize(markerArray[0].num, 1);

    RMSerror = ICP::RunICP(markerPointLocArray[0].transpose(), ICP.refV_aligned.transpose(), ICP.Rmat, ICP.Tmat, ICP.matchIdx);
    UpdateMarkerMesh();
    logger().debug("RunICP RMS error {}", RMSerror);
}


void GUI::UpdateMarkerMesh() {

    int row, col, i, count;
    std::map<int, int> indexMap;
    const int N = ICP.matchIdx.size();
    Eigen::MatrixXi tempMarkerMeshArray(N*2, 3);

    for (i=0; i<N; i++) {
        indexMap.insert({ICP.matchIdx(i), i});  // pattern idx -> marker idx
    }

    count = 0;
    // iterate triangles pointing up
        //    it3
        //  it1  it2
    for (row=0; row<=ICP.patternRows-2; row++)
        for (col=0; col<=ICP.patternCols-2; col++) {

              int refIdx = row * ICP.patternCols + col;
            auto it1 = indexMap.find(refIdx);
            if (it1 == indexMap.end()) continue;
            auto it2 = indexMap.find(refIdx + 1);
            if (it2 == indexMap.end()) continue;
              refIdx += ICP.patternCols;
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
    for (row=1; row<=ICP.patternRows-1; row++)
        for (col=0; col<=ICP.patternCols-2; col++) {

              int refIdx = row * ICP.patternCols + col;
            auto it1 = indexMap.find(refIdx);
            if (it1 == indexMap.end()) continue;
            auto it2 = indexMap.find(refIdx + 1);
            if (it2 == indexMap.end()) continue;
              refIdx -= ICP.patternCols;
              if (row & 1) refIdx++;  // odd row
            auto it3 = indexMap.find(refIdx);
            if (it3 == indexMap.end()) continue;

            // we find a valid triangle
            tempMarkerMeshArray(count, 0) = it1->second;
            tempMarkerMeshArray(count, 1) = it3->second;
            tempMarkerMeshArray(count, 2) = it2->second;  // oriented!
            count++;
    }

    // update mesh array
    markerMeshArray.resize(count, 3);
    markerMeshArray = tempMarkerMeshArray.block(0, 0, count, 3);
    // DEBUG PURPOSE
    // std::cout << "mesh array" << std::endl;
    // std::cout << markerMeshArray << std::endl;
    // std::cout << ICP.matchIdx.transpose() << std::endl;
}


void GUI::PreprocessPatternLoc() {

    assert(ICP.refV.size() > 0);

    ICP.refV.col(0).array() -= ICP.refV.col(0).minCoeff();  // x
    ICP.refV.col(1).array() -= ICP.refV.col(1).minCoeff();  // y
    ICP.refV.col(2) = Eigen::MatrixXd::Ones(ICP.refV.rows(), 1);  // force "z" to be 1

    // save a copy to "ICP.refV_aligned"
    ICP.refV_aligned = ICP.refV;
}


void GUI::UpdateRefPointManualAlignment() {

    assert(ICP.refV.size() > 0);

    static RMat_t stretchMat;
    stretchMat << ICP.xStretch+1, 0.0, 0.0, 
                  0.0, ICP.yStretch+1, 0.0, 
                  0.0, 0.0, 1.0;

    static RMat_t manualAlignTransMat;
    manualAlignTransMat << std::cos(-ICP.angleRot), - std::sin(-ICP.angleRot), ICP.xDisp, 
                         std::sin(-ICP.angleRot),   std::cos(-ICP.angleRot), ICP.yDisp, 
                         0.0,          0.0,            1.0;

    ICP.refV_aligned = ICP.refV.array() * ICP.scale;
    ICP.refV_aligned = ( stretchMat * ICP.refV_aligned.transpose() ).transpose();
    ICP.refV_aligned = ( manualAlignTransMat * ICP.refV_aligned.transpose() ).transpose();
}


void GUI::UpdateRefPointLoc() {

    /// NOTE: ICP.Rmat is a unitary matrix
    ///       ICP.Rmat.transpose() * ICP.Rmat = Identity
    refPointLoc = (ICP.Rmat.transpose() * (ICP.refV_aligned.transpose().colwise() - ICP.Tmat)).transpose();

    // std::cout << ICP.Rmat << std::endl;
    // std::cout << ICP.Tmat << std::endl;
    // std::cout << refPointLoc << std::endl;
}


void GUI::DrawICPLines() {

    static Eigen::MatrixXd refTemp;
    static Eigen::MatrixXd lineColor(1, 3);
    lineColor << 0.77, 0.28, 0.24;
    viewer.data().line_width = 1.0f;

    const int N = markerArray[0].num;
    if (ICP.matchIdx.rows() != N) return;  // not ready yet
    refTemp.resize(N, 3);

    for (int i=0; i<N; i++) {
        refTemp(i, 0) = refPointLoc(ICP.matchIdx(i), 0);
        refTemp(i, 1) = refPointLoc(ICP.matchIdx(i), 1);
        refTemp(i, 2) = refPointLoc(ICP.matchIdx(i), 2);
    }

    viewer.data().add_edges(markerPointLocArray[0], refTemp, lineColor);
}

}  // namespace zebrafish
