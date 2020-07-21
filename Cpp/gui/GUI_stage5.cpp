#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>


namespace zebrafish {

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 5: Filter & Cluster

void GUI::DrawStage5() {

    if (stage4to5Flag) {

        UpdateCylEnergyHist();
        UpdateCylRadiusHist();
        UpdateCylIterHist();
        stage4to5Flag = false;
    }

    // Visualize filtered points
    static int filterPointSize = 7;
    if (showCylFilterPoints) {

        static float cylinderEnergyThres_cache = -1;
        static float cylinderRadiusThres_cache = -1;
        static int   cylinderIterThres_cache = -1;
        if (cylinderEnergyThres != cylinderEnergyThres_cache ||
            cylinderRadiusThres != cylinderRadiusThres_cache ||
            cylinderIterThres   != cylinderIterThres_cache) {
            // Update the points to visualize when the three thresholds have changed

            CylinderFilter();
            UpdateCylPointLoc();

            cylinderEnergyThres_cache = cylinderEnergyThres;
            cylinderRadiusThres_cache = cylinderRadiusThres;
            cylinderIterThres_cache   = cylinderIterThres;
        }

        viewer.data().point_size = filterPointSize;
        Eigen::MatrixXd pointColor(1, 3);
        pointColor << 0.87, 0.33, 0.33;

        if (cylPointLoc.rows() > 0) {
            // show optimized points
            viewer.data().add_points(
                cylPointLoc,
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

    ImGui::Checkbox("Show filtered locations", &showCylFilterPoints);

    ImGui::PushItemWidth(zebrafishWidth / 3.0);
    ImGui::SliderInt("Point Size", &filterPointSize, 1, 30);
    ImGui::PopItemWidth();

    ImGui::Separator(); /////////////////////////////////////////
    
    // ----------------------------------------------------------

    if (ImGui::CollapsingHeader("Cylinder Filter", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImVec2 before, after;
        float ratio;
        ImDrawList *drawList = ImGui::GetWindowDrawList();
        // Histogram of cylinder energy
        ImGui::Text("Histogram of energy");

        const float width = ImGui::GetWindowWidth() * 0.75f - 2;
        ImGui::PushItemWidth(width + 2);

        before = ImGui::GetCursorScreenPos();
        ImGui::PlotHistogram("", cylEnergyHist.hist.data(), cylEnergyHist.hist.size(), 0, NULL, 0, cylEnergyHist.hist.maxCoeff(), ImVec2(0, 80));
        after = ImGui::GetCursorScreenPos();
        after.y -= ImGui::GetStyle().ItemSpacing.y;

        ratio = ((cylinderEnergyThres-cylEnergyHist.minValue)/(cylEnergyHist.maxValue-cylEnergyHist.minValue));
        ratio = std::min(1.0f, std::max(0.0f, ratio));
        drawList->PushClipRectFullScreen();
        drawList->AddLine(
            ImVec2(before.x + width * ratio, before.y), 
            ImVec2(before.x + width * ratio, after.y), 
            IM_COL32(50, 205, 50, 255), 
            2.0f
        );
        drawList->PopClipRect();
        ImGui::PopItemWidth();
        ImGui::SliderFloat("energy threshold", &cylinderEnergyThres, cylEnergyHist.minValue, cylEnergyHist.maxValue);

        ImGui::Separator(); /////////////////////////////////////////

        // Histogram of cylinder radius
        ImGui::Text("Histogram of cylinder radius");

        ImGui::PushItemWidth(width + 2);

        before = ImGui::GetCursorScreenPos();
        ImGui::PlotHistogram("", cylRadiusHist.hist.data(), cylRadiusHist.hist.size(), 0, NULL, 0, cylRadiusHist.hist.maxCoeff(), ImVec2(0, 80));
        after = ImGui::GetCursorScreenPos();
        after.y -= ImGui::GetStyle().ItemSpacing.y;

        ratio = ((cylinderRadiusThres-cylRadiusHist.minValue)/(cylRadiusHist.maxValue-cylRadiusHist.minValue));
        ratio = std::min(1.0f, std::max(0.0f, ratio));
        drawList->PushClipRectFullScreen();
        drawList->AddLine(
            ImVec2(before.x + width * ratio, before.y), 
            ImVec2(before.x + width * ratio, after.y), 
            IM_COL32(50, 205, 50, 255), 
            2.0f
        );
        drawList->PopClipRect();
        ImGui::PopItemWidth();
        ImGui::SliderFloat("maximum radius", &cylinderRadiusThres, cylRadiusHist.minValue, cylRadiusHist.maxValue, "%.2f pixels");
        
        ImGui::Separator(); /////////////////////////////////////////

        // Histogram of cylinder radius
        ImGui::Text("Histogram of optimization iterations");

        ImGui::PushItemWidth(width + 2);

        before = ImGui::GetCursorScreenPos();
        ImGui::PlotHistogram("", cylIterHist.hist.data(), cylIterHist.hist.size(), 0, NULL, 0, cylIterHist.hist.maxCoeff(), ImVec2(0, 80));
        after = ImGui::GetCursorScreenPos();
        after.y -= ImGui::GetStyle().ItemSpacing.y;

        ratio = ((cylinderIterThres-cylIterHist.minValue)/(cylIterHist.maxValue-cylIterHist.minValue));
        ratio = std::min(1.0f, std::max(0.0f, ratio));
        drawList->PushClipRectFullScreen();
        drawList->AddLine(
            ImVec2(before.x + width * ratio, before.y), 
            ImVec2(before.x + width * ratio, after.y), 
            IM_COL32(50, 205, 50, 255), 
            2.0f
        );
        drawList->PopClipRect();
        ImGui::PopItemWidth();
        ImGui::SliderInt("maximum iteration", &cylinderIterThres, 1, optimMaxIt, "%d iterations");
    }

    ImGui::Separator(); /////////////////////////////////////////

    // ----------------------------------------------------------

    if (ImGui::CollapsingHeader("Cluster Filter", ImGuiTreeNodeFlags_DefaultOpen)) {

        if (ImGui::Button("Cluster")) {
            
            Cluster();
        }
        ImGui::Text("TBD");
    }

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::Text("Stage 5: Filter & Cluster");
}


////////////////////////////////////////////////////////////////////////////////////////
// Cylidner Filter

void GUI::CylinderFilter() {

    const int N = pointRecord.num;

    for (int i=0; i<N; i++) {

        pointRecord.alive(i) = 
            (pointRecord.optimization(i, 4) < cylinderEnergyThres) &&
            (pointRecord.optimization(i, 3) < cylinderRadiusThres) &&
            (pointRecord.optimization(i, 5) < cylinderIterThres);
    }
}


void GUI::UpdateCylPointLoc() {

    const int N = pointRecord.num;
    int M = 0, i, count;
    Eigen::MatrixXd tempLoc;

    for (i=0; i<N; i++)
        if (pointRecord.alive(i)) M++;
    if (M == 0) return;

    cylPointLoc.resize(M, 3);
    tempLoc.resize(M, 3);

    count = 0;
    for (i=0; i<N; i++)
        if (pointRecord.alive(i)) {
            tempLoc(count, 0) = pointRecord.optimization(i, 0);
            tempLoc(count, 1) = pointRecord.optimization(i, 1);
            tempLoc(count, 2) = pointRecord.optimization(i, 2);
            count++;
        }

    assert(count == M);

    cylPointLoc.col(0) = tempLoc.col(1).array() + 0.5;
    cylPointLoc.col(1) = (imgRows-0.5) - tempLoc.col(0).array();
    cylPointLoc.col(2) = tempLoc.col(2);

    logger().info("[Visualization] Filtered points updated: total number = {}", M);
}


void GUI::UpdateCylEnergyHist() {
// Note: only count cylinders with energy smaller than 0
//       some cylinders failed to converge and thus have energy 1.0

    Eigen::MatrixXd energyCol = pointRecord.optimization.col(4);
    double minValue = energyCol.minCoeff();
    double maxValue = 0.0;
    const int N = pointRecord.num;
    assert(N > 0);

    const double epsilon = 0.0001;  // to make sure every number lies inside
    minValue -= epsilon;
    const double gap = (maxValue - minValue) / double(histBars);

    cylEnergyHist.hist = Eigen::MatrixXf::Zero(histBars, 1);
    cylEnergyHist.minValue = minValue;
    cylEnergyHist.maxValue = maxValue;

    int idx;
    for (int i=0; i<N; i++) {
        idx = std::floor((energyCol(i) - minValue)/gap);
        if (idx >= histBars) continue;  // energy > 0.0
        cylEnergyHist.hist(idx)++;
    }
}


void GUI::UpdateCylRadiusHist() {

    Eigen::MatrixXd radiusCol = pointRecord.optimization.col(3);
    double minValue = 0.0;
    double maxValue = 12.0;
    const int N = pointRecord.num;
    assert(N > 0);

    const double gap = (maxValue - minValue) / double(histBars);

    cylRadiusHist.hist = Eigen::MatrixXf::Zero(histBars, 1);
    cylRadiusHist.minValue = minValue;
    cylRadiusHist.maxValue = maxValue;

    int idx;
    for (int i=0; i<N; i++) {
        idx = std::floor((radiusCol(i) - minValue)/gap);
        if (idx >= histBars) idx = histBars - 1;
        if (idx < 0) idx = 0;
        cylRadiusHist.hist(idx)++;
    }
}


void GUI::UpdateCylIterHist() {

    Eigen::MatrixXd iterCol = pointRecord.optimization.col(5);
    double minValue = iterCol.minCoeff();
    double maxValue = iterCol.maxCoeff();
    const int N = pointRecord.num;
    assert(N > 0);

    const double epsilon = 0.0001;  // to make sure every number lies inside
    maxValue += epsilon;
    minValue -= epsilon;
    const double gap = (maxValue - minValue) / double(histBars);

    cylIterHist.hist = Eigen::MatrixXf::Zero(histBars, 1);
    cylIterHist.minValue = minValue;
    cylIterHist.maxValue = maxValue;

    for (int i=0; i<N; i++) {
        cylIterHist.hist( std::floor((iterCol(i) - minValue)/gap) )++;
    }
}

// ---------------------------------------------------------------

void GUI::Cluster() {


}

}  // namespace zebrafish
