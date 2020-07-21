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

    ImGui::Separator(); /////////////////////////////////////////

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

    if (ImGui::CollapsingHeader("Cluster Filter", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::Text("TBD");
    }

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::Text("Stage 5: Filter & Cluster");
}


////////////////////////////////////////////////////////////////////////////////////////
// 

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
    double minValue = radiusCol.minCoeff() - 0.001;
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

}  // namespace zebrafish
