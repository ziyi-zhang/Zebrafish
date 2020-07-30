#include <zebrafish/GUI.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/Quantile.h>


namespace zebrafish {

void NormalizeImage(image_t &image, double thres) {

    // normalize & trim all layers
    for (auto it=image.begin(); it!=image.end(); it++) {
        Eigen::MatrixXd &slice = *it;
        for (int r=0; r<slice.rows(); r++)
            for (int c=0; c<slice.cols(); c++) {
                slice(r, c) = (slice(r, c)>=thres) ? 1.0f : slice(r, c)/thres;
            }
    }
}

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 2: Pre-process & B-spline

void GUI::DrawStage2() {

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Preprocess", ImGuiTreeNodeFlags_DefaultOpen)) {

        // Histogram
        ImGui::Text("Histogram of pixel brightness");

        const float width = ImGui::GetWindowWidth() * 0.75f - 2;
        ImGui::PushItemWidth(width + 2);

        ImVec2 before = ImGui::GetCursorScreenPos();
        ImGui::PlotHistogram("", imgHist.hist.data(), imgHist.hist.size(), 0, NULL, 0, imgHist.hist.maxCoeff(), ImVec2(0, 80));
        ImVec2 after = ImGui::GetCursorScreenPos();
        after.y -= ImGui::GetStyle().ItemSpacing.y;

            // logger().debug("normalizeQuantile = {}   normalizeQuantileRes = {}", normalizeQuantile, normalizeQuantileRes);
        float ratio = ((normalizeQuantileRes-imgHist.minValue)/(imgHist.maxValue-imgHist.minValue));
        ratio = std::min(1.0f, std::max(0.0f, ratio));
        ImDrawList *drawList = ImGui::GetWindowDrawList();
        drawList->PushClipRectFullScreen();
        drawList->AddLine(
            ImVec2(before.x + width * ratio, before.y), 
            ImVec2(before.x + width * ratio, after.y), 
            IM_COL32(50, 205, 50, 255), 
            2.0f
        );
        drawList->PopClipRect();
        ImGui::PopItemWidth();

        ImGui::PushItemWidth(zebrafishWidth * 0.75);
        ImGui::Text("Quantile Thres");
        if (ImGui::SliderFloat("", &normalizeQuantile, 0.95, 0.999, "%.3f") || stage1to2Flag) {

            if (stage1to2Flag) ComputeImgHist(imgData[0]);
            stage1to2Flag = false;

            if (normalizeQuantile>0.999) normalizeQuantile = 0.999;
            if (normalizeQuantile<0.95) normalizeQuantile = 0.95;
            normalizeQuantileRes = QuantileImage(imgData[0], normalizeQuantile);
            image_t img_ = imgData[0];
            NormalizeImage(img_, normalizeQuantileRes);  // Do not modify "imgData[0]" now
            logger().info("Trial: normalizeQuantile =  {:.4f}  thres =  {:.4f}", normalizeQuantile, normalizeQuantileRes);

            ComputeCompressedImg(img_, 0);
        }
        ImGui::PopItemWidth();
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("B-spline Config", ImGuiTreeNodeFlags_DefaultOpen)) {

        if (ImGui::TreeNode("Advanced B-spline")) {
            ImGui::Text("Not implemented yet...");
            
            ImGui::TreePop();
            ImGui::Separator();
        }

        if (ImGui::Button("Compute B-spline")) {

            // IMPORTANT: The effect of quantile does not affect the raw data until now
            // put into effect
            normalizeQuantileRes = QuantileImage(imgData[0], normalizeQuantile);
            NormalizeImage(imgData[0], normalizeQuantileRes);
            ComputeCompressedImg(imgData[0], 0);
            logger().info("[Finalized] Image normalized with normalizeQuantile =  {:.4f}  thres =  {:.4f}", normalizeQuantile, normalizeQuantileRes);

            // Compute B-spline
            logger().info("Computing Bspine for the first frame");
            const int bsplineDegree = 2;
            bsplineSolver.SetResolution(resolutionX, resolutionY, resolutionZ);
            bsplineSolver.CalcControlPts(imgData[0], 0.7, 0.7, 0.7, bsplineDegree);
        }
    }

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::Text("Stage 2: Pre-process & B-spline");
}


////////////////////////////////////////////////////////////////////////////////////////
// Hist


void GUI::ComputeImgHist(const image_t &img_) {

    double minValue = std::numeric_limits<double>::max();
    double maxValue = std::numeric_limits<double>::min();

    for (int i=layerBegin; i<=layerEnd; i++) {

        if (img_[i].maxCoeff() > maxValue) maxValue = img_[i].maxCoeff();
        if (img_[i].minCoeff() < minValue) minValue = img_[i].minCoeff();
    }
    const double epsilon = 0.0001;  // to make sure every number lies inside
    maxValue += epsilon;
    minValue -= epsilon;
    const double gap = (maxValue - minValue) / double(histBars);
    imgHist.hist = Eigen::MatrixXf::Zero(histBars, 1);
    imgHist.minValue = minValue;
    imgHist.maxValue = maxValue;
    assert(gap > 0);

    for (int i=layerBegin; i<=layerEnd; i++) {
        const Eigen::MatrixXd &mat = img_[i];
        for (int j=0; j<mat.size(); j++) {
            imgHist.hist( std::floor((mat(j) - minValue)/gap) )++;
        }
    }
}

}  // namespace zebrafish
