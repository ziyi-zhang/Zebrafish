#include <zebrafish/GUI.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/Quantile.h>


namespace zebrafish {

void NormalizeImage(image_t &image, double thres) {

    // normalize & trim all layers
    for (auto it=image.begin(); it!=image.end(); it++) {
        Eigen::MatrixXd &img = *it;
        for (int r=0; r<img.rows(); r++)
            for (int c=0; c<img.cols(); c++) {
                img(r, c) = (img(r, c)>=thres) ? 1.0f : img(r, c)/thres;
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

        ImGui::Text("Histogram of pixel brightness");
        ImGui::PlotHistogram("", imgHist.data(), imgHist.size(), 0, NULL, 0, imgHist.maxCoeff(), ImVec2(0, 80));
        ImGui::PushItemWidth(zebrafishWidth / 3.0);
        if (ImGui::InputDouble("Quantile Thres", &normalizeQuantile)) {
            
            if (normalizeQuantile>1.0) normalizeQuantile = 1.0;
            if (normalizeQuantile<0.95) normalizeQuantile = 0.95;
            double thres = QuantileImage(img, normalizeQuantile);
            image_t img_ = img;
            NormalizeImage(img_, thres);  // Do not modify "img" now
            logger().info("Image normalized with normalizeQuantile =  {:.4f}  thres =  {:.4f}", normalizeQuantile, thres);

            // Note: This has been calculated once in stage 1: "Reload image" button
            ComputeImgHist(img_);
        }
        ImGui::PopItemWidth();
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("B-spline Config"), ImGuiTreeNodeFlags_DefaultOpen) {

        if (ImGui::TreeNode("Advanced")) {
            ImGui::Text("TBD");
            
            ImGui::TreePop();
            ImGui::Separator();
        }

        if (ImGui::Button("Compute Bspline")) {

            logger().info("Computing Bspine for the first frame");
            const int bsplineDegree = 2;
            bsplineSolver.SetResolution(resolutionX, resolutionY, resolutionZ);
            bsplineSolver.CalcControlPts(img, 0.7, 0.7, 0.7, bsplineDegree);
        }
    }

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::Text("Stage 2: Pre-process & B-spline");
}


////////////////////////////////////////////////////////////////////////////////////////
// Hist


void GUI::ComputeImgHist(image_t &img_) {

    float minValue = std::numeric_limits<float>::max();
    float maxValue = std::numeric_limits<float>::min();
    imgHist = Eigen::MatrixXf::Zero(histBars, 1);

    for (int i=layerBegin; i<=layerEnd; i++) {

        if (img_[i].maxCoeff() > maxValue) maxValue = img_[i].maxCoeff();
        if (img_[i].minCoeff() < minValue) minValue = img_[i].minCoeff();
    }
    const float gap = (maxValue - minValue) / float(histBars);

    for (int i=layerBegin; i<=layerEnd; i++) {
        const Eigen::MatrixXd &mat = img_[i];
        for (int j=0; j<mat.size(); j++) {
            imgHist( std::floor((mat(j) - minValue)/gap) )++;
        }
    }
}

}  // namespace zebrafish
