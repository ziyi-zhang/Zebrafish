#include <zebrafish/GUI.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/Quantile.h>


namespace zebrafish {

void NormalizeImage(image_t &image) {

    // normalize all layers
    double quantile = zebrafish::QuantileImage(image, 0.995);
    logger().info("Quantile of image with q=0.995 is {}", quantile);
    for (auto it=image.begin(); it!=image.end(); it++) {
        Eigen::MatrixXd &img = *it;
        img.array() /= quantile;
    }
    logger().info("Image normalized: most pixels will have value between 0 and 1");
}

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 2: Pre-process & B-spline

void GUI::DrawStage2() {

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::Text("Preprocess");
    ImGui::Text("Histogram of pixel brightness");
    ImGui::PlotHistogram("", imgHist.data(), imgHist.size(), 0, NULL, 0, imgHist.maxCoeff(), ImVec2(0, 80));
    ImGui::PushItemWidth(zebrafishWidth / 3.0);
    ImGui::InputDouble("Quantile Threshold", &normalizeQuantile);
    ImGui::PopItemWidth();

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::Text("B-spline Config");
    if (ImGui::Button("Compute Bspline")) {

        logger().info("Computing Bspine for the first frame");
        const int bsplineDegree = 2;
        bsplineSolver.SetResolution(resolutionX, resolutionY, resolutionZ);
        bsplineSolver.CalcControlPts(img, 0.7, 0.7, 0.7, bsplineDegree);
    }

    ImGui::Text("Stage 2: Pre-process & B-spline");
}


////////////////////////////////////////////////////////////////////////////////////////
// Hist


void GUI::ComputeImgHist() {

    static float minValue = std::numeric_limits<float>::max();
    static float maxValue = std::numeric_limits<float>::min();
    imgHist = Eigen::MatrixXf::Zero(histBars, 1);

    for (int i=layerBegin; i<=layerEnd; i++) {

        if (img[i].maxCoeff() > maxValue) maxValue = img[i].maxCoeff();
        if (img[i].minCoeff() < minValue) minValue = img[i].minCoeff();
    }
    const float gap = (maxValue - minValue) / float(histBars);

    for (int i=layerBegin; i<=layerEnd; i++) {
        const Eigen::MatrixXd &mat = img[i];
        for (int j=0; j<mat.size(); j++) {
            imgHist( std::floor((mat(j) - minValue)/gap) )++;
        }
    }
}

}  // namespace zebrafish
