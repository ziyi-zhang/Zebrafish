#include <zebrafish/GUI.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/Quantile.h>
#include <zebrafish/TiffReader.h>


namespace zebrafish {

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 2: Pre-process & B-spline

void GUI::DrawStage2() {

    static bool updateNormalizedTexture;
    static std::string bsplineStr = "";
    if (stage1to2Flag) {

        imageCrop.showCropArea = false;
        ComputeImgHist(imgData[0]);
        updateNormalizedTexture = true;
        bsplineStr = "";
        stage1to2Flag = false;
    }

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

        ImGui::PushItemWidth(UIsize.zebrafishWidth * 0.75);
        ImGui::Text("Outlier brightness threshold (in quantile)");
        if (ImGui::SliderFloat("", &normalizeQuantile, 0.95, 0.999, "%.3f") || updateNormalizedTexture) {

            updateNormalizedTexture = false;
            normalizeQuantileRes = QuantileImage(imgData[0], normalizeQuantile, layerBegin, layerEnd);
            image_t img_ = imgData[0];
            NormalizeImage(img_, normalizeQuantileRes);  // Do not modify "imgData[0]" now
            logger().info("Slider trial: normalizeQuantile =  {:.4f}  thres =  {:.4f}", normalizeQuantile, normalizeQuantileRes);

            switch (imageViewerCompressType) {
                case COMPRESS_AVG:
                    ComputeCompressedTextureAvg(img_, 0);
                    break;
                case COMPRESS_MAX:
                    ComputeCompressedTextureMax(img_, 0);
                    break;
                default:
                    assert(false);
                    break;
            }
        }
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Set the pixel brightness threshold. All pixels brighter than this value will be trimmed.\nAn appropriate value would move the green line at the tail of the histogram, leaving only negligible bars to its right.");
        }
        ImGui::PopItemWidth();
    }

    ImGui::Separator(); /////////////////////////////////////////

    static float bsp_xratio = 0.7, bsp_yratio = 0.7, bsp_zratio = 0.7;
    if (ImGui::CollapsingHeader("B-spline Config", ImGuiTreeNodeFlags_DefaultOpen)) {

        if (ImGui::TreeNode("Advanced B-spline")) {
            
            if (ImGui::SliderInt("B-spline degree", &bsplineDegree, 2, 3, "Degree %d")) {
                bsplineArray[0].Set_degree(bsplineDegree);
            }
            if (ImGui::InputDouble("B-spline solver tolerance", &bsplineSolverTol, 0.0, 0.0, "%.2E")) {
                bsplineArray[0].Set_solverTol(bsplineSolverTol);
            }
            ImGui::SliderFloat("X ratio", &bsp_xratio, 0.2, 1.0, "%.2f");
            ImGui::SliderFloat("Y ratio", &bsp_yratio, 0.2, 1.0, "%.2f");
            ImGui::SliderFloat("Z ratio", &bsp_zratio, 0.2, 1.0, "%.2f");

            ImGui::TreePop();
            ImGui::Separator();
        }

        if (ImGui::Button("Compute B-spline")) {

            // IMPORTANT: The effect of quantile does not affect the raw data until now
            // put into effect
            normalizeQuantileRes = QuantileImage(imgData[0], normalizeQuantile, layerBegin, layerEnd);
            NormalizeImage(imgData[0], normalizeQuantileRes);
            switch (imageViewerCompressType) {
                case COMPRESS_AVG:
                    ComputeCompressedTextureAvg(imgData[0], 0);
                    break;
                case COMPRESS_MAX:
                    ComputeCompressedTextureMax(imgData[0], 0);
                    break;
                default:
                    assert(false);
                    break;
            }
            logger().info("[Finalized] Image normalized with normalizeQuantile =  {:.4f}  thres =  {:.4f}", normalizeQuantile, normalizeQuantileRes);

            // Compute B-spline
            logger().info("Computing Bspine for the first frame");
            const int bsplineDegree = 2;
            bsplineArray[0].SetResolution(resolutionX, resolutionY, resolutionZ);
            bsplineArray[0].CalcControlPts(imgData[0], bsp_xratio, bsp_yratio, bsp_zratio, bsplineDegree);
            bsplineStr = "Done";
            stage2Lock = true;
            logger().debug("   <button> Compute B-spline");
        }
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Compute B-spline control points for the first frame. This may take some time.\nNote: If the brightness threshold changes, the B-spline must be re-computed.");
        }
        ImGui::SameLine();
        ImGui::Text("%s", bsplineStr.c_str());
    }

    /*
    if (ImGui::Button("Interp Img")) {
        std::vector<Eigen::MatrixXd> img;
        Eigen::VectorXd xArray = Eigen::VectorXd::LinSpaced(imgRows-4, 2, imgRows-3);
        Eigen::VectorXd yArray = Eigen::VectorXd::LinSpaced(imgCols-4, 2, imgCols-3);
        // Eigen::VectorXd zArray = Eigen::VectorXd::LinSpaced(layerPerImg-4, layerBegin+2, layerEnd-2);
        Eigen::VectorXd zArray = Eigen::VectorXd::LinSpaced(1, 4, 4);

        const auto InterpImage = [&xArray, &yArray, &zArray, this, &img]() {
            img.clear();
            for (int iz=0; iz<zArray.size(); iz++) {
                Eigen::MatrixXd slice(xArray.size(), yArray.size());
                for (int ix=0; ix<xArray.size(); ix++)
                    for (int iy=0; iy<yArray.size(); iy++) {
                        slice(ix, iy) = bsplineArray[0].Interp3D(xArray(ix), yArray(iy), zArray(iz));
                    }
                img.push_back(slice);
            }
        };

        InterpImage();
        std::string outName = "/Users/ziyizhang/Projects/tmp/interp_zebra-" + std::to_string(bsp_xratio) + "-" + std::to_string(bsp_yratio) + ".tif";
        WriteTif(outName, img, 0, img.size()-1);
        logger().info("Interp image saved to {}", outName);
    }
    */

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Load All Frames", ImGuiTreeNodeFlags_DefaultOpen)) {

        const float inputWidth = ImGui::GetWindowWidth() / 3.0;
        ImGui::PushItemWidth(inputWidth);

        ImGui::SliderInt("Desired #frames", &desiredFrames, 1, ttlFrames, "%d frames");
        if (ImGui::Button("Prepare all frames")) {
            LoadSubsequentFrames();
            ComputeBsplineForAllFrames();
            // do not really change "currentLoadedFrames" now. Will crash something
            currentLoadedFrames_temp = currentLoadedFrames;
            currentLoadedFrames = 1;

            preLoadAllFrames = true;
            logger().debug("   <button> Prepare all frames");
        }
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("[Optional] Load desired number of frames and prepare them for optimization. This may take some time.\nThis is optional as it can be postponed to stage-7.");
        }

        ImGui::PopItemWidth();
    }
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
