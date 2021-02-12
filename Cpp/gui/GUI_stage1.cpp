#include <zebrafish/GUI.h>
#include <zebrafish/TiffReader.h>
#include <zebrafish/FileDialog.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/Quantile.h>

#include <igl/unproject_onto_mesh.h>
#include <algorithm>
#include <limits>


namespace zebrafish {

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 1: image read

void GUI::DrawStage1() {

    // cropActive
    if (imageCrop.cropActive && imageCrop.downClicked) {
        // crop activated & has been updated (not default value)
        Eigen::MatrixXd lineColor(1, 3);
        lineColor << 0.77, 0.28, 0.24;
        viewer.data().line_width = 2.0f;

        // upper-left corner (x0, y0)
        // lower-right corner (x1, y1)
        float x0 = imageCrop.baseLoc(0);
        float x1 = std::max(x0, imageCrop.currentLoc(0));
        float y0 = imageCrop.baseLoc(1);
        float y1 = std::min(y0, imageCrop.currentLoc(1));
        DrawRect(x0, y0, x1, y1, lineColor);
    }

    // showCropArea
    if (imageCrop.showCropArea && currentLoadedFrames > 0) {
        // show the area specified by current [r0, c0] x [r1, c1]
        Eigen::MatrixXd lineColor(1, 3);
        lineColor << 0.77, 0.28, 0.24;
        viewer.data().line_width = 2.0f;

        // upper-left corner (x0, y0)
        // lower-right corner (x1, y1)
        float x0 = (imageCrop.c0 == -1) ? 0 : imageCrop.c0;
        float x1 = (imageCrop.c1 == -1) ? imgCols : imageCrop.c1;
        float y0 = (imageCrop.r0 == -1) ? imgRows : imgRows - imageCrop.r0;
        float y1 = (imageCrop.r1 == -1) ? 0 : imgRows - imageCrop.r1;
        DrawRect(x0, y0, x1, y1, lineColor);
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Tiff Image Info", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::PushItemWidth(UIsize.zebrafishWidth / 3.0);
        ImGui::InputInt("Layers Per Image", &layerPerImg);
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("The number of layers in one frame.\nFor example there can be ten frames, and each frame is consisted of 40 2D images of size 1024x1024.");
        }
        ImGui::InputInt("Channels Per Slice", &channelPerSlice);
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("The number of channels this image has.\nZebrafish will only use the data from one channel, but this information is required to properly load the image.");
        }
        ImGui::InputInt("Frames", &ttlFrames);
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("The number of frames or timestamps");
        }
        if (ImGui::TreeNode("Advanced channel")) {

            ImGui::SliderInt("Which channel to load?", &channelToLoad, 0, channelPerSlice-1, "channel %d");
            if (showTooltip && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Only the chosen channel will be used");
            }
            ImGui::TreePop();
            ImGui::Separator();
        }
        ImGui::PopItemWidth();
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Depth Crop", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::SliderInt("Slice index start", &layerBegin, 0, layerEnd);
        ImGui::SliderInt("Slice index end", &layerEnd, layerBegin, layerPerImg-1);
    }
    if (showTooltip && ImGui::IsItemHovered()) {
        ImGui::SetTooltip("Crop the 3D images for all frames in depth.\nOnly slices in this closed interval will be used.");
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Area Crop", ImGuiTreeNodeFlags_DefaultOpen)) {

        if (ImGui::Checkbox("Mouse crop", &imageCrop.cropActive)) {
            if (!imageCrop.cropActive) 
                logger().debug("Mouse crop: de-activated.");
            else {
                logger().debug("Mouse crop: activated.");
                imageCrop.showCropArea = true;
            }
        }
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Use mouse to select an area of interest.\nClick this once and move the mouse to the image. Click and hold. Drag the mouse to select an area. Release.\nThe area will be highlighted in a red box.");
        }
        if (ImGui::Button("Reset crop area")) {
            imageCrop.r0 = -1;
            imageCrop.c0 = -1;
            imageCrop.r1 = -1;
            imageCrop.c1 = -1;
            logger().debug("   <button> Reset crop area");
        }
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Reset the area to be the entire 2D image");
        }

        if (ImGui::TreeNode("Advanced area crop")) {

            ImGui::Checkbox("Show cropped area", &imageCrop.showCropArea);
            ImGui::PushItemWidth(80);
            ImGui::InputInt("R0", &imageCrop.r0);
            ImGui::SameLine();
            ImGui::InputInt("C0", &imageCrop.c0);
            ImGui::InputInt("R1", &imageCrop.r1);
            ImGui::SameLine();
            ImGui::InputInt("C1", &imageCrop.c1);
            ImGui::PopItemWidth();
            ImGui::SliderFloat("Line width", &lineWidth, 1, 16);

            ImGui::TreePop();
            ImGui::Separator();
        }
    }

    ImGui::Separator(); /////////////////////////////////////////

    static bool overrideMaskCheck = false;
    if (ImGui::CollapsingHeader("Mask Crop")) {
        
        ImGui::Checkbox("Override mask check", &overrideMaskCheck);
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Override mask value validity check. But the size still has to match the image.");
        }
        ImGui::SliderFloat("Mask Thres", &maskThres, 0.0, 1.0);

        static std::string loadMaskStr = "";
        if (ImGui::Button("Load mask TIFF image")) {

            loadMaskStr = "failed: see log";
            if (layerPerImg == 1) {
                logger().info("Please first load the image before loading the mask");
            } else {
                if (maskPath == "")
                    maskPath = FileDialog::openFileName("./.*", {"*.tif", "*.tiff"});  // debug mode, pass in through '-m'
                if (!maskPath.empty()) {

                    int layerPerImg_mask, channelPerSlice_mask, ttlFrames_mask;
                    GetDescription(maskPath, layerPerImg_mask, channelPerSlice_mask, ttlFrames_mask);

                    if (layerPerImg_mask != layerPerImg) {
                        logger().error("Z-stack size in the mask does not match with the image");
                        std::cerr << "Z-stack size in the mask does not match with the image" << std::endl;
                    } else {
                        if (ReadTifFirstFrame(maskPath, layerPerImg, 1, membraneMask, imageCrop.r0, imageCrop.c0, imageCrop.r1, imageCrop.c1, 0)) {
                            if (imgRows == membraneMask[0].rows() && imgCols == membraneMask[0].cols()) {

                                double minValue = std::numeric_limits<double>::max();
                                double maxValue = std::numeric_limits<double>::min();
                                for (int i=0; i<membraneMask.size(); i++) {
                                    if (membraneMask[i].maxCoeff() > maxValue) maxValue = membraneMask[i].maxCoeff();
                                    if (membraneMask[i].minCoeff() < minValue) minValue = membraneMask[i].minCoeff();
                                }
                                bool ok = false;
                                if (!overrideMaskCheck) {
                                    if (minValue < 0 || maxValue > 255) {
                                        logger().error("mask value not in [0, 255]: min={} max={}", minValue, maxValue);
                                    } else if (maxValue < 128) {
                                        logger().warn("mask max value < 128. Very likely to reject all markers.");
                                    } else {
                                        ok = true;
                                    }
                                }
                                if (overrideMaskCheck || ok) {
                                    logger().info("mask value min={} max={}", minValue, maxValue);
                                    loadMaskStr = "success";
                                    membraneMaskLoad = true;
                                    membraneMaskCylApply = true;
                                    membraneMaskClusterApply = true;
                                    maskMax = maxValue;
                                }
                            } else {
                                logger().error("Row/Col number in each z-slice does not match the image's number");
                                std::cerr << "Row/Col number in each z-slice does not match the image's number" << std::endl;
                            }
                        } else {
                            logger().error("Error open mask tiff image");
                            std::cerr << "Error open mask tiff image" << std::endl;
                        }
                    }
                }
            }
            logger().debug("   <button> Load mask");
        }
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("[Optional] Load a TIFF image as a binary mask of the membrane area.\nThe mask must be of the same dimension as the background image.");
        }
        ImGui::SameLine();
        ImGui::Text("%s", loadMaskStr.c_str());
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Microscope", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::PushItemWidth(UIsize.zebrafishWidth / 3.0);
        ImGui::InputDouble("Resolution X (um)", &resolutionX, 0.0f, 1000.0f, "%.4f");
        ImGui::InputDouble("Resolution Y (um)", &resolutionY, 0.0f, 1000.0f, "%.4f");
        ImGui::InputDouble("Resolution Z (um)", &resolutionZ, 0.0f, 1000.0f, "%.4f");
        ImGui::PopItemWidth();
    }
    if (showTooltip && ImGui::IsItemHovered()) {
        ImGui::SetTooltip("The physical distance of one pixel (in micrometers).\nResolution X is the distance between two rows; resolution Y is the distance between two columns; resolution Z is the vertical distance between two z-slices.");
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::TreeNode("Advanced visualization")) {

        ImGui::PushItemWidth(UIsize.zebrafishWidth / 2.0);
        ImGui::SliderFloat("Enhance contrast", &stage1contrast, 1.0, 4.0, "%.2f");
        ImGui::PopItemWidth();
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Enhance the contrast of the image. This will only take effect in stage one.");
        }
        ImGui::TreePop();
        ImGui::Separator();
    }

    ImGui::Separator(); /////////////////////////////////////////

    static std::string applyStr = "";
    if (ImGui::Button("Reset")) {
        ImageReadReset();
        logger().debug("   <button> Reset");
    }
    if (showTooltip && ImGui::IsItemHovered()) {
        ImGui::SetTooltip("Reset all parameters entered in this stage");
    }
    ImGui::SameLine();
    if (ImGui::Button("Apply")) {
        if (imagePath.empty()) {
            logger().error("Error: ImagePath is empty");
            return;  // invalid operation
                     // only entry is through "File" - "Open"
        }
        logger().info("Re-load image {}", imagePath);

        try {
            // NOTE: we crop the area [r0, c0]x[r1, c1] when loading the image, but keep the entire z-stack despite the z-crop
            if (ReadTifFirstFrame(imagePath, layerPerImg, channelPerSlice, imgData[0], imageCrop.r0, imageCrop.c0, imageCrop.r1, imageCrop.c1, channelToLoad)) {
                imgRows = imgData[0][0].rows();
                imgCols = imgData[0][0].cols();
                currentLoadedFrames = 1;
                desiredFrames = ttlFrames;
                sliceToShow = layerBegin;  // do not start with slice 0 anymore

                // brightness
                // re-compute compressed image texture
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

                stage1Lock = true;  // allowed to proceed to stage 2
                applyStr = "Done";
                logger().info("Image reloaded");
            } else {
                logger().error("Error open tiff image (reload)");
                std::cerr << "Error open tiff image (reload)" << std::endl;
                applyStr = "Failed";
            }
        } catch (const std::exception &e) {
            logger().error("   <button> [Apply] Fatal error when trying to re-load the image. This is often due to wrong metadata like incorrect number of channels or slices.");
            std::cerr << "   <button> [Apply] Fatal error when trying to re-load the image. This is often due to wrong metadata like incorrect number of channels or slices." << std::endl;
            applyStr = "Failed";
        }
        logger().debug("   <button> Apply");
    }
    if (showTooltip && ImGui::IsItemHovered()) {
        ImGui::SetTooltip("Preview the image with the parameters in this page.\nThe result will be the image used by upcoming stages.");
    }
    ImGui::SameLine();
    ImGui::Text("%s", applyStr.c_str());
}


void GUI::ImageReadReset() {

    layerPerImg = 1;
    channelPerSlice = 1;
    channelToLoad = 0;
    layerBegin = 0;
    layerEnd = layerPerImg - 1;
    resolutionX = 0;
    resolutionY = 0;
    resolutionZ = 0;
    imageCrop.r0 = -1;
    imageCrop.c0 = -1;
    imageCrop.r1 = -1;
    imageCrop.c1 = -1;
    imageCrop.showCropArea = true;
}


////////////////////////////////////////////////////////////////////////////////////////
// Crop image


void GUI::CropImage(const Eigen::Vector2f &mouse, MOUSE_TYPE mousetype, crop_t &cropStruct) {

    static int fid;
    static Eigen::Vector3f bc, loc;
    static bool hit;

    hit = igl::unproject_onto_mesh(
        Eigen::Vector2f(mouse(0), viewer.core().viewport(3)-mouse(1)), 
        viewer.core().view, 
        viewer.core().proj,
        viewer.core().viewport, 
        V, 
        F, 
        fid, 
        bc);

    if (hit) {
        // has hit

        loc(0) = (V(F(fid, 0), 0) * bc(0) + V(F(fid, 1), 0) * bc(1) + V(F(fid, 2), 0) * bc(2));
        loc(1) = (V(F(fid, 0), 1) * bc(0) + V(F(fid, 1), 1) * bc(1) + V(F(fid, 2), 1) * bc(2));
        loc(2) = (V(F(fid, 0), 2) * bc(0) + V(F(fid, 1), 2) * bc(1) + V(F(fid, 2), 2) * bc(2));

        if (mousetype == MOUSEDOWN) {

            logger().debug("Mouse down: fid = {} | bc = {:.2f} x {:.2f} x {:.2f} | loc = {:.2f} x {:.2f} x {:.2f}", fid, bc(0), bc(1), bc(2), loc(0), loc(1), loc(2));
            cropStruct.downClicked = true;

            cropStruct.baseLoc = loc;
            cropStruct.currentLoc = loc;
        } else if (mousetype == MOUSEUP) {
            logger().debug("Mouse up:   fid = {} | bc = {:.2f} x {:.2f} x {:.2f} | loc = {:.2f} x {:.2f} x {:.2f}", fid, bc(0), bc(1), bc(2), loc(0), loc(1), loc(2));
            cropStruct.cropActive = false;
            cropStruct.downClicked = false;

            cropStruct.c0 = std::round(cropStruct.baseLoc(0));
            cropStruct.c1 = std::round(cropStruct.currentLoc(0));
            cropStruct.r0 = imgRows - std::round(cropStruct.baseLoc(1));
            cropStruct.r1 = imgRows - std::round(cropStruct.currentLoc(1));

            // [r0, c0] must be upper-left, [r1, c1] must be lower-right
            if (cropStruct.r0 > cropStruct.r1 && cropStruct.c0 > cropStruct.c1) {
                // cropping from lower-right to upper-left, need to swap
                int t;
                t = cropStruct.c0; cropStruct.c0 = cropStruct.c1; cropStruct.c1 = t;
                t = cropStruct.r0; cropStruct.r0 = cropStruct.r1; cropStruct.r1 = t;
            }
            if (!(cropStruct.r1 > cropStruct.r0 && cropStruct.c1 > cropStruct.c0)) {
                // invalid
                cropStruct.r0 = -1;
                cropStruct.c0 = -1;
                cropStruct.r1 = -1;
                cropStruct.c1 = -1;
                logger().info("Mouse up invalid: has been reset to -1");
            }
        } else if (mousetype == MOUSEMOVE) {
            cropStruct.currentLoc = loc;
        }
    } else {
        // no hit

        if (mousetype == MOUSEDOWN) {
            logger().debug("Mouse down: no hit");
            cropStruct.downClicked = true;
        } else if (mousetype == MOUSEUP) {
            logger().debug("Mouse up:   no hit");
            cropStruct.cropActive = false;
            cropStruct.downClicked = false;
        }
    }
}


void GUI::DrawRect(double x0, double y0, double x1, double y1, const Eigen::MatrixXd &lineColor) {

    static Eigen::MatrixXd lineX, lineY;
    lineX.resize(4, 3);
    lineY.resize(4, 3);
    // upper-left corner (x0, y0)
    // lower-right corner (x1, y1)
    lineX << x0, y0, 0.1, 
             x0, y1, 0.1, 
             x1, y1, 0.1, 
             x1, y0, 0.1;
    lineY << x0, y1, 0.1, 
             x1, y1, 0.1, 
             x1, y0, 0.1, 
             x0, y0, 0.1;
    viewer.data().line_width = lineWidth;
    viewer.data().add_edges(lineX, lineY, lineColor);
}

}  // namespace zebrafish
