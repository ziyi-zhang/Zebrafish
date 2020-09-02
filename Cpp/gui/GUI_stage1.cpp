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
    if (cropActive && downClicked) {
        // crop activated & has been updated (not default value)
        Eigen::MatrixXd lineColor(1, 3);
        lineColor << 0.77, 0.28, 0.24;
        viewer.data().line_width = 1.0f;

        // upper-left corner (x0, y0)
        // lower-right corner (x1, y1)
        float x0 = baseLoc(0);
        float x1 = std::max(x0, currentLoc(0));
        float y0 = baseLoc(1);
        float y1 = std::min(y0, currentLoc(1));
        DrawRect(x0, y0, x1, y1, lineColor);
    }

    // showCropArea
    if (showCropArea && currentLoadedFrames > 0) {
        // show the area specified by current [r0, c0] x [r1, c1]
        Eigen::MatrixXd lineColor(1, 3);
        lineColor << 0.77, 0.28, 0.24;
        viewer.data().line_width = 1.0f;

        // upper-left corner (x0, y0)
        // lower-right corner (x1, y1)
        float x0 = (c0 == -1) ? 0 : c0;
        float x1 = (c1 == -1) ? imgCols : c1;
        float y0 = (r0 == -1) ? imgRows : imgRows - r0;
        float y1 = (r1 == -1) ? 0 : imgRows - r1;
        DrawRect(x0, y0, x1, y1, lineColor);
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Tiff Image Info", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::PushItemWidth(zebrafishWidth / 3.0);
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

        if (ImGui::Checkbox("Mouse crop", &cropActive)) {
            if (!cropActive) 
                logger().debug("Mouse crop: de-activated.");
            else {
                logger().debug("Mouse crop: activated.");
                showCropArea = true;
            }
        }
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Use mouse to select an area of interest.\nClick this once and move the mouse to the image. Click and hold. Drag the mouse to select an area. Release.\nThe area will be highlighted in a red box.");
        }
        if (ImGui::Button("Reset crop area")) {
            r0 = -1;
            c0 = -1;
            r1 = -1;
            c1 = -1;
            logger().debug("   <button> Reset crop area");
        }
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Reset the area to be the entire 2D image");
        }

        if (ImGui::TreeNode("Advanced area crop")) {

            ImGui::Checkbox("Show cropped area", &showCropArea);
            ImGui::PushItemWidth(80);
            ImGui::InputInt("R0", &r0);
            ImGui::SameLine();
            ImGui::InputInt("C0", &c0);
            ImGui::InputInt("R1", &r1);
            ImGui::SameLine();
            ImGui::InputInt("C1", &c1);
            ImGui::PopItemWidth();
            ImGui::SliderFloat("Line width", &lineWidth, 1, 16);

            ImGui::TreePop();
            ImGui::Separator();
        }
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Microscope", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::PushItemWidth(zebrafishWidth / 3.0);
        ImGui::InputDouble("Resolution X (um)", &resolutionX);
        ImGui::InputDouble("Resolution Y (um)", &resolutionY);
        ImGui::InputDouble("Resolution Z (um)", &resolutionZ);
        ImGui::PopItemWidth();
    }
    if (showTooltip && ImGui::IsItemHovered()) {
        ImGui::SetTooltip("The physical distance of one pixel (in micrometers).\nResolution X is the distance between two rows; resolution Y is the distance between two columns; resolution Z is the vertical distance between two z-slices.");
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::TreeNode("Advanced visualization")) {

        ImGui::PushItemWidth(zebrafishWidth / 2.0);
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
            if (ReadTifFirstFrame(imagePath, layerPerImg, channelPerSlice, imgData[0], r0, c0, r1, c1, channelToLoad)) {
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
    r0 = -1;
    c0 = -1;
    r1 = -1;
    c1 = -1;
    showCropArea = true;
}


////////////////////////////////////////////////////////////////////////////////////////
// Crop image


void GUI::CropImage(const Eigen::Vector2f &mouse, MOUSE_TYPE mousetype) {

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
            downClicked = true;

            baseLoc = loc;
            currentLoc = loc;
        } else if (mousetype == MOUSEUP) {
            logger().debug("Mouse up:   fid = {} | bc = {:.2f} x {:.2f} x {:.2f} | loc = {:.2f} x {:.2f} x {:.2f}", fid, bc(0), bc(1), bc(2), loc(0), loc(1), loc(2));
            cropActive = false;
            downClicked = false;

            c0 = std::round(baseLoc(0));
            c1 = std::round(currentLoc(0));
            r0 = imgRows - std::round(baseLoc(1));
            r1 = imgRows - std::round(currentLoc(1));

            // [r0, c0] must be upper-left, [r1, c1] must be lower-right
            if (r0 > r1 && c0 > c1) {
                // cropping from lower-right to upper-left, need to swap
                int t;
                t = c0; c0 = c1; c1 = t;
                t = r0; r0 = r1; r1 = t;
            }
            if (!(r1 > r0 && c1 > c0)) {
                // invalid
                r0 = -1;
                c0 = -1;
                r1 = -1;
                c1 = -1;
                logger().info("Mouse up invalid: has been reset to -1");
            }
        } else if (mousetype == MOUSEMOVE) {
            currentLoc = loc;
        }
    } else {
        // no hit

        if (mousetype == MOUSEDOWN) {
            logger().debug("Mouse down: no hit");
            downClicked = true;
        } else if (mousetype == MOUSEUP) {
            logger().debug("Mouse up:   no hit");
            cropActive = false;
            downClicked = false;
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
