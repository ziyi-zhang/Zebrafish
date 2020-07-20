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

        // upper-left corner (x0, y0)
        // lower-right corner (x1, y1)
        float x0 = baseLoc(0);
        float x1 = std::max(x0, currentLoc(0));
        float y0 = imgRows - baseLoc(1);
        float y1 = std::min(y0, imgRows - currentLoc(1));
        DrawRect(x0, y0, x1, y1, lineColor);
    }

    // showCropArea
    if (showCropArea) {
        // show the area specified by current [r0, c0] x [r1, c1]
        Eigen::MatrixXd lineColor(1, 3);
        lineColor << 0.77, 0.28, 0.24;

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
        ImGui::InputInt("Channels Per Slice", &channelPerSlice);
        ImGui::PopItemWidth();
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Depth Cropping", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::SliderInt("Slice index start", &layerBegin, 0, layerEnd);
        ImGui::SliderInt("Slice index end", &layerEnd, layerBegin, layerPerImg-1);
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Area Cropping", ImGuiTreeNodeFlags_DefaultOpen)) {

        if (ImGui::Checkbox("Enable Mouse Crop", &cropActive)) {
            if (!cropActive) 
                logger().info("Mouse crop: de-activated.");
            else
                logger().info("Mouse crop: activated.");
        }
        ImGui::Checkbox("Show cropped area", &showCropArea);
        if (ImGui::Button("Reset crop area")) {
            r0 = -1;
            c0 = -1;
            r1 = -1;
            c1 = -1;
        }

        if (ImGui::TreeNode("Detailed")) {

            ImGui::PushItemWidth(80);
            ImGui::InputInt("R0", &r0);
            ImGui::SameLine();
            ImGui::InputInt("C0", &c0);
            ImGui::InputInt("R1", &r1);
            ImGui::SameLine();
            ImGui::InputInt("C1", &c1);
            ImGui::PopItemWidth();

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

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::Button("Save")) {
        if (imagePath.empty()) {
            logger().error("Error: ImagePath is empty");
            return;  // invalid operation
                     // only entry is through "File" - "Open"
            // imagePath = FileDialog::openFileName("./.*", {"*.tif", "*.tiff"});
        }
        logger().info("Re-load image {}", imagePath);
        if (ReadTifFirstImg(imagePath, layerPerImg, channelPerSlice, img, r0, c0, r1, c1)) {
            imgRows = img[0].rows();
            imgCols = img[0].cols();

            ComputeCompressedImg();  // re-compute compressed image texture
            showCropArea = false;  // turn this into false

            logger().info("Image reloaded");
        } else {
            logger().error("Error open tiff image (reload)");
            std::cerr << "Error open tiff image (reload)" << std::endl;
        }
    }

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::Text("Stage 1: image read");
}


////////////////////////////////////////////////////////////////////////////////////////
// Crop image


void GUI::CropImage(const Eigen::Vector2f &mouse, MOUSE_TYPE mousetype) {

    static int fid;
    static Eigen::Vector3f bc, loc;
    static bool hit;

    hit = igl::unproject_onto_mesh(
        mouse, 
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
            r0 = std::round(baseLoc(1));
            r1 = std::round(currentLoc(1));
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
    viewer.data().add_edges(lineX, lineY, lineColor);
}

}  // namespace zebrafish
