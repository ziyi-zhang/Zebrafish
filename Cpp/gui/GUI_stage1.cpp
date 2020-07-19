#include <zebrafish/GUI.h>
#include <zebrafish/TiffReader.h>
#include <zebrafish/FileDialog.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/Quantile.h>

#include <igl/unproject_onto_mesh.h>
#include <algorithm>

namespace zebrafish {

namespace {

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

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 1: image read

void GUI::DrawStage1() {

    ImGui::SliderInt("Slice", &slice, 0, img.size()-1);

    viewer.data().clear();

    if (!img.empty()) {
        texture = (img[slice].array() * 255).cast<unsigned char>();
        texture.transposeInPlace();
        V << 0, 0, 0, imgCols, 0, 0, imgCols, imgRows, 0, 0, imgRows, 0;
        F << 0, 1, 2, 2, 3, 0;
        viewer.core().align_camera_center(V);
        Eigen::MatrixXd UV(4, 2);
        UV << 0, 1, 1, 1, 1, 0, 0, 0;
        viewer.data().set_mesh(V, F);
        viewer.data().set_uv(UV);
        Eigen::Vector3d ambient = Eigen::Vector3d(146./255., 172./255., 178./255.);
        Eigen::Vector3d diffuse = Eigen::Vector3d(146./255., 172./255., 178./255.);
        Eigen::Vector3d specular = Eigen::Vector3d(0., 0., 0.);
        viewer.data().uniform_colors(ambient,diffuse,specular);
        viewer.data().show_faces = true;
        viewer.data().show_lines = false;
        viewer.data().show_texture = true;
        viewer.data().set_texture(texture, texture, texture);

        Eigen::MatrixXd points(1, 3);
        points << 0, 1, 5;
        Eigen::MatrixXd pointColor(1, 3);
        pointColor << 1, 0, 0;
        viewer.data().add_points(points, pointColor);

        if (cropActive && baseLoc(0) != 0) {
            // crop activated & has been updated (not default value)
            Eigen::MatrixXd lineColor(1, 3);
            pointColor << 1, 0, 1;
            Eigen::MatrixXd lineX, lineY;
            lineX.resize(4, 3);
            lineY.resize(4, 3);
            // upper-left corner (x1, y0)
            // lower-right corner (x1, y1)
            int x0 = std::round(baseLoc(0));
            int x1 = std::max(x0, int(std::round(currentLoc(0))));
            int y0 = std::round(imgRows - baseLoc(1));
            int y1 = std::min(y0, int(std::round(imgRows - currentLoc(1))));
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
    }

    if (ImGui::Button("Print some log")) {
        logger().info("test");
    }

    ImGui::PushItemWidth(zebrafishWidth / 3.0);
    ImGui::InputInt("Layers Per Image", &layerPerImg);
    ImGui::InputInt("Channels Per Slice", &channelPerSlice);
    ImGui::PopItemWidth();

    ImGui::Separator();

    if (ImGui::Checkbox("Mouse Crop", &cropActive)) {
        if (!cropActive) 
            logger().info("Mouse crop: de-activated.");
        else
            logger().info("Mouse crop: activated.");
    }

    ImGui::PushItemWidth(80);
    ImGui::InputInt("R0", &r0);
    ImGui::SameLine();
    ImGui::InputInt("C0", &c0);

    ImGui::InputInt("R1", &r1);
    ImGui::SameLine();
    ImGui::InputInt("C1", &c1);
    ImGui::PopItemWidth();

    ImGui::Separator();

    ImGui::Text("Histogram of pixel brightness");
    ImGui::Text("Not yet implemented");
    ImGui::PushItemWidth(zebrafishWidth / 3.0);
    ImGui::InputDouble("Quantile Threshold", &normalizeQuantile);
    ImGui::PopItemWidth();

    ImGui::Separator();

    ImGui::PushItemWidth(zebrafishWidth / 3.0);
    ImGui::InputDouble("Resolution X (um)", &resolutionX);
    ImGui::InputDouble("Resolution Y (um)", &resolutionY);
    ImGui::InputDouble("Resolution Z (um)", &resolutionZ);
    ImGui::PopItemWidth();

    ImGui::Separator();

    ImGui::Text("Bspline config ... ...");
    ImGui::Text("Not yet implemented");

    ImGui::Separator();

    //////////////////////////////////////////////////////////////////////////

    if (ImGui::Button("(Re)load image")) {
        if (imagePath.empty()) {
            imagePath = FileDialog::openFileName("./.*", {"*.tif", "*.tiff"});
        }
        logger().info("(Re)load image {}", imagePath);
        if (ReadTifFirstImg(imagePath, layerPerImg, channelPerSlice, img, r0, c0, r1, c1)) {
            imgRows = img[0].rows();
            imgCols = img[0].cols();
            // In case the tiff image is very small
            layerPerImg = img.size();
            logger().info("Image reloaded");
            NormalizeImage(img);
        } else {
            logger().error("Error open tiff image (reload)");
            std::cerr << "Error open tiff image (reload)" << std::endl;
        }
    }

    ImGui::Text("Stage 1");
    // viewer.core().set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_TRACKBALL);
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
            baseLoc = loc;
            currentLoc = loc;
        } else if (mousetype == MOUSEUP) {
            logger().debug("Mouse up:   fid = {} | bc = {:.2f} x {:.2f} x {:.2f} | loc = {:.2f} x {:.2f} x {:.2f}", fid, bc(0), bc(1), bc(2), loc(0), loc(1), loc(2));
            cropActive = false;

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
        } else if (mousetype == MOUSEUP) {
            logger().debug("Mouse up:   no hit");
            cropActive = false;
        }
    }
}

}  // namespace zebrafish
