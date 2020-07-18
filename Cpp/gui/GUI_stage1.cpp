#include <zebrafish/GUI.h>
#include <zebrafish/TiffReader.h>
#include <zebrafish/FileDialog.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/Quantile.h>

#include <igl/unproject_on_plane.h>


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
        int xMax = img[slice].cols();
        int yMax = img[slice].rows();
        Eigen::MatrixXd V(4, 3);
        V << 0, 0, 0, xMax, 0, 0, xMax, yMax, 0, 0, yMax, 0;
        viewer.core().align_camera_center(V);

        Eigen::MatrixXi F(2, 3);
        F << 0, 1, 2, 2, 3, 0;

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
        points << 0, 1, 0;
        Eigen::MatrixXd pointColor(1, 3);
        pointColor << 1, 0, 0;
        viewer.data().add_points(points, pointColor);
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

}  // namespace zebrafish
