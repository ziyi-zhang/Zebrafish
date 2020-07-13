#include <zebrafish/GUI.h>
#include <zebrafish/FileDialog.h>


namespace zebrafish {

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 1: 

void GUI::DrawStage1() {

    // if(ImGui::SliderInt2("Slide"))
    viewer.data().clear();
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_NO_ROTATION);

    texture = (img[slice].array()).cast<unsigned char>();

    int xMax = img[slice].cols();
    int yMax = img[slice].rows();
    Eigen::MatrixXd V(4, 3);
    V << 0, 0, 0, yMax, 0, 0, yMax, xMax, 0, 0, xMax, 0;
    viewer.core().align_camera_center(V);

    Eigen::MatrixXi F(2, 3);
    F << 0, 1, 2, 2, 3, 0;

    Eigen::MatrixXd UV(4, 2);
    UV << 0, 1, 1, 1, 1, 0, 0, 0;
    viewer.data().set_mesh(V, F);
    viewer.data().set_uv(UV);
    viewer.data().show_faces = true;
    viewer.data().show_lines = false;
    viewer.data().show_texture = true;
    viewer.data().set_texture(texture, texture, texture);


    ImGui::Text("Stage 1");
}

}  // namespace zebrafish
