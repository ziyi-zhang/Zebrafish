#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

class GUI : public igl::opengl::glfw::imgui::ImGuiMenu {

private:
    igl::opengl::glfw::Viewer viewer;
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> texture;
    std::vector<Eigen::MatrixXd> img;
    int stage;
    int slice;

public:
    GUI();
    void init();

protected:
    virtual void draw_menu() override;

    void draw_menu_stage0();
    void draw_menu_stage1();
    void draw_menu_stage2();
};

}  // namespace zebrafish
