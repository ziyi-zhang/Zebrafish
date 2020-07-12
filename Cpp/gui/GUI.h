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

    void DrawMainMenuBar();
    void DrawMenuFile();
    void DrawMenuWindow();
    void DrawWindowGraphics();

    bool show_graphics = false;

public:
    GUI();
    void init();

protected:
    virtual void draw_menu() override;

    void DrawStage0();
    void DrawStage1();
    void DrawStage2();
};

}  // namespace zebrafish
