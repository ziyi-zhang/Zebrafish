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
    int stage;  // stage in zebrafish_panel
    int slice;  // which slice in the 3D image to show

    //////////////////////////////////////////////////
    int windowWidth = 1600;
    int windowHeight = 900;
    int zebrafishWidth = 300;

    double mainMenuHeight = 0;

    // bool flag indicating whether the panel is being rendered
    bool show_log = false;
    bool show_graphics = false;

public:
    GUI();
    void init();

protected:
    void draw_menu() override;
    void post_resize(int w, int h) override;

    void DrawStage0();
    void DrawStage1();
    void DrawStage2();

private:
    const int stageMax = 3;

    void DrawMainMenuBar();
    void DrawZebrafishPanel();
    void DrawMenuFile();
    void DrawMenuWindow();
    
    void DrawWindowGraphics();
    void DrawWindowLog();
};

}  // namespace zebrafish
