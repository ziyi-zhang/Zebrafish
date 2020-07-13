#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>

#include <string>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

class GUI : public igl::opengl::glfw::imgui::ImGuiMenu {

private:
    // Core private variables
    igl::opengl::glfw::Viewer viewer;

    int stage;  // stage in zebrafish_panel

    // image (imageData)
    std::string imagePath;
    imageData_t imgData;
    image_t img;
    int layerPerImg, channelPerSlice;
    int slice;  // which slice in the 3D image to show

    //////////////////////////////////////////////////
    // visualization
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> texture;
    int windowWidth = 1600;
    int windowHeight = 900;
    int zebrafishWidth = 300;
    int logHeight = 150;
    int RHSPanelWidth = 300;
    double mainMenuHeight;

    // bool flag indicating whether the panel is being rendered
    bool show_log = false;
    bool show_graphics = false;

public:
    GUI();
    void init(std::string imagePath);

protected:
    void draw_menu() override;
    void post_resize(int w, int h) override;

    void DrawStage0();
    void DrawStage1();
    void DrawStage2();

private:
    const int stageMax = 2;  // 3 stages

    void DrawMainMenuBar();
    void DrawZebrafishPanel();
    void DrawMenuFile();
    void DrawMenuWindow();
    
    void DrawWindowGraphics();
    void DrawWindowLog();
};

}  // namespace zebrafish
