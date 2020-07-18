#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/Cylinder.h>
#include <zebrafish/autodiff.h>

#include <string>
#include <sstream>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

typedef struct pointRecord_t {

    int num;
    Eigen::Matrix<bool,   Eigen::Dynamic, 1> alive;
    Eigen::Matrix<double, Eigen::Dynamic, 5> grid_search;
    Eigen::Matrix<double, Eigen::Dynamic, 6> optimization;

    pointRecord_t() : num(0) {}
} pointRecord_t;

////////////////////////////////////////////////////////
// GUI

class GUI : public igl::opengl::glfw::imgui::ImGuiMenu {

private:
    igl::opengl::glfw::Viewer viewer;
    int stage;  // stage in zebrafish_panel

    //////////////////////////////////////////////////
    // core
    bspline bsplineSolver;

    double gridEnergyThres;
    pointRecord_t pointRecord;
    ///       |   Grid Search  |   Optimization
    /// alive | x y z r energy | x y z r energy iter

    //////////////////////////////////////////////////
    // image (imageData)
    std::string imagePath;
    imageData_t imgData;
    image_t img;
    int layerPerImg, channelPerSlice;
    double resolutionX, resolutionY, resolutionZ;
    int slice;  // which slice in the 3D image to show
    double normalizeQuantile;

    //////////////////////////////////////////////////
    // crop image
    bool cropActive;
    int clickCount;
    int r0, c0, r1, c1;  // upper-left [r0, c0]

    //////////////////////////////////////////////////
    // property editor
    int propertyListType;

    //////////////////////////////////////////////////
    // visualization
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> texture;
    int windowWidth;
    int windowHeight;
    int zebrafishWidth;
    int logHeight;
    int Image3DViewerHeight;
    int RHSPanelWidth;
    double mainMenuHeight;

    //////////////////////////////////////////////////
    // bool flag indicating whether the panel is being rendered
    bool show_log;
    bool show_3DImage_viewer;
    bool show_property_editor;
    bool show_graphics;

public:
    std::ostringstream oss;

    GUI();
    void init(std::string imagePath);

protected:
    void draw_menu() override;
    void post_resize(int w, int h) override;
    bool MouseDownCallback(igl::opengl::glfw::Viewer &viewer, int button, int modifier);

    void DrawStage1();
    void DrawStage2();
    void DrawStage3();

private:
    const int stageMax = 3;  // 3 stages

    void DrawMainMenuBar();
    void DrawZebrafishPanel();
    void DrawMenuFile();
    void DrawMenuWindow();
    
    void DrawWindowLog();
    void DrawWindow3DImageViewer();
    void DrawWindowPropertyEditor();
    void DrawWindowGraphics();

    //////////////////////////////////////////////////
    // Stage 2
    void GridSearch();
    void Optimization();
    void UpdateSampleNewton(const Eigen::MatrixXd &gridSampleInput, const Eigen::MatrixXd &gridSampleOutput);
    /// This function will clear "pointRecord" and intialize it 
    /// with grid search results
};

}  // namespace zebrafish
