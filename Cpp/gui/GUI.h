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
///       |   Grid Search  |   Optimization
/// alive | x y z r energy | x y z r energy iter

    int num;
    Eigen::Matrix<bool,   Eigen::Dynamic, 1> alive;
    Eigen::Matrix<double, Eigen::Dynamic, 5> grid_search;
    Eigen::Matrix<double, Eigen::Dynamic, 6> optimization;

    pointRecord_t() : num(0) {}
} pointRecord_t;


typedef struct hist_t {

    Eigen::MatrixXf hist;
    double minValue, maxValue;
    // Note: not the min/max of "hist"
    //       this is the min/max of the data used to calculate the hist
    //       In other words, the x-axis min and max
} hist_t;

////////////////////////////////////////////////////////
// GUI

class GUI : public igl::opengl::glfw::imgui::ImGuiMenu {

private:
    igl::opengl::glfw::Viewer viewer;
    int stage;  // stage in zebrafish_panel
    int histBars;  // number of bars in histogram

    //////////////////////////////////////////////////
    // core
    bspline bsplineSolver;

    pointRecord_t pointRecord;
    ///       |   Grid Search  |   Optimization
    /// alive | x y z r energy | x y z r energy iter

    //////////////////////////////////////////////////
    // image (imageData)
    std::string imagePath;
    imageData_t imgData;
    image_t img;
    int imgRows, imgCols;
    int layerPerImg, channelPerSlice;
    int layerBegin, layerEnd;  // only slices in this interval in each 3D image will be computed
    double resolutionX, resolutionY, resolutionZ;
    int slice;  // which slice in the 3D image to show
    float normalizeQuantile, normalizeQuantileRes;
    // Hist
    hist_t imgHist;

    //////////////////////////////////////////////////
    // Grid Search
    Eigen::MatrixXd gridSampleInput, gridSampleOutput;
    double gapX_grid, gapY_grid, gapZ_grid;
    double rArrayMin_grid, rArrayMax_grid, rArrayGap_grid;
    float gridEnergyThres;  // energy threshold to decide whether a starting point is worth optimizing
    bool showPromisingPoints;
    Eigen::MatrixXd promisingPointLoc;  // visualization purpose
    // Hist
    hist_t gridEnergyHist;

    //////////////////////////////////////////////////
    // Optimization
    float optimEnergyThres;  // only cylinders with energy smaller than this value would have the chance to go to "filter" stage
    bool showOptimizedPoints;
    Eigen::MatrixXd optimPointLoc;  // visualization purpose
    double optimEpsilon, optimMaxIt;
    // Hist
    hist_t optimEnergyHist;

    //////////////////////////////////////////////////
    // Cylinder Filter
    float cylinderEnergyThres, cylinderRadiusThres;
    int cylinderIterThres;
    bool showCylFilterPoints;
    Eigen::MatrixXd cylPointLoc;  // visualization purpose
    // Hist
    hist_t cylEnergyHist, cylRadiusHist, cylIterHist;

    //////////////////////////////////////////////////
    // 3D image viewer
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    int imageViewerType;

    //////////////////////////////////////////////////
    // crop image
    bool cropActive;  // user use mouse to crop a rectangular area
    bool downClicked;  // helper variable used by "cropActive"
    bool showCropArea;  // visualize the current [r0, c0] x [r1, c1] area
    Eigen::Vector3f baseLoc, currentLoc;
    int r0, c0, r1, c1;  // upper-left [r0, c0]

    //////////////////////////////////////////////////
    // property editor
    int propertyListType;  // 0 for compressed, 1 for per slice

    //////////////////////////////////////////////////
    // visualization
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> texture, compressedImgTexture;
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

    //////////////////////////////////////////////////
    // bool flag indicating moving from a stage to another
    bool stage1to2Flag;
    bool stage4to5Flag;

public:
    std::ostringstream oss;

    GUI();
    void init(std::string imagePath);

protected:
    void draw_menu() override;
    void post_resize(int w, int h) override;
    bool MouseDownCallback(igl::opengl::glfw::Viewer &viewer, int button, int modifier);
    bool MouseUpCallback(igl::opengl::glfw::Viewer &viewer, int button, int modifier);
    bool MouseMoveCallback(igl::opengl::glfw::Viewer &viewer, int mouse_x, int mouse_y);

    void DrawStage1();
    void DrawStage2();
    void DrawStage3();
    void DrawStage4();
    void DrawStage5();
    void DrawStage6();

private:
    const int stageMax = 6;  // 6 stages

    void Draw3DImage();
    void ComputeCompressedImg(const image_t &img_);
    void DrawMainMenuBar();
    void DrawZebrafishPanel();
    void DrawMenuFile();
    void DrawMenuWindow();
    
    void DrawWindowLog();
    void DrawWindow3DImageViewer();
    void DrawWindowPropertyEditor();
    void DrawWindowGraphics();

    //////////////////////////////////////////////////
    // Stage 1
    void CropImage(const Eigen::Vector2f &mouse, MOUSE_TYPE mousetype);
    void DrawRect(double x0, double y0, double x1, double y1, const Eigen::MatrixXd &lineColor);

    //////////////////////////////////////////////////
    // Stage 2
    void ComputeImgHist(const image_t &img_);

    //////////////////////////////////////////////////
    // Stage 3
    void GridSearch();
    void UpdateSampleNewton(const Eigen::MatrixXd &gridSampleInput, const Eigen::MatrixXd &gridSampleOutput);
        /// This function will clear "pointRecord" and intialize it 
        /// with grid search results
    void UpdatePromisingPointLoc();
    void UpdateGridEnergyHist();

    //////////////////////////////////////////////////
    // Stage 4
    void Optimization();
    void UpdateOptimPointLoc();
    void UpdateOptimEnergyHist();

    //////////////////////////////////////////////////
    // Stage 5.1 Cylinder Filter
    void CylinderFilter();
    void UpdateCylPointLoc();
    void UpdateCylEnergyHist();
    void UpdateCylRadiusHist();
    void UpdateCylIterHist();

    //////////////////////////////////////////////////
    // Stage 5.2 Cluster Filter
    void Cluster();

    //////////////////////////////////////////////////
    // Stage 6

};

}  // namespace zebrafish
