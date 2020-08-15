#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/Cylinder.h>
#include <zebrafish/autodiff.h>
#include <zebrafish/ICP.h>

#include <string>
#include <vector>
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


typedef struct clusterRecord_t {
///       |         Loc             |   energy   |  size
/// alive | meanX meanY meanZ meanR | meanEnergy | size(count)

    int num;
    Eigen::Matrix<bool,   Eigen::Dynamic, 1> alive;
    Eigen::Matrix<double, Eigen::Dynamic, 4> loc;
    Eigen::Matrix<double, Eigen::Dynamic, 1> energy;
    Eigen::Matrix<int, Eigen::Dynamic, 1>    size;

    clusterRecord_t() : num(0) {}
} clusterRecord_t;


typedef struct markerRecord_t {
///   Loc   | energy |   size
/// X Y Z R | energy | size(count)

    int num;
    Eigen::Matrix<double, Eigen::Dynamic, 4> loc;
    Eigen::Matrix<double, Eigen::Dynamic, 1> energy;
    Eigen::Matrix<int, Eigen::Dynamic, 1>    size;

    markerRecord_t() : num(0) {}
} markerRecord_t;


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
    // shared
    igl::opengl::glfw::Viewer viewer;
    imageData_t imgData;
    int stage;  // stage in zebrafish_panel
    int histBars;  // number of bars in histogram
    bool showBackgroundImage;
    bool showTooltip;
    float lineWidth;

    //////////////////////////////////////////////////
    // core
    std::vector<bspline> bsplineArray;

    pointRecord_t pointRecord;
    ///       |   Grid Search  |   Optimization
    /// alive | x y z r energy | x y z r energy iter

    clusterRecord_t clusterRecord;
    ///       |         Loc             |   size
    /// alive | meanX meanY meanZ meanR | size(count)

    std::vector<markerRecord_t> markerArray;
    ///   Loc   | energy |   size
    /// X Y Z R | energy | size(count)

    //////////////////////////////////////////////////
    // image (image metadata)
    std::string imagePath;
    int ttlFrames;
    int imgRows, imgCols;
    int layerPerImg, channelPerSlice;
    int channelToLoad;
    int layerBegin, layerEnd;  // only slices in this interval in each 3D image will be computed and visualized
    double resolutionX, resolutionY, resolutionZ;
    float normalizeQuantile, normalizeQuantileRes;
    // Hist
    hist_t imgHist;

    //////////////////////////////////////////////////
    // B-spline
    int bsplineDegree;
    double bsplineSolverTol;

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
    // Cluster Filter
    float clusterDistThres, finalizeClusterDistThres;
    int clusterSizeThres;
    bool showClusterFilterPoints;
    Eigen::MatrixXd clusterPointLoc;  // visualization purpose
    // Hist
    hist_t clusterSizeHist;

    //////////////////////////////////////////////////
    // ICP
    bool showMarkerPoints, showReferencePoints, showICPLines, showMarkerMesh;
    std::string patternFilename;
    int ICP_patternRef, ICP_patternRows, ICP_patternCols;
    double ICP_patternSpacing;
    float ICP_xDisp, ICP_yDisp, ICP_angleRot, ICP_scale;
    RMat_t ICP_Rmat;  // rotation matrix
    TMat_t ICP_Tmat;  // translation matrix
    Eigen::MatrixXd refPointLoc;  // visualization purpose
    Eigen::MatrixXd refV, refV_aligned;  // #refV * 3 reference point locations
    Eigen::MatrixXi markerMeshArray;  // F matrix for mesh
    int meshID;
    Eigen::VectorXi ICP_matchIdx;  // p[i] corresponds to q[ ICP_matchIdx[i] ]

    //////////////////////////////////////////////////
    // Optical Flow
    int desiredFrames;  // only load the first "desiredFrames" frames
    double opticalFlowAlpha;  // weight factor
    int opticalFlowIter;
    bool showOpticalFlow;
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, 3> > opticalFlowCorrection;
        /// opticalFlowCorrection[prevFrameIdx] is an #N x 3 matrix
        /// i-th matrix correspond to the correction from i-th frame to (i+1)-th frame

    //////////////////////////////////////////////////
    // Displacement
    int depthCorrectionNum;
    float depthCorrectionGap;

    //////////////////////////////////////////////////
    // 3D image viewer
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    int imageViewerType;
    int imageViewerCompressType;
    float imageViewerDarkenFactor_max, imageViewerDarkenFactor_avg;

    //////////////////////////////////////////////////
    // [mouse pick] crop image
    bool cropActive;  // user use mouse to crop a rectangular area
    bool downClicked;  // helper variable used by "cropActive"
    bool showCropArea;  // visualize the current [r0, c0] x [r1, c1] area
    Eigen::Vector3f baseLoc, currentLoc;
    int r0, c0, r1, c1;  // upper-left [r0, c0]

    //////////////////////////////////////////////////
    // [mouse pick] manually reject clusters
    bool rejectActive;
    bool rejectHit;  // checked when "MOUSEMOVE", return true if there is a hit with clusters
    int rejectMode;
    double mousePickDistSquareThres;  // in pixels
    Eigen::VectorXi rejectHitIndex;

    //////////////////////////////////////////////////
    // property editor
    int propertyListType;

    //////////////////////////////////////////////////
    // visualization
    texture_t texture;
    textureArray_t compressedImgTextureArray;
    std::vector<Eigen::MatrixXd> markerPointLocArray;  // visualization purpose
    Eigen::Matrix<bool, Eigen::Dynamic, 1> markerPointStatusArray;  // visualization purpose
    bool manualOverrideMarkerVis, showAllMarkers;
    int sliceToShow;  // which slice in the 3D image to show
    int frameToShow;  // which frame to show
    int currentLoadedFrames;
    int windowWidth;
    int windowHeight;
    int zebrafishWidth;
    int logHeight;
    int Image3DViewerHeight;
    int RHSPanelWidth;
    double mainMenuHeight;
    // color
    Eigen::MatrixXd markerPointColor;

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
    bool stage5to6Flag;

public:
    std::ostringstream oss;

    GUI();
    void init(std::string imagePath, int debugMode);

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
    void DrawStage7();
    void DrawStage8();

private:
    const int stageMax = 8;  // 8 stages (1..8)

    void Draw3DImage();
    void DrawMarkerMesh();
    void DrawMainMenuBar();
    void DrawZebrafishPanel();
    void DrawMenuFile();
    void DrawMenuWindow();
    void DrawMenuHelp();
    
    void DrawWindowLog();
    void DrawWindow3DImageViewer();
    void DrawWindowPropertyEditor();
    void DrawWindowGraphics();

    // shared
    void ComputeCompressedTextureAvg(const image_t &img_, int index);
    void ComputeCompressedTextureMax(const image_t &img_, int index);
    void ComputeCompressedTextureForAllLoadedFrames();
    void UpdateMarkerPointLocArray();
    void MarkerDepthCorrection(int frameIdx, int num = 2, double gap = 0.5);
    void NormalizeImage(image_t &image, double thres);

    //////////////////////////////////////////////////
    // Stage 1 Image Read
    void ImageReadReset();
    void CropImage(const Eigen::Vector2f &mouse, MOUSE_TYPE mousetype);
    void DrawRect(double x0, double y0, double x1, double y1, const Eigen::MatrixXd &lineColor);

    //////////////////////////////////////////////////
    // Stage 2 Pre-process & B-spline
    void ComputeImgHist(const image_t &img_);

    //////////////////////////////////////////////////
    // Stage 3 Grid Search
    void GridSearch();
    void UpdateSampleNewton(const Eigen::MatrixXd &gridSampleInput, const Eigen::MatrixXd &gridSampleOutput);
        /// This function will clear "pointRecord" and intialize it 
        /// with grid search results
    void UpdatePromisingPointLoc();
    void UpdateGridEnergyHist();

    //////////////////////////////////////////////////
    // Stage 4 Optimization
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
    void ClusterFilter();
    void UpdateClusterPointLoc();
    void UpdateClusterSizeHist();
    void FinalizeClusterLoc();
    void MouseSelectCluster(const Eigen::Vector2f &mouse);
    void MouseRejectCluster();

    //////////////////////////////////////////////////
    // Stage 6 Iterative Closest Point
    void InitializeICPPattern();
    void GenerateICPPattern();
    void ResetICPTransformation();
    void SearchICP();
    void UpdateMarkerMesh();
    void PreprocessPatternLoc();
    void UpdateRefPointManualAlignment();
    void UpdateRefPointLoc();
    void DrawICPLines();

    //////////////////////////////////////////////////
    // Stage 7 Optical Flow
    void LoadSubsequentFrames();
    void ComputeBsplineForAllFrames();
    void RunOpticalFlow();

    //////////////////////////////////////////////////
    // Stage 8 Displacement
    void OptimizeAllFrames();
    void OptimizeOneFrame(int prevFrameIdx);
    bool SaveMeshToVTU();
    void SaveImageToTIFF();
};

}  // namespace zebrafish
