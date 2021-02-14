#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/Cylinder.h>
#include <zebrafish/autodiff.h>
#include <zebrafish/ICP.h>
#include <zebrafish/Padding.h>

#include <string>
#include <vector>
#include <sstream>
#include <map>
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

////////////////////////////////////////////////////////////////////////

typedef struct hist_t {

    Eigen::MatrixXf hist;
    double minValue, maxValue;
    // Note: not the min/max of "hist"
    //       this is the min/max of the data used to calculate the hist
    //       In other words, the x-axis min and max
} hist_t;


typedef struct crop_t {

    bool cropActive;  // user use mouse to crop a rectangular area
    bool downClicked;  // helper variable used by "cropActive"
    bool showCropArea;  // visualize the current [r0, c0] x [r1, c1] area
    Eigen::Vector3f baseLoc, currentLoc;
    int r0, c0, r1, c1;  // upper-left [r0, c0]

    crop_t() : cropActive(false), downClicked(false), showCropArea(true), r0(-1), c0(-1), r1(-1), c1(-1) {
        baseLoc << 0.0f, 0.0f, 0.0f;
        currentLoc << 0.0f, 0.0f, 0.0f;
    }
} crop_t;


typedef struct grid_t {

    double gapX, gapY, gapZ;
    double rArrayMin, rArrayMax, rArrayGap;

    bool skipMembrane;
    double membraneThres;  // brightness threshold to decide whether a point is in the membrane area

    float energyThres;  // energy threshold to decide whether a starting point is worth optimizing

    grid_t() : gapX(1.0), gapY(1.0), gapZ(1.0), rArrayMin(3.0), rArrayMax(6.0), rArrayGap(1.0), skipMembrane(false), membraneThres(), energyThres(-0.1) {}
} grid_t;


typedef struct ICP_t {

    int patternRows, patternCols, patternRef;
    double patternSpacing;
    float xDisp, yDisp, angleRot, scale, xStretch, yStretch, xdistort, ydistort;
    RMat_t Rmat;  // rotation matrix
    TMat_t Tmat;  // translation matrix
    Eigen::VectorXi matchIdx;  // p[i] corresponds to q[ matchIdx[i] ]

    // reference pattern
    Eigen::MatrixXd refV, refV_aligned;  // #ICP.refV * 3 reference point locations
    Eigen::MatrixXi refV_RC;  // row & col index of corresponding vertex

    ICP_t() : patternRows(0), patternCols(0), patternSpacing(0.0), 
              xDisp(0.0), yDisp(0.0), angleRot(0.0), scale(1.0), 
              xStretch(0.0), yStretch(0.0), xdistort(0.0), ydistort(0.0) {
        Rmat = Eigen::MatrixXd::Identity(3, 3);
        Tmat = Eigen::MatrixXd::Zero(3, 1);
    }
} ICP_t;


/*
typedef struct padding_t {

    bool enable, success;
    
    // Eigen::MatrixXd padV
} padding_t;
*/


typedef struct analysis_t {

    double offset;  // Diagonal multiplier for box mesh
    double radius_edge_ratio;  // Radius edge ratio used by tetgen
    double max_tet_vol;  // Maximum tet vol used by tetgen
    double E;  // Young's modulus (default 566.7Pa)
    double nu;  // Poisson's ratio
    bool is_linear;  // Use non-linear material
    int discr_order;  // Analysis discretization order
    int n_refs;  // Number of mesh uniform refinements
    double vismesh_rel_area;  // Desnsity of the output visualization
    int upsample;  // upsample for a denser mesh

    std::vector<Eigen::MatrixXd> V;  // only used by analysis input file
    Eigen::MatrixXi F;  // only used by analysis input file
} analysis_t;


typedef struct UIsize_t {
    int windowWidth;
    int windowHeight;
    int zebrafishWidth;
    int logHeight;
    int Image3DViewerHeight;
    int RHSPanelWidth;
    double mainMenuHeight;
} UTsize_t;


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
    bool reverseColor;

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
    double previewQuantileBrightness;  // used in stage-1 preview
    std::string maskPath;
    std::string analysisInputPath;
    // Hist
    hist_t imgHist;

    //////////////////////////////////////////////////
    // membrane mask
    image_t membraneMask;
    bool membraneMaskLoad;
    bool membraneMaskCylApply;
    bool membraneMaskClusterApply;
    float maskMax, maskThres;

    //////////////////////////////////////////////////
    // B-spline
    int bsplineDegree;
    double bsplineSolverTol;

    //////////////////////////////////////////////////
    // Grid Search
    grid_t grid;
    Eigen::MatrixXd gridSampleInput, gridSampleOutput;
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
    float clusterSizeThres;
    bool showClusterFilterPoints;
    Eigen::MatrixXd clusterPointLoc;  // visualization purpose
    // Hist
    hist_t clusterSizeHist;

    //////////////////////////////////////////////////
    // ICP
    bool showMarkerPoints, showReferencePoints, showICPLines, showMarkerMesh;
    std::string patternFilename;
    ICP_t ICP;
    Eigen::MatrixXd refPointLoc;  // visualization purpose
    Eigen::MatrixXi markerMeshArray;  // F matrix for mesh
    int meshID;

    //////////////////////////////////////////////////
    // Padding
    RCMap_t markerRCMap;  // map marker index in "markerMeshArray" to (row, col) in regular triangular mesh

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
    bool secondRoundDepthCorrection;
    int depthCorrectionNum;
    float depthCorrectionGap;
    float optimMaxXYDisp;  // max XY displacement during optimization

    //////////////////////////////////////////////////
    // Analysis
    analysis_t analysisPara;

    //////////////////////////////////////////////////
    // 3D image viewer
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    int imageViewerType;
    int imageViewerCompressType;
    float imageViewerDarkenFactor_max, imageViewerDarkenFactor_avg;

    //////////////////////////////////////////////////
    // [mouse pick] crop image
    crop_t imageCrop;  // stage-1 crop image in XY
    crop_t meanCrop;  // stage-8 crop image for average displacement area

    //////////////////////////////////////////////////
    // [mouse pick] manually reject clusters
    bool rejectActive;
    bool rejectHit;  // checked when "MOUSEMOVE", return true if there is a hit with clusters
    int rejectMode;
    double mousePickDistSquareThres;  // in pixels
    Eigen::VectorXi rejectHitIndex;

    //////////////////////////////////////////////////
    // [mouse pick] manually drag markers
    bool markerDragActive;
    bool markerDragHit;  // whether there is a hit with a marker
    bool markerDragFocused;  // whether a marker has been chosen
    int markerDragHitIndex;
    Eigen::Vector2f markerDragLoc;

    //////////////////////////////////////////////////
    // property editor
    int propertyListType;

    //////////////////////////////////////////////////
    // visualization
    texture_t texture;
    textureArray_t compressedImgTextureArray;
    std::vector<std::vector<int> > markerDepthCorrectionSuccess;
    std::vector<Eigen::MatrixXd> markerPointLocArray;  // visualization purpose
    Eigen::Matrix<bool, Eigen::Dynamic, 1> markerPointStatusArray;  // visualization purpose
    bool manualOverrideMarkerVis, showAllMarkers;
    int sliceToShow;  // which slice in the 3D image to show
    int frameToShow;  // which frame to show
    int currentLoadedFrames;
    UIsize_t UIsize;
    bool UIsize_redraw;  // redraw after window resize
    bool show_refPoints;
    bool show_allMarkerIndex;
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
    bool stage2to3Flag;
    bool stage4to5Flag;
    bool stage5to6Flag;

    //////////////////////////////////////////////////
    // stage lock
    bool stage1Lock;
    bool stage2Lock;
    bool stage3Lock;
    bool stage4Lock;
public:
    std::ostringstream oss;

    GUI();
    void init(std::string imagePath, std::string maskPath, std::string analysisInputPath, int debugMode, bool NoGUI);

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
    bool MarkerRecursiveDepthCorrection(int frameIdx, int depthNum, double depthGap, bool logEnergy = false, bool forceSecondRound = false);
    bool MarkerDepthCorrection(int frameIdx, int num, double gap, bool logEnergy = false, bool logSuccess = false);
    void NormalizeImage(image_t &image, double thres);
    void LoadPreviewImage(std::string path);
    void StateChangeReset();

    // visualization
    void DrawText(Eigen::Vector3d pos, const std::string &text, const Eigen::Vector4f color);  // DO NOT USE THIS
    void DrawReferenceDots();
    void ShowAllMarkerIndex();

    // marker drag
    void RenderMarkerDragGUI();
    void MouseSelectMarker(const Eigen::Vector2f &mouse);
    void MarkerDragSelect();
    void MarkerDragSetNewLoc();
    void MarkerDragReset();
    void MarkerDragVisualization();

    // mask
    bool PointInMaskArea(double x, double y, double z);

    //////////////////////////////////////////////////
    // Stage 1 Image Read
    void ImageReadReset();
    void CropImage(const Eigen::Vector2f &mouse, MOUSE_TYPE mousetype, crop_t &cropStruct);
    void DrawRect(double x0, double y0, double x1, double y1, const Eigen::MatrixXd &lineColor);

    //////////////////////////////////////////////////
    // Stage 2 Pre-process & B-spline
    void ComputeImgHist(const image_t &img_);

    //////////////////////////////////////////////////
    // Stage 3 Grid Search
    void GridSearch();
    bool InMembraneArea(const image_t &image, const double thres, double x, double y, double z, double r);
    bool ValidGridSearchPoint(const image_t &image, const bspline &bsp, bool skipMembrane, double membraneThres, double x, double y, double z, double r);
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
    bool ClusterNearBorderWarn();

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
    void GetMarkersInAvgDispArea(std::vector<bool> &markerInAvgDispArea);
    bool OptimizeAllFrames(bool logEnergy);
    void ApplyOpticalFlow(int prevFrameIdx);
    void OptimizeOneFrame(int prevFrameIdx);
    bool SaveMeshToVTU_point(bool onlySaveFirstFrameMesh, bool saveAccumulativeDisplacement, bool saveAccumulativeDisplacement_relative, bool saveIncrementalDisplacement, bool saveIncrementalDisplacement_relative);
    bool SaveMeshToVTU_cell(bool onlySaveFirstFrameMesh);
    bool SaveImageToTIFF(bool saveMarkerImage, bool saveCellImage, int cellChannel);
    void SaveMeshToOBJ();
    void SaveDisplacementToTXT();
};

}  // namespace zebrafish
