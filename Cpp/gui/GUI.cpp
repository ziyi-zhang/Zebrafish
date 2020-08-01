#include <zebrafish/GUI.h>
#include <zebrafish/TiffReader.h>
#include <zebrafish/FileDialog.h>
#include <zebrafish/Logger.hpp>

#include <string>
#include <sstream>
#include <algorithm>


namespace zebrafish {

namespace {

struct PropertyEditorItem {
/// Used by "Property Editor"

    static void AppendPointRecordItem(const char* prefix, int uid, const pointRecord_t &pointRecord) {

        ImGui::PushID(uid);
        ImGui::AlignTextToFramePadding();
        bool nodeOpen = ImGui::TreeNode("Object", "%s %u", prefix, uid);
        ImGui::NextColumn();
        ImGui::AlignTextToFramePadding();

        if (pointRecord.optimization(0, 0) == 0)
            ;
        else if (pointRecord.alive(uid))
            ImGui::Text("valid");
        else 
            ImGui::Text("invalid");

        ImGui::NextColumn();

        if (nodeOpen) {
            static const std::vector<std::string> itemName{"Energy", "x", "y", "z", "r", "Energy", "x", "y", "z", "r", "Iter"};
            for (int i=0; i<5; i++) {
                ImGui::PushID(i);
                ImGui::AlignTextToFramePadding();
                ImGui::TreeNodeEx(itemName[i].c_str(), ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_Bullet);
                ImGui::NextColumn();
                if (i == 0) {
                    ImGui::Text("%.4f", pointRecord.grid_search(uid, 4));
                } else {
                    ImGui::Text("%.3f", pointRecord.grid_search(uid, i-1));
                }
                ImGui::NextColumn();
                ImGui::PopID();
            }

            ImGui::Separator(); ///////////////////////

            for (int i=0; i<6; i++) {
                ImGui::PushID(i);
                ImGui::AlignTextToFramePadding();
                ImGui::TreeNodeEx(itemName[i+5].c_str(), ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_Bullet);
                ImGui::NextColumn();
                if (i == 0) {
                    ImGui::Text("%.4f", pointRecord.optimization(uid, 4));
                } else if (i == 5) {
                    ImGui::Text("%.0f", pointRecord.optimization(uid, 5));
                } else {
                    ImGui::Text("%.3f", pointRecord.optimization(uid, i-1));
                }
                ImGui::NextColumn();
                ImGui::PopID();
            }

            ImGui::Separator(); ///////////////////////
            ImGui::TreePop();
        }

        ImGui::PopID();
    }

    // ----------------------------------------------------------------------------------------------

    static bool AppendClusterRecordItem(const char* prefix, int uid, clusterRecord_t &clusterRecord) {

        bool res = false;

        ImGui::PushID(uid);
        ImGui::AlignTextToFramePadding();
        bool nodeOpen = ImGui::TreeNode("Object", "%s %u", prefix, uid);
        ImGui::NextColumn();
        ImGui::AlignTextToFramePadding();
        if (clusterRecord.alive(uid)) {
            if (ImGui::Checkbox(" valid", &clusterRecord.alive(uid))) {
                res = true;
            }
        } else {
            if (ImGui::Checkbox(" inalid", &clusterRecord.alive(uid))) {
                res = true;
            }
        }

        ImGui::NextColumn();

        if (nodeOpen) {
            static const std::vector<std::string> itemName{"meanX", "meanY", "meanZ", "meanR", "energy", "size"};
            for (int i=0; i<4; i++) {
                ImGui::PushID(i);
                ImGui::AlignTextToFramePadding();
                ImGui::TreeNodeEx(itemName[i].c_str(), ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_Bullet);
                ImGui::NextColumn();
                    
                ImGui::Text("%.3f", clusterRecord.loc(uid, i));

                ImGui::NextColumn();
                ImGui::PopID();
            }

            ImGui::Separator(); ///////////////////////

            for (int i=4; i<=5; i++) {
                ImGui::PushID(i);
                ImGui::AlignTextToFramePadding();
                ImGui::TreeNodeEx(itemName[i].c_str(), ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_Bullet);
                ImGui::NextColumn();
                if (i == 4) {
                    ImGui::Text("%.4f", clusterRecord.energy(uid));
                } else {
                    ImGui::Text("%d", clusterRecord.size(uid));
                }
                ImGui::NextColumn();
                ImGui::PopID();
            }

            ImGui::Separator(); ///////////////////////
            ImGui::TreePop();
        }

        ImGui::PopID();

        return res;
    }

    // ----------------------------------------------------------------------------------------------

    static void AppendMarkerRecordItem(const char* prefix, int uid, const markerRecord_t &markerRecord) {

        ImGui::PushID(uid);
        ImGui::AlignTextToFramePadding();
        bool nodeOpen = ImGui::TreeNode("Object", "%s %u", prefix, uid);
        ImGui::NextColumn();
        ImGui::AlignTextToFramePadding();

        // no text here

        ImGui::NextColumn();

        if (nodeOpen) {
            static const std::vector<std::string> itemName{"X", "Y", "Z", "R", "energy", "size"};
            for (int i=0; i<4; i++) {
                ImGui::PushID(i);
                ImGui::AlignTextToFramePadding();
                ImGui::TreeNodeEx(itemName[i].c_str(), ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_Bullet);
                ImGui::NextColumn();
                    
                ImGui::Text("%.3f", markerRecord.loc(uid, i));

                ImGui::NextColumn();
                ImGui::PopID();
            }

            ImGui::Separator(); ///////////////////////

            for (int i=4; i<=5; i++) {
                ImGui::PushID(i);
                ImGui::AlignTextToFramePadding();
                ImGui::TreeNodeEx(itemName[i].c_str(), ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_Bullet);
                ImGui::NextColumn();
                if (i == 4) {
                    ImGui::Text("%.4f", markerRecord.energy(uid));
                } else {
                    ImGui::Text("%d", markerRecord.size(uid));
                }
                ImGui::NextColumn();
                ImGui::PopID();
            }

            ImGui::Separator(); ///////////////////////
            ImGui::TreePop();
        }

        ImGui::PopID();
    }
};

}  // anonymous namespace

//////////////////////
// Function Decleration
// static void DrawWindowGraphics(bool* p_open);

void GUI::post_resize(int w, int h) {

    const double dpiScale = hidpi_scaling();
    windowWidth = w / dpiScale;
    windowHeight = h / dpiScale;
}


bool GUI::MouseDownCallback(igl::opengl::glfw::Viewer &viewer, int button, int modifier) {

    if (cropActive) {

        Eigen::Vector2f mouse;
        mouse << viewer.down_mouse_x, viewer.down_mouse_y;
        logger().info("Mouse down call back");
        CropImage(mouse, MOUSEDOWN);

        // disable ligigl default mouse_down
        return true;
    }

    if (rejectActive) {

        MouseRejectCluster();

        // do not block default mouse_down
    }

    return false;
}


bool GUI::MouseUpCallback(igl::opengl::glfw::Viewer &viewer, int button, int modifier) {

    if (cropActive) {

        Eigen::Vector2f mouse;
        mouse << viewer.down_mouse_x, viewer.down_mouse_y;
        CropImage(mouse, MOUSEUP);

        // disable ligigl default mouse_up
        return true;
    }

    return false;
}


bool GUI::MouseMoveCallback(igl::opengl::glfw::Viewer &viewer, int mouse_x, int mouse_y) {

    if (cropActive) {

        Eigen::Vector2f mouse;
        mouse << mouse_x, mouse_y;
        CropImage(mouse, MOUSEMOVE);

        // disable ligigl default mouse_move
        return true;
    }

    if (rejectActive) {

        Eigen::Vector2f mouse;
        mouse << mouse_x, mouse_y;
        MouseSelectCluster(mouse);

        // do not block default mouse_move
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////////////
/// This is the main starting point
/// We override the libigl function "draw_menu" as the main GUI function

void GUI::draw_menu() {

    Draw3DImage();
    DrawMainMenuBar();
    DrawZebrafishPanel();

    if (show_log) DrawWindowLog();
    if (show_3DImage_viewer) DrawWindow3DImageViewer();
    if (show_property_editor) DrawWindowPropertyEditor();
    if (show_graphics) DrawWindowGraphics();
}


////////////////////////////////////////////////////////////////////////////////////////
// Draw 3D image


void GUI::Draw3DImage() {

    viewer.data().clear();

    // do not draw if there is no image or this is not needed
    if (currentLoadedFrames == 0) return;
    if (!showBackgroundImage) return;

    static const Eigen::Vector3d ambient = Eigen::Vector3d(146./255., 172./255., 178./255.);
    static const Eigen::Vector3d diffuse = Eigen::Vector3d(146./255., 172./255., 178./255.);
    static const Eigen::Vector3d specular = Eigen::Vector3d(0., 0., 0.);
    static Eigen::MatrixXd UV(4, 2);
    static int layerBegin_cache = -1, layerEnd_cache = -1;
    static int imgCols_cache = -1, imgRows_cache = -1;
    UV << 0, 1, 1, 1, 1, 0, 0, 0;
    V << 0, 0, 0, imgCols, 0, 0, imgCols, imgRows, 0, 0, imgRows, 0;
    F << 0, 1, 2, 2, 3, 0;

    if (imageViewerType == 0) {
        // compressed view

        if (layerBegin_cache != layerBegin || layerEnd_cache != layerEnd) {
            // if depth crop has updated
            ComputeCompressedTexture(imgData[0], 0);  // only frame 0 could reach here
            layerBegin_cache = layerBegin;
            layerEnd_cache = layerEnd;
        }
        texture = compressedImgTextureArray[frameToShow];
    } else {
        // per slice view

        image_t &imgChosen = imgData[frameToShow];
        texture = (imgChosen[sliceToShow].array() * 255).cast<unsigned char>();
        texture.transposeInPlace();
    }

    viewer.core().align_camera_center(V);
    viewer.data().set_mesh(V, F);
    viewer.data().set_uv(UV);
    viewer.data().uniform_colors(ambient, diffuse, specular);
    viewer.data().show_faces = true;
    viewer.data().show_lines = false;
    viewer.data().show_texture = true;
    viewer.data().set_texture(texture, texture, texture);
}


////////////////////////////////////////////////////////////////////////////////////////
// zebrafish panel


void GUI::DrawZebrafishPanel() {
// This panel cannot be closed

    ImGui::SetNextWindowPos(ImVec2(0.0, mainMenuHeight), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(zebrafishWidth, windowHeight-mainMenuHeight), ImGuiCond_FirstUseEver);
    ImGui::Begin("Zebrafish Config", NULL, 
        ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoMove);

    // Stage info
    ImGui::Separator();
    if (ImGui::Button("Prev Stage", ImVec2(zebrafishWidth / 2.0, 0))) {
        stage--;
        stage = std::max(1, stage);
    }
    if (ImGui::Button("Next Stage", ImVec2(zebrafishWidth / 2.0, 0))) {
        stage++;
        stage = std::min(stage, stageMax);

        if (stage == 2) stage1to2Flag = true;
        if (stage == 5) stage4to5Flag = true;
        if (stage == 6) stage5to6Flag = true;
    }
    ImGui::Text("Stage %d", stage);

    // Stage specific GUI
    ImGui::Separator();
    switch (stage) {
        case 1:
            DrawStage1();
            break;
        case 2:
            DrawStage2();
            break;
        case 3:
            DrawStage3();
            break;
        case 4:
            DrawStage4();
            break;
        case 5:
            DrawStage5();
            break;
        case 6:
            DrawStage6();
            break;
        case 7:
            DrawStage7();
            break;
        case 8:
            DrawStage8();
            break;
        default:
            assert(false);
    }

    ImGui::End();
}


////////////////////////////////////////////////////////////////////////////////////////
// main menu


void GUI::DrawMainMenuBar() {
// [ File ] [ Window ]
// This is the main menu

    if (ImGui::BeginMainMenuBar()) {

        if (ImGui::BeginMenu("File")) {
            DrawMenuFile();
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Window")) {
            DrawMenuWindow();
            ImGui::EndMenu();
        }
        mainMenuHeight = ImGui::GetWindowHeight();
        ImGui::EndMainMenuBar();
    }
}


void GUI::DrawMenuFile() {
// [ New ] [ Open ] [Close]
// Accessed from [ Main menu - File ]

    ImGui::MenuItem("New", NULL, false, false);
    if (ImGui::MenuItem("Open", "Ctrl+O")) { 
        // Only accept tif/tiff files
        std::string filename = FileDialog::openFileName("./.*", {"*.tif", "*.tiff"});
        if (!filename.empty()) {
            imagePath = filename;
            GetDescription(imagePath, layerPerImg, channelPerSlice, ttlFrames);
            if (ReadTifFirstFrame(filename, layerPerImg, channelPerSlice, imgData[0])) {
                imgRows = imgData[0][0].rows();
                imgCols = imgData[0][0].cols();
                currentLoadedFrames = 1;
                // In case the tiff image is very small
                layerPerImg = imgData[0].size();
                layerEnd = layerPerImg - 1;
            } else {
                logger().error("Error open tiff image");
                std::cerr << "Error open tiff image" << std::endl;
            }
        } 
    }

    ImGui::Separator();

    ImGui::MenuItem("Close", NULL, false, false);
}


void GUI::DrawMenuWindow() {
// [ Graphics ]
// Accessed from [ Main menu - Window ]

    ImGui::MenuItem("Log", NULL, &show_log);
    ImGui::MenuItem("3D Image Viewer", NULL, &show_3DImage_viewer);
    ImGui::MenuItem("Property Editor", NULL, &show_property_editor);

    ImGui::Separator();

    ImGui::MenuItem("Graphics", NULL, &show_graphics);
}


////////////////////////////////////////////////////////////////////////////////////////
// window: log


void GUI::DrawWindowLog() {

    ImGui::SetNextWindowPos(ImVec2(zebrafishWidth, windowHeight-logHeight), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(windowWidth-zebrafishWidth-RHSPanelWidth, logHeight), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Log", &show_log)) {
        ImGui::End();
        return;
    }

    ImGui::Separator();
    ImGui::BeginChild("scrolling", ImVec2(0, 0), false, ImGuiWindowFlags_HorizontalScrollbar);

    std::string log = oss.str();

    ImGui::TextUnformatted(log.c_str());
    ImGui::SetScrollHere(1.0f);

    ImGui::EndChild();
    ImGui::End();
}


////////////////////////////////////////////////////////////////////////////////////////
// window: 3D image viewer


void GUI::DrawWindow3DImageViewer() {

    ImGui::SetNextWindowPos(ImVec2(windowWidth-RHSPanelWidth, windowHeight-Image3DViewerHeight), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(RHSPanelWidth, Image3DViewerHeight), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("3D Image Viewer", &show_3DImage_viewer)) {
        ImGui::End();
        return;
    }

    // Plot "imgData"
    if (currentLoadedFrames > 0) {

        ImGui::PushItemWidth(RHSPanelWidth/2.0);
        std::vector<std::string> typeName{"Compressed", "Per Slice"};
        ImGui::Combo("3D Image Viewer Type", &imageViewerType, typeName);
        ImGui::PopItemWidth();

        ImGui::Separator(); ////////////////////////

        if (imageViewerType == 0) {
            // compressed viewer
            ImGui::SliderInt("Frame", &frameToShow, 0, currentLoadedFrames-1);
        } else {
            // per slice view
            ImGui::SliderInt("Frame", &frameToShow, 0, currentLoadedFrames-1);
            ImGui::SliderInt("Slice", &sliceToShow, layerBegin, layerEnd);
        }

        ImGui::Separator(); ////////////////////////

        ImGui::Text("Image path = %s", imagePath.c_str());
        ImGui::Text("Current loaded frames = %d", currentLoadedFrames);
        ImGui::Text("Using slices (top-down index) %d to %d", layerBegin, layerEnd);
        ImGui::Text("layers per image = %d", layerPerImg);
        ImGui::Text("channels per slice = %d", channelPerSlice);

    } else {

        ImGui::Text("No 3D image registered.");
    }

    ImGui::End();
}


////////////////////////////////////////////////////////////////////////////////////////
// window: property editor


void GUI::DrawWindowPropertyEditor() {

    ImGui::SetNextWindowPos(ImVec2(windowWidth-RHSPanelWidth, mainMenuHeight), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(RHSPanelWidth, windowHeight-Image3DViewerHeight-mainMenuHeight), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Property Editor", &show_3DImage_viewer)) {
        ImGui::End();
        return;
    }

    ImGui::PushItemWidth(RHSPanelWidth/2.0);
    std::vector<std::string> typeName{"Grid Search & Opt", "Clusters", "Markers"};
    ImGui::Combo("Property List Type", &propertyListType, typeName);
    ImGui::PopItemWidth();

    ImGui::Separator();
    ImGui::BeginChild("scrolling", ImVec2(0, 0), false, ImGuiWindowFlags_HorizontalScrollbar);

    switch (propertyListType) {
    case 0:
        // Grid Search & Optimization
        if (pointRecord.num == 0) {
            ImGui::Text("Optimization cylinder list is empty");
        } else {
            ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(2, 2));
            ImGui::Columns(2);

            const int maxNumItemDisplayed = 500;
            const int ttlItem = pointRecord.num;
            const int numItemToDisplay = std::min(maxNumItemDisplayed, ttlItem);
            for (int i=0; i<numItemToDisplay; i++) {

                PropertyEditorItem::AppendPointRecordItem("Point", i, pointRecord);
            }

            ImGui::Columns(1);
            if (ttlItem >= maxNumItemDisplayed) {
                ImGui::Text("Only the first %d items will be displayed", maxNumItemDisplayed);
            }
            ImGui::PopStyleVar();
        }
        break;

    case 1:
        // Clustered Cylinders
        if (clusterRecord.num == 0) {
            ImGui::Text("Cluster cylinder list is empty");
        } else {
            ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(2, 2));
            ImGui::Columns(2);

            const int maxNumItemDisplayed = 500;
            const int ttlItem = clusterRecord.num;
            const int numItemToDisplay = std::min(maxNumItemDisplayed, ttlItem);
            for (int i=0; i<numItemToDisplay; i++) {

                if (PropertyEditorItem::AppendClusterRecordItem("Cluster", i, clusterRecord)) {
                    UpdateClusterPointLoc();
                }
            }

            ImGui::Columns(1);
            if (ttlItem >= maxNumItemDisplayed) {
                ImGui::Text("Only the first %d items will be displayed", maxNumItemDisplayed);
            }
            ImGui::PopStyleVar();
        }
        break;

    case 2:
        // markers (finalized clusters)
        if (markerArray.empty()) {
            ImGui::Text("Marker cluster list is empty");
        } else {
            ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(2, 2));
            ImGui::Columns(2);

            const int maxNumItemDisplayed = 500;
            const int ttlItem = markerArray[frameToShow].num;
            const int numItemToDisplay = std::min(maxNumItemDisplayed, ttlItem);
            for (int i=0; i<numItemToDisplay; i++) {

                PropertyEditorItem::AppendMarkerRecordItem("Marker", i, markerArray[frameToShow]);
            }

            ImGui::Columns(1);
            if (ttlItem >= maxNumItemDisplayed) {
                ImGui::Text("Only the first %d items will be displayed", maxNumItemDisplayed);
            }
            ImGui::PopStyleVar();
        }
        break;

    default:
        assert(false);
        break;
    }

    ImGui::EndChild();
    ImGui::End();
}


////////////////////////////////////////////////////////////////////////////////////////
// window: graphics


void GUI::DrawWindowGraphics() {

    if (!ImGui::Begin("Graphics", &show_graphics)) {
        ImGui::End();
        return;
    }
    igl::opengl::glfw::imgui::ImGuiMenu::draw_viewer_menu();
    ImGui::End();
}


////////////////////////////////////////////////////////////////////////////////////////
// shared


void GUI::ComputeCompressedTexture(const image_t &img_, int index) {
/// Compress "img_" and store the result to "compressedImgTextureArray [index] "

    const int num = img_.size();
    assert(num > 0);
    assert(layerBegin >= 0 && layerBegin < num);
    assert(layerEnd >=0 && layerEnd < num);
    assert(layerBegin <= layerEnd);
    const int imgRows_ = img_[0].rows();
    const int imgCols_ = img_[0].cols();

    Eigen::MatrixXd compressed;
    compressed = Eigen::MatrixXd::Zero(imgRows_, imgCols_);
    for (int i=layerBegin; i<=layerEnd; i++) {
        compressed += img_[i];
    }

    compressedImgTextureArray[index] = (compressed.array() * (255.0 / double(layerEnd-layerBegin+1))).cast<unsigned char>();
    compressedImgTextureArray[index].transposeInPlace();

    logger().info("Compressed image texture (index = {}) re-computed: slice index {} to {}", index, layerBegin, layerEnd);
}


void GUI::NormalizeImage(image_t &image, double thres) {
/// This function modifies "image"

    // normalize & trim all layers
    for (auto it=image.begin(); it!=image.end(); it++) {
        Eigen::MatrixXd &slice = *it;
        for (int r=0; r<slice.rows(); r++)
            for (int c=0; c<slice.cols(); c++) {
                slice(r, c) = (slice(r, c)>=thres) ? 1.0f : slice(r, c)/thres;
            }
    }
}


////////////////////////////////////////////////////////////////////////////////////////
// maintenance methods


GUI::GUI() : pointRecord(), clusterRecord() {

    // shared
    bsplineArray.resize(1);
    imgData.resize(1);
    stage = 1;
    histBars = 50;
    showBackgroundImage = true;

    // image (imageData)
    layerPerImg = 40;  // a random guess to preview the image file
    channelPerSlice = 2;  // a random guess to preview the image file
    channelToLoad = 0;
    layerBegin = 0;
    layerEnd = layerPerImg - 1;
    resolutionX = 0;
    resolutionY = 0;
    resolutionZ = 0;
    normalizeQuantile = 0.995;
    imgHist.hist = Eigen::MatrixXf::Zero(histBars, 1);

    // grid search
    gapX_grid = 1.0;
    gapY_grid = 1.0;
    gapZ_grid = 1.0;
    rArrayMin_grid = 3.0;
    rArrayMax_grid = 6.0;
    rArrayGap_grid = 1.0;
    showPromisingPoints = true;
    gridEnergyThres = -0.1;
    promisingPointLoc.resize(1, 3);
    gridEnergyHist.hist = Eigen::MatrixXf::Zero(histBars, 1);

    // optimization
    showOptimizedPoints = true;
    optimEnergyThres = -0.1;
    optimEpsilon = 1e-4;
    optimMaxIt = 50;
    optimPointLoc.resize(1, 3);

    // cylinder filter
    cylinderEnergyThres = -0.15;
    cylinderRadiusThres = 5.0;
    cylinderIterThres = optimMaxIt;
    showCylFilterPoints = true;
    cylPointLoc.resize(1, 3);
    cylEnergyHist.hist = Eigen::MatrixXf::Zero(histBars, 1);
    cylRadiusHist.hist = Eigen::MatrixXf::Zero(histBars, 1);
    cylIterHist.hist = Eigen::MatrixXf::Zero(histBars, 1);

    // cluster filter
    clusterDistThres = 0.01;
    finalizeClusterDistThres = 2.0;
    clusterSizeThres = 10;
    showClusterFilterPoints = false;
    clusterPointLoc.resize(1, 3);
    clusterSizeHist.hist = Eigen::MatrixXf::Zero(histBars, 1);

    // ICP
    showMarkerPoints = true;
    showReferencePoints = true;
    markerPointLoc.resize(1, 3);
    refPointLoc.resize(1, 3);

    // Optical Flow
    desiredFrames = 0;

    // 3D image viewer
    V.resize(4, 3);
    F.resize(2, 3);
    imageViewerType = 0;

    // crop image
    cropActive = false;
    downClicked = false;
    showCropArea = true;
    baseLoc << 0.0f, 0.0f, 0.0f;
    currentLoc << 0.0f, 0.0f, 0.0f;
    r0 = -1;
    c0 = -1; 
    r1 = -1;
    c1 = -1;

    // manually reject clusters
    rejectActive = false;
    rejectHit = false;
    rejectMode = REJECT_AREA;
    rejectHitIndex.resize(1, 1);
    mousePickDistSquareThres = 2.0 * 2.0;  // 2 pixels by default

    // property editor
    propertyListType = 0;

    //////////////////////////////////////////////////
    // visualization
    compressedImgTextureArray.resize(1);
    sliceToShow = 0;
    frameToShow = 0;
    currentLoadedFrames = 0;
    windowWidth = 1600;
    windowHeight = 900;
    zebrafishWidth = 300;
    logHeight = 150;
    Image3DViewerHeight = 320;
    RHSPanelWidth = 300;

    // bool flag indicating whether the panel is being rendered
    show_log = false;
    show_3DImage_viewer = false;
    show_property_editor = false;
    show_graphics = false;

    // bool flag indicating moving from a stage to another
    stage1to2Flag = false;
    stage4to5Flag = false;
    stage5to6Flag = false;
}


void GUI::init(std::string imagePath_) {

    // Debug purpose
    if (!imagePath_.empty()) {
        // only true in debug mode
        imagePath = imagePath_;
        GetDescription(imagePath_, layerPerImg, channelPerSlice, ttlFrames);
        ReadTifFirstFrame(imagePath_, layerPerImg, channelPerSlice, imgData[0]);
        imgRows = imgData[0][0].rows();
        imgCols = imgData[0][0].cols();
        currentLoadedFrames = 1;
        // In case the tiff image is very small
        layerPerImg = imgData[0].size();
        layerEnd = layerPerImg - 1;

        // debug helper
        show_log = true;
        show_3DImage_viewer = true;
        show_property_editor = true;
        layerBegin = 15;
        layerEnd = 60;
        r0 = 419;
        c0 = 516;
        r1 = 499;
        c1 = 589;
    }

    // callback
    viewer.callback_mouse_down = [this](igl::opengl::glfw::Viewer &viewer, int button, int modifier) {
        return this->MouseDownCallback(viewer, button, modifier);
    };
    viewer.callback_mouse_up   = [this](igl::opengl::glfw::Viewer &viewer, int button, int modifier) {
        return this->MouseUpCallback(viewer, button, modifier);
    };
    viewer.callback_mouse_move = [this](igl::opengl::glfw::Viewer &viewer, int mouse_x, int mouse_y) {
        return this->MouseMoveCallback(viewer, mouse_x, mouse_y);
    };

    // libigl viewer
    viewer.core().orthographic = true;
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_NO_ROTATION);
    viewer.core().background_color << 0.7f, 0.7f, 0.75f, 1.0f;
        // viewer.core().is_animating = true;
    viewer.plugins.push_back(this);
    viewer.launch(true, false, "Zebrafish GUI", windowWidth, windowHeight);
}

}  // namespace zebrafish
