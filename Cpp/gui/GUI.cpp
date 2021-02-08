#include <zebrafish/GUI.h>
#include <zebrafish/TiffReader.h>
#include <zebrafish/FileDialog.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/Quantile.h>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

#include <LBFGS.h>
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
            if (ImGui::Checkbox(" invalid", &clusterRecord.alive(uid))) {
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

    if (imageCrop.cropActive) {

        Eigen::Vector2f mouse;
        mouse << viewer.down_mouse_x, viewer.down_mouse_y;
        CropImage(mouse, MOUSEDOWN, imageCrop);

        // disable ligigl default mouse_down
        return true;
    }

    if (meanCrop.cropActive) {

        Eigen::Vector2f mouse;
        mouse << viewer.down_mouse_x, viewer.down_mouse_y;
        CropImage(mouse, MOUSEDOWN, meanCrop);

        // disable ligigl default mouse_down
        return true;
    }

    if (rejectActive) {

        MouseRejectCluster();

        // do not block default mouse_down
    }

    if (markerDragActive) {

        if (markerDragHit)
            MarkerDragSelect();

        // disable ligigl default mouse_down
        return true;
    }

    return false;
}


bool GUI::MouseUpCallback(igl::opengl::glfw::Viewer &viewer, int button, int modifier) {

    if (imageCrop.cropActive) {

        Eigen::Vector2f mouse;
        mouse << viewer.down_mouse_x, viewer.down_mouse_y;  // this will not be used
        CropImage(mouse, MOUSEUP, imageCrop);

        // do not block default mouse_up
    }

    if (meanCrop.cropActive) {

        Eigen::Vector2f mouse;
        mouse << viewer.down_mouse_x, viewer.down_mouse_y;  // this will not be used
        CropImage(mouse, MOUSEUP, meanCrop);

        // do not block default mouse_up
    }

    if (markerDragActive && markerDragFocused) {

        MarkerDragSetNewLoc();

        // do not block default mouse_up
    }

    return false;
}


bool GUI::MouseMoveCallback(igl::opengl::glfw::Viewer &viewer, int mouse_x, int mouse_y) {

    if (imageCrop.cropActive) {

        Eigen::Vector2f mouse;
        mouse << mouse_x, mouse_y;
        CropImage(mouse, MOUSEMOVE, imageCrop);

        // disable ligigl default mouse_move
        return true;
    }

    if (meanCrop.cropActive) {

        Eigen::Vector2f mouse;
        mouse << mouse_x, mouse_y;
        CropImage(mouse, MOUSEMOVE, meanCrop);

        // disable ligigl default mouse_move
        return true;
    }

    if (rejectActive) {

        Eigen::Vector2f mouse;
        mouse << mouse_x, mouse_y;
        MouseSelectCluster(mouse);

        // do not block default mouse_move
    }

    if (markerDragActive) {

        Eigen::Vector2f mouse;
        mouse << mouse_x, mouse_y;
        MouseSelectMarker(mouse);

        // do not block default mouse_move
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////////////
/// This is the main starting point
/// We override the libigl function "draw_menu" as the main GUI function

void GUI::draw_menu() {

    Draw3DImage();
    DrawMarkerMesh();
    DrawMainMenuBar();
    DrawZebrafishPanel();

    if (show_log) DrawWindowLog();
    if (show_3DImage_viewer) DrawWindow3DImageViewer();
    if (show_property_editor) DrawWindowPropertyEditor();
    if (show_graphics) DrawWindowGraphics();
}


////////////////////////////////////////////////////////////////////////////////////////
// Draw 3D image & marker mesh


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
            // Note: only frame 0 in stage 1 could reach here

            const int num = imgData[0].size();
            // the user could provide false metadata
            try {
                ComputeCompressedTextureMax(imgData[0], 0);
            } catch (const std::exception &e) {
                logger().warn("Failed when drawing 3D image. This is often due to wrong metadata like incorrect number of channels or slices.");
                std::cerr << "Failed when drawing 3D image. This is often due to wrong metadata like incorrect number of channels or slices." << std::endl;
            }
            layerBegin_cache = layerBegin;
            layerEnd_cache = layerEnd;
        }
        texture = compressedImgTextureArray[frameToShow];
        if (stage == 1) {
            // special contrast enhencement for stage 1
            for (int r=0; r<texture.rows(); r++)
                for (int c=0; c<texture.cols(); c++)
                    if (texture(r, c) * stage1contrast < 255) 
                        texture(r, c) *= stage1contrast;
                    else 
                        texture(r, c) = 255;
        }
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


void GUI::DrawMarkerMesh() {

    viewer.data(meshID).clear();

    // do not draw if not needed
    if (!showMarkerMesh) return;
    if (markerMeshArray.size() == 0 || markerPointLocArray.empty()) return;
    if (stage < 6) return;  // do not visualize mesh before stage 6: mesh

    // show mesh
    static Eigen::MatrixXd tempV;
    tempV = markerPointLocArray[frameToShow];
    tempV.col(2).array() -= 0.3;  // colliding with marker points
    viewer.data(meshID).show_faces = false;
    viewer.data(meshID).show_lines = true;
    viewer.data(meshID).show_texture = false;
    // viewer.data(meshID).set_colors(Eigen::RowVector3d(0.78, 0.82, 0.83));
    viewer.data(meshID).line_width = lineWidth;
    viewer.data(meshID).shininess = 0;
    viewer.data(meshID).set_mesh(
        tempV,  // mesh V
        markerMeshArray  // mesh F
    );
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
        StateChangeReset();
    }
    if (ImGui::Button("Next Stage", ImVec2(zebrafishWidth / 2.0, 0))) {

        // stage lock
        bool lock = true;
        if (stage == 1) lock = stage1Lock;
        if (stage == 2) lock = stage2Lock;
        if (stage == 3) lock = stage3Lock;
        if (stage == 4) lock = stage4Lock;

        if (lock) {
            stage++;
            stage = std::min(stage, stageMax);

            if (stage == 2) stage1to2Flag = true;
            if (stage == 3) stage2to3Flag = true;
            if (stage == 5) stage4to5Flag = true;
            if (stage == 6) stage5to6Flag = true;
        } else {
            logger().warn("[Stage lock] Cannot proceed to the next stage. At least one mandatory step missing.");
        }
        StateChangeReset();
    }
    switch (stage) {
    case 1:
        ImGui::Text("Stage 1: Image Data");
        break;
    case 2:
        ImGui::Text("Stage 2: Pre-process & B-spline");
        break;
    case 3:
        ImGui::Text("Stage 3: Grid Search");
        break;
    case 4:
        ImGui::Text("Stage 4: Optimization");
        break;
    case 5:
        ImGui::Text("Stage 5: Filter & Cluster");
        break;
    case 6:
        ImGui::Text("Stage 6: Mesh");  // iterative closest point
        break;
    case 7:
        ImGui::Text("Stage 7: Estimate Movement");  // optical flow
        break;
    case 8:
        ImGui::Text("Stage 8: Displacement & Export");
        break;
    default:
        assert(false);
        break;
    }

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
        if (ImGui::BeginMenu("Help")) {
            DrawMenuHelp();
            ImGui::EndMenu();
        }
        mainMenuHeight = ImGui::GetWindowHeight();
        ImGui::EndMainMenuBar();
    }
}


void GUI::DrawMenuFile() {
// [ New ] [ Open ] [ Close ]
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


void GUI::DrawMenuHelp() {
// [ tooltip ]
// Accessed from [ Main menu - Help ]

    ImGui::MenuItem("Tooltip", NULL, &showTooltip);
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

    static ImGuiTextBuffer buf;
    
    // std::string log = oss.str();
    // ImGui::TextUnformatted(log.c_str());
    std::string log = oss.str();
    oss.str("");
    oss.clear();
    // AddLog(log.c_str(), buf);
    buf.appendf("%s", log.c_str());


    ImGui::TextUnformatted(buf.begin(), buf.end());

    if (ImGui::GetScrollY() >= ImGui::GetScrollMaxY())
        ImGui::SetScrollHereY(1.0f);

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

        if (currentLoadedFrames > 1) {
            // advanced visualization of markers in different frames

            ImGui::Separator(); ////////////////////////

            if (ImGui::TreeNode("Advanced multi-frame marker visualization")) {

                ImGui::Checkbox("Manual override", &manualOverrideMarkerVis);
                if (ImGui::Checkbox("Show all markers", &showAllMarkers)) {
                    if (showAllMarkers) {
                        for (int i=0; i<markerPointStatusArray.rows(); i++)
                            markerPointStatusArray(i) = true;
                    } else {
                        for (int i=0; i<markerPointStatusArray.rows(); i++)
                            markerPointStatusArray(i) = false;
                        markerPointStatusArray(frameToShow) = true;
                    }
                }

                for (int i=0; i<markerPointStatusArray.rows(); i++) {
                    std::string label = "Frame ";
                    label += std::to_string(i);
                    ImGui::Checkbox(label.c_str(), &markerPointStatusArray(i));
                }

                ImGui::TreePop();
                ImGui::Separator();
            }

            ImGui::Separator(); ////////////////////////

            if (ImGui::TreeNode("Advanced optical flow visualization")) {

                ImGui::Checkbox("Show optical flow", &showOpticalFlow);

                ImGui::TreePop();
                ImGui::Separator();
            }
        }

        ImGui::Separator(); ////////////////////////

        ImGui::Text("Image path = %s", imagePath.c_str());
        ImGui::Text("Current loaded frames = %d", currentLoadedFrames);
        ImGui::Text("Using slices (top-down index) %d to %d", layerBegin, layerEnd);
        ImGui::Text("layers per image = %d", layerPerImg);
        ImGui::Text("channels per slice = %d", channelPerSlice);
        ImGui::Text("Rows = %d  Cols = %d", imgRows, imgCols);

        ImGui::Separator(); ////////////////////////

        if (ImGui::TreeNode("Advanced viewer")) {

            ImGui::PushItemWidth(RHSPanelWidth/3.0);
            std::vector<std::string> typeName{"Max", "Mean"};
            if (ImGui::Combo("Compress (flatten) method", &imageViewerCompressType, typeName)) {
                ComputeCompressedTextureForAllLoadedFrames();
            }

            if (imageViewerCompressType == COMPRESS_AVG) {
                ImGui::SliderFloat("Darken factor", &imageViewerDarkenFactor_avg, 1.0, 3.0);
            } else if (imageViewerCompressType == COMPRESS_MAX) {
                ImGui::SliderFloat("Darken factor", &imageViewerDarkenFactor_max, 1.0, 3.0);
            }
            ImGui::PopItemWidth();
            ImGui::TreePop();
            ImGui::Separator();
        }

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

                PropertyEditorItem::AppendPointRecordItem("Cylinder", i, pointRecord);
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


void GUI::ComputeCompressedTextureAvg(const image_t &img_, int index) {
/// Compress "img_" and store the result to "compressedImgTextureArray[index]"
/// Flatten by taking the average of all slices

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

    compressedImgTextureArray[index] = (compressed.array() * (255.0 / double(layerEnd-layerBegin+1) / imageViewerDarkenFactor_avg)).cast<unsigned char>();
    compressedImgTextureArray[index].transposeInPlace();

    logger().info("Compressed (avg) image texture (index = {}) re-computed: slice index {} to {}", index, layerBegin, layerEnd);
}


void GUI::ComputeCompressedTextureMax(const image_t &img_, int index) {
/// Compress "img_" and store the result to "compressedImgTextureArray[index]"
/// Flatten by taking the max of all slices

    const int num = img_.size();
    /*
    assert(num > 0);
    assert(layerBegin >= 0 && layerBegin < num);
    assert(layerEnd >=0 && layerEnd < num);
    assert(layerBegin <= layerEnd);
    */
    if (!(num > 0)) {
        std::cerr << "ERROR: assert(num > 0)" << std::endl;
        return;
    }
    if (!(layerBegin >= 0 && layerBegin < num)) {
        std::cerr << "ERROR: assert(layerBegin >= 0 && layerBegin < num)" << std::endl;
        return;
    }
    if (!(layerEnd >=0 && layerEnd < num)) {
        std::cerr << "ERROR: assert(layerEnd >=0 && layerEnd < num)" << std::endl;
        return;
    }
    const int imgRows_ = img_[0].rows();
    const int imgCols_ = img_[0].cols();

    Eigen::MatrixXd compressed;
    compressed = Eigen::MatrixXd::Zero(imgRows_, imgCols_);
    for (int i=layerBegin; i<=layerEnd; i++) {
        compressed = compressed.cwiseMax(img_[i]);
    }

    compressedImgTextureArray[index] = (compressed.array() * 255.0 / imageViewerDarkenFactor_max).cast<unsigned char>();
    compressedImgTextureArray[index].transposeInPlace();

    logger().info("Compressed (max) image texture (index = {}) re-computed: slice index {} to {}", index, layerBegin, layerEnd);
}


void GUI::ComputeCompressedTextureForAllLoadedFrames() {

    for (int i=0; i<currentLoadedFrames; i++) {
        
        switch (imageViewerCompressType) {
            case COMPRESS_AVG:
                ComputeCompressedTextureAvg(imgData[i], i);
                break;
            case COMPRESS_MAX:
                ComputeCompressedTextureMax(imgData[i], i);
                break;
            default:
                assert(false);
                break;
        }
    }
}


void GUI::UpdateMarkerPointLocArray() {
// Critical marker location visualization array

    markerPointLocArray.resize(currentLoadedFrames);
    // only show the markers in one frame
    markerPointStatusArray.resize(currentLoadedFrames, 1);
    for (int i=0; i<markerPointStatusArray.rows(); i++)
        markerPointStatusArray(i) = false;
    markerPointStatusArray(frameToShow) = true;

    for (int i=0; i<currentLoadedFrames; i++) {

        markerRecord_t &markerRecord = markerArray[i];
        const int N = markerRecord.num;
        Eigen::MatrixXd tempLoc;

        markerPointLocArray[i].resize(N, 3);
        
        tempLoc.resize(N, 3);
        tempLoc.col(0) = markerRecord.loc.col(0);  // x
        tempLoc.col(1) = markerRecord.loc.col(1);  // y
        tempLoc.col(2) = markerRecord.loc.col(2);  // z

        markerPointLocArray[i].col(0) = tempLoc.col(1).array() + 0.5;
        markerPointLocArray[i].col(1) = (imgRows-0.5) - tempLoc.col(0).array();
        markerPointLocArray[i].col(2) = tempLoc.col(2);
    }

    logger().trace("   [Vis & Critical] Marker point location array updated: frames = {}", currentLoadedFrames);
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


void GUI::StateChangeReset() {

    MarkerDragReset();
}


////////////////////////////////////////////////////////////////////////////////////////
// maintenance methods

// definition of "cylinder" class static variables
double cylinder::alpha;
double cylinder::K;
double cylinder::H;

GUI::GUI() : pointRecord(), clusterRecord() {

    // shared
    bsplineArray.resize(1);
    imgData.resize(1);
    cylinder::alpha = 0.5;
    cylinder::K = std::sqrt(2);
    cylinder::H = 2.5;
    stage = 1;
    histBars = 50;
    showBackgroundImage = true;
    showTooltip = true;
    lineWidth = 4.0;
    reverseColor = false;

    // image (imageData)
    layerPerImg = 1;  // a random guess to preview the image file
    channelPerSlice = 1;  // a random guess to preview the image file
    ttlFrames = 1;
    channelToLoad = 0;
    layerBegin = 0;
    layerEnd = layerPerImg - 1;
    resolutionX = 0;
    resolutionY = 0;
    resolutionZ = 0;
    normalizeQuantile = 0.995;
    stage1contrast = 1.0;
    imgHist.hist = Eigen::MatrixXf::Zero(histBars, 1);

    // membrane mask
    membraneMask.resize(1);
    membraneMaskLoad = false;
    membraneMaskCylApply = false;
    membraneMaskClusterApply = false;
    maskMax = 0.0;
    maskThres = 0.3;

    // B-spline
    bsplineDegree = 2;
    bsplineSolverTol = 1e-7;
    bsplineArray[0].Set_degree(bsplineDegree);
    bsplineArray[0].Set_solverTol(bsplineSolverTol);

    // grid search
    gapX_grid = 1.0;
    gapY_grid = 1.0;
    gapZ_grid = 1.0;
    rArrayMin_grid = 3.0;
    rArrayMax_grid = 6.0;
    rArrayGap_grid = 1.0;
    showPromisingPoints = true;
    skipMembrane = false;
    gridEnergyThres = -0.1;
    promisingPointLoc.resize(0, 3);
    gridEnergyHist.hist = Eigen::MatrixXf::Zero(histBars, 1);

    // optimization
    showOptimizedPoints = true;
    optimEnergyThres = -0.1;
    optimEpsilon = 1e-4;
    optimMaxIt = 50;
    optimPointLoc.resize(0, 3);

    // cylinder filter
    cylinderEnergyThres = -0.1;
    cylinderRadiusThres = 6.0;
    cylinderIterThres = optimMaxIt;
    showCylFilterPoints = true;
    cylFilterMembraneCheck = true;
    cylPointLoc.resize(0, 3);
    cylEnergyHist.hist = Eigen::MatrixXf::Zero(histBars, 1);
    cylRadiusHist.hist = Eigen::MatrixXf::Zero(histBars, 1);
    cylIterHist.hist = Eigen::MatrixXf::Zero(histBars, 1);

    // cluster filter
    clusterDistThres = 0.01;
    finalizeClusterDistThres = 2.0;
    clusterSizeThres = 1;
    showClusterFilterPoints = false;
    clusterPointLoc.resize(0, 3);
    clusterSizeHist.hist = Eigen::MatrixXf::Zero(histBars, 1);

    // ICP
    showMarkerPoints = true;
    showReferencePoints = true;
    showICPLines = false;
    showMarkerMesh = false;
    ICP_patternRows = 0;
    ICP_patternCols = 0;
    ICP_patternSpacing = 18.0;
    ICP_xDisp = 0.0f;
    ICP_yDisp = 0.0f;
    ICP_angleRot = 0.0f;
    ICP_scale = 1.0f;
    refPointLoc.resize(0, 3);
    refV.resize(0, 3);
    refV_aligned.resize(0, 3);
    ICP_matchIdx.resize(0, 1);
    markerMeshArray.resize(0, 3);
    ICP_Rmat = Eigen::MatrixXd::Identity(3, 3);
    ICP_Tmat = Eigen::MatrixXd::Zero(3, 1);

    // Optical Flow
    desiredFrames = 0;
    opticalFlowAlpha = 0.1;
    opticalFlowIter = 30;
    showOpticalFlow = false;

    // Displacement
    depthCorrectionNum = 30;
    depthCorrectionGap = 0.1;
    optimMaxXYDisp = 5.0;

    // 3D image viewer
    V.resize(4, 3);
    F.resize(2, 3);
    imageViewerType = 0;
    imageViewerCompressType = COMPRESS_MAX;
    imageViewerDarkenFactor_avg = 1.0;
    imageViewerDarkenFactor_max = 1.4;

    // [mouse pick] manually reject clusters
    rejectActive = false;
    rejectHit = false;
    rejectMode = REJECT_AREA;
    rejectHitIndex.resize(0, 1);
    mousePickDistSquareThres = 3.0 * 3.0;  // 3 pixels by default

    // [mouse pick] manually drag markers
    MarkerDragReset();

    // property editor
    propertyListType = 0;

    //////////////////////////////////////////////////
    // visualization
    compressedImgTextureArray.resize(1);
    markerPointLocArray.clear();
    markerPointStatusArray.resize(1, 1);
    manualOverrideMarkerVis = false;
    showAllMarkers = false;
    sliceToShow = 0;
    frameToShow = 0;
    currentLoadedFrames = 0;
    windowWidth = 1600;
    windowHeight = 900;
    zebrafishWidth = 300;
    logHeight = 150;
    Image3DViewerHeight = 320;
    RHSPanelWidth = 300;
    // color
    markerPointColor.resize(1, 3);
    markerPointColor << 0.93, 0.32, 0.15;

    // bool flag indicating whether the panel is being rendered
    show_log = true;
    show_3DImage_viewer = true;
    show_property_editor = false;
    show_graphics = false;

    // bool flag indicating moving from a stage to another
    stage1to2Flag = false;
    stage2to3Flag = false;
    stage4to5Flag = false;
    stage5to6Flag = false;

    // stage lock
    stage1Lock = false;
    stage2Lock = false;
    stage3Lock = false;
    stage4Lock = false;
}


void GUI::init(std::string imagePath_, std::string maskPath_, int debugMode) {

    // Debug purpose
    if (!imagePath_.empty()) {
        // only true in debug mode
        imagePath = imagePath_;
        maskPath = maskPath_;
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
        show_property_editor = false;
    }

    // Debug mode
    if (debugMode == -1) {
        // very small area used by "build"
        layerBegin = 24;
        layerEnd = 40;
        imageCrop.r0 = 419;
        imageCrop.c0 = 516;
        imageCrop.r1 = 469;
        imageCrop.c1 = 556;
    } else if (debugMode == 1) {
        // large area
        layerBegin = 24;
        layerEnd = 46;
        imageCrop.r0 = 356;
        imageCrop.c0 = 448;
        imageCrop.r1 = 507;
        imageCrop.c1 = 596;
    } else if (debugMode == 2) {
        // large area
        layerBegin = 24;
        layerEnd = 46;
        imageCrop.r0 = 356;
        imageCrop.c0 = 448;
        imageCrop.r1 = 437;
        imageCrop.c1 = 556;
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
    int defaultMeshID = viewer.selected_data_index;
    meshID = viewer.append_mesh();
    viewer.selected_data_index = defaultMeshID;
    viewer.plugins.push_back(this);
    viewer.launch(true, false, "Zebrafish GUI", windowWidth, windowHeight);
}

}  // namespace zebrafish
