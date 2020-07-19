#include <zebrafish/GUI.h>
#include <zebrafish/TiffReader.h>
#include <zebrafish/FileDialog.h>
#include <zebrafish/Logger.hpp>

#include <igl/unproject_onto_mesh.h>
#include <string>
#include <sstream>
#include <algorithm>


namespace zebrafish {

namespace {

struct PropertyEditorItem {

    static void AppendPointRecordItem(const char* prefix, int uid, const pointRecord_t &pointRecord) {

        ImGui::PushID(uid);
        ImGui::AlignTextToFramePadding();
        bool nodeOpen = ImGui::TreeNode("Object", "%s %u", prefix, uid);
        ImGui::NextColumn();
        ImGui::AlignTextToFramePadding();
        static bool enabled = true; //
        if (pointRecord.optimization(0, 0) == 0)
            ImGui::Checkbox("to be optimized", &enabled);
        else
            ImGui::Text("has been optimized");

        ImGui::NextColumn();

        if (nodeOpen) {
            static const std::vector<std::string> itemName{"Energy", "x", "y", "z", "h"};
            for (int i=0; i<5; i++) {
                ImGui::PushID(i);
                ImGui::AlignTextToFramePadding();
                ImGui::TreeNodeEx(itemName[i].c_str(), ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_Bullet);
                ImGui::NextColumn();
                if (i == 0) {
                    ImGui::Text("%.4f", pointRecord.grid_search(uid, 4));
                } else {
                    ImGui::Text("%.2f", pointRecord.grid_search(uid, i-1));
                }
                ImGui::NextColumn();
                ImGui::PopID();
            }
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

        static int fid;
        static Eigen::Vector3f bc, loc;
        static bool hit;
        hit = igl::unproject_onto_mesh(
            Eigen::Vector2f(viewer.down_mouse_x, viewer.down_mouse_y), 
            viewer.core().view, 
            viewer.core().proj,
            viewer.core().viewport, 
            V, 
            F, 
            fid, 
            bc);

        logger().debug("Mouse = [{} x {}]", viewer.down_mouse_x, viewer.down_mouse_y);
        if (hit) {
            loc(0) = (V(F(fid, 0), 0) * bc(0) + V(F(fid, 1), 0) * bc(1) + V(F(fid, 2), 0) * bc(2));
            loc(1) = (V(F(fid, 0), 1) * bc(0) + V(F(fid, 1), 1) * bc(1) + V(F(fid, 2), 1) * bc(2));
            loc(2) = (V(F(fid, 0), 2) * bc(0) + V(F(fid, 1), 2) * bc(1) + V(F(fid, 2), 2) * bc(2));
            logger().debug("fid = {} | bc = {:.2f} x {:.2f} x {:.2f} | loc = {:.2f} x {:.2f} x {:.2f}", fid, bc(0), bc(1), bc(2), loc(0), loc(1), loc(2));
        } else {
            logger().debug("No hit on image");
        }

        
        // disable ligigl default mouse_down
        return true;
    }
    return false;
}


////////////////////////////////////////////////////////////////////////////////////////
/// This is the main starting point
/// We override the libigl function "draw_menu" as the main GUI function

void GUI::draw_menu() {

    DrawMainMenuBar();
    DrawZebrafishPanel();

    if (show_log) DrawWindowLog();
    if (show_3DImage_viewer) DrawWindow3DImageViewer();
    if (show_property_editor) DrawWindowPropertyEditor();
    if (show_graphics) DrawWindowGraphics();
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
            if (ReadTifFirstImg(filename, layerPerImg, channelPerSlice, img)) {
                // In case the tiff image is very small
                layerPerImg = img.size();
            } else {
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

    // Plot "img"
    if (!img.empty()) {

        ImGui::Text("Not implemented yet");
        /*
        ImGui::Image()
        viewer.data().clear();
        viewer.core().set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_NO_ROTATION);

        texture = (img[slice].array()).cast<unsigned char>();

        int xMax = img[slice].cols();
        int yMax = img[slice].rows();
        Eigen::MatrixXd V(4, 3);
        V << 0, 0, 0, yMax, 0, 0, yMax, xMax, 0, 0, xMax, 0;
        viewer.core().align_camera_center(V);

        Eigen::MatrixXi F(2, 3);
        F << 0, 1, 2, 2, 3, 0;

        Eigen::MatrixXd UV(4, 2);
        UV << 0, 1, 1, 1, 1, 0, 0, 0;
        viewer.data().set_mesh(V, F);
        viewer.data().set_uv(UV);
        viewer.data().show_faces = true;
        viewer.data().show_lines = false;
        viewer.data().show_texture = true;
        viewer.data().set_texture(texture, texture, texture);
        */
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
    std::vector<std::string> typeName{"Grid Search & Opt", "Cylinders"};
    ImGui::Combo("Property List Type", &propertyListType, typeName);
    ImGui::PopItemWidth();

    ImGui::Separator();
    ImGui::BeginChild("scrolling", ImVec2(0, 0), false, ImGuiWindowFlags_HorizontalScrollbar);

    switch (propertyListType) {
    case 0:
        // Grid Search
        if (pointRecord.num == 0) {
            ImGui::Text("Optimization point list is empty");
        } else {
            ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(2, 2));
            ImGui::Columns(2);

            const int maxNumItemDisplayed = 300;
            const int ttlItem = pointRecord.num;
            const int numItemToDisplay = std::min(maxNumItemDisplayed, ttlItem);
            for (int i=0; i<numItemToDisplay; i++) {

                PropertyEditorItem::AppendPointRecordItem("Point", i, pointRecord);
            }

            ImGui::Columns(1);
            if (ttlItem >= maxNumItemDisplayed) {
                ImGui::Text("Only the first %d items will be displayed", maxNumItemDisplayed);
            }
            ImGui::PopID();
        }
        break;

    case 1:
        // Optimized Cylinders
        ImGui::Text("Candidate cylinder list is empty");
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
// maintenance methods


GUI::GUI() : bsplineSolver(), pointRecord() {

    stage = 1;
    slice = 0;

    // core
    gridEnergyThres = -0.05;


    // image (imageData)
    layerPerImg = 40;  // a random guess to preview the image file
    channelPerSlice = 2;  // a random guess to preview the image file
    resolutionX = 0;
    resolutionY = 0;
    resolutionZ = 0;
    normalizeQuantile = 0.995;

    // texture image
    V.resize(4, 3);
    F.resize(2, 3);

    // crop image
    cropActive = false;
    clickCount = 0;
    r0 = -1;
    c0 = -1; 
    r1 = -1;
    c1 = -1;

    // property editor
    propertyListType = 0;

    //////////////////////////////////////////////////
    // visualization
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
}


void GUI::init(std::string imagePath) {

    // Debug purpose
    if (!imagePath.empty()) {
        // only true in debug mode
        ReadTifFirstImg(imagePath, layerPerImg, channelPerSlice, img);
        // In case the tiff image is very small
        layerPerImg = img.size();
    }

    // callback
    viewer.callback_mouse_down = [this](igl::opengl::glfw::Viewer &viewer, int button, int modifier) {
        return this->MouseDownCallback(viewer, button, modifier);
    };

    // libigl viewer
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_NO_ROTATION);
    viewer.core().background_color << 0.7f, 0.7f, 0.75f, 1.0f;
        // viewer.core().is_animating = true;
    viewer.plugins.push_back(this);
    viewer.launch(true, false, "Zebrafish GUI", windowWidth, windowHeight);
}

}  // namespace zebrafish
