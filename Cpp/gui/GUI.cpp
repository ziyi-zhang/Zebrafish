#include <zebrafish/GUI.h>
#include <zebrafish/TiffReader.h>
#include <zebrafish/FileDialog.h>

#include <string>
#include <fstream>


namespace zebrafish {

namespace {

}  // anonymous namespace

//////////////////////
// Function Decleration
// static void DrawWindowGraphics(bool* p_open);

void GUI::post_resize(int w, int h) {

    const double dpiScale = 2.0;  // FIXME
    windowWidth = w / dpiScale;
    windowHeight = h / dpiScale;
}


bool GUI::MouseDownCallback(igl::opengl::glfw::Viewer &viewer, int button, int modifier) {

    if (cropActive) {
        std::cout << viewer.down_mouse_x << " " << viewer.down_mouse_y << std::endl;
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
        stage = std::max(0, stage);
    }
    if (ImGui::Button("Next Stage", ImVec2(zebrafishWidth / 2.0, 0))) {
        stage++;
        stage = std::min(stage, stageMax);
    }
    ImGui::Text("Stage %d", stage);

    // Stage specific GUI
    ImGui::Separator();
    switch (stage) {
        case 0:
            DrawStage0();
            break;
        case 1:
            DrawStage1();
            break;
        case 2:
            DrawStage2();
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

    std::ifstream logFile("Zebrafish_gui.log", std::ios::binary | std::ios::ate);
    std::streamsize size = logFile.tellg();
    logFile.seekg(0, std::ios::beg);
    char logBuffer[size];

    if (logFile.read(logBuffer, size)) {
        ImGui::TextUnformatted(logBuffer);
        ImGui::SetScrollHere(1.0f);
    } else {
        ImGui::Text("Failed to load log file.");
    }

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
    std::vector<std::string> typeName{"Grid Starting Points", "Cylinders"};
    ImGui::Combo("Property List Type", &propertyListType, typeName);
    ImGui::PopItemWidth();

    ImGui::Separator();

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


GUI::GUI() : bsplineSolver() {

    stage = 0;
    slice = 0;

    // image (imageData)
    layerPerImg = 40;  // a random guess to preview the image file
    channelPerSlice = 2;  // a random guess to preview the image file
    resolutionX = 0;
    resolutionY = 0;
    resolutionZ = 0;
    normalizeQuantile = 0.995;

    // clip image
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
    viewer.core().background_color << 0.7f, 0.7f, 0.75f, 1.0f;
        // viewer.core().is_animating = true;
    viewer.plugins.push_back(this);
    viewer.launch(true, false, "Zebrafish GUI", windowWidth, windowHeight);
}

}  // namespace zebrafish
