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

    windowWidth = w;
    windowHeight = h;
}


////////////////////////////////////////////////////////////////////////////////////////
/// This is the main starting point
/// We override the libigl function "draw_menu" as the main GUI function

void GUI::draw_menu() {

    DrawMainMenuBar();
    DrawZebrafishPanel();

    if (show_log) DrawWindowLog();
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
    ImGui::MenuItem("Graphics", NULL, &show_graphics);
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
// window: log


void GUI::DrawWindowLog() {

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
    } else {
        ImGui::Text("Failed to load log file.");
    }

    ImGui::EndChild();
    ImGui::End();
}


////////////////////////////////////////////////////////////////////////////////////////
// maintenance methods


GUI::GUI() {

    stage = 0;
    slice = 0;

    layerPerImg = 40;  // a random guess to preview the image file
    channelPerSlice = 2;  // a random guess to preview the image file
}


void GUI::init(std::string imagePath) {

    // Debug purpose
    if (!imagePath.empty()) {
        // only true in debug mode
        ReadTifFirstImg(imagePath, layerPerImg, channelPerSlice, img);
        // In case the tiff image is very small
        layerPerImg = img.size();
    }

    // libigl viewer
    viewer.core().background_color << 0.7f, 0.7f, 0.75f, 1.0f;
        // viewer.core().is_animating = true;
    viewer.plugins.push_back(this);
    viewer.launch(true, false, "Zebrafish GUI", windowWidth, windowHeight);
}

}  // namespace zebrafish
