#include <zebrafish/GUI.h>
#include <zebrafish/TiffReader.h>
#include <zebrafish/FileDialog.h>


namespace zebrafish {

// Function Decleration
// static void DrawWindowGraphics(bool* p_open);


void GUI::DrawStage0() {

    if(ImGui::Button("Load Image")) {
        
    }

    ImGui::Text("Stage 0");
    // viewer.core().set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_TRACKBALL);
}


void GUI::DrawStage1() {

    // if(ImGui::SliderInt2("Slide"))
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


    ImGui::Text("Stage 1");
}


void GUI::DrawStage2() {

    ImGui::Text("Stage 2");
}


void GUI::draw_menu() {

    if (show_graphics) DrawWindowGraphics();

    DrawMainMenuBar();

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
    }

    if (ImGui::Button("Next", ImVec2(-1, 0))) {
        if(stage == 0) {
            read_tif_image("../data/test_img.tif", img);
        }
        stage++;
        stage = std::min(stage, 2);
    }
    if (ImGui::Button("Prev", ImVec2(-1, 0))) {
        stage--;
        stage = std::max(0, stage);
    }
}


////////////////////////////////////////////////////////////////////////////////////////
// main menu


void GUI::DrawMainMenuBar() {
// File, Window

    if (ImGui::BeginMainMenuBar()) {

        if (ImGui::BeginMenu("File")) {
            DrawMenuFile();
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Window")) {
            DrawMenuWindow();
            ImGui::EndMenu();
        }
        ImGui::EndMainMenuBar();
    }
}


void GUI::DrawMenuFile() {
// New, Open
// Accessed from [ Main menu - File ]

    ImGui::MenuItem("New", NULL, false, false);
    if (ImGui::MenuItem("Open", "Ctrl+O")) { 
        std::string fname = FileDialog::openFileName("./.*", {"*.tif", "*.tiff"});
        if (!fname.empty())
        {
            std::cout << fname << std::endl;
        } 
    }
}


void GUI::DrawMenuWindow() {
// Graphics
// Accessed from [ Main menu - Window ]

    if (ImGui::MenuItem("Graphics", NULL, &show_graphics)) {
    }
}


void GUI::DrawWindowGraphics() {

    if (!ImGui::Begin("Graphics", &show_graphics)) {
        ImGui::End();
        return;
    }
    igl::opengl::glfw::imgui::ImGuiMenu::draw_viewer_menu();
    ImGui::Text("Test");
    ImGui::End();
}


////////////////////////////////////////////////////////////////////////////////////////
// maintenance methods


GUI::GUI() : stage(0), slice(0) {

}


void GUI::init() {

    viewer.core().background_color << 0.7f, 0.7f, 0.75f, 1.0f;
    // viewer.core().is_animating = true;
    viewer.plugins.push_back(this);
    viewer.launch(true, false, "Zebrafish GUI");
}

}  // namespace zebrafish
