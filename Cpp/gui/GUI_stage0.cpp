#include <zebrafish/GUI.h>
#include <zebrafish/TiffReader.h>
#include <zebrafish/FileDialog.h>


namespace zebrafish {

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 0: image read

void GUI::DrawStage0() {

    if(ImGui::Button("Load Image")) {
        
    }

    ImGui::Text("Stage 0");
    // viewer.core().set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_TRACKBALL);
}

}  // namespace zebrafish
