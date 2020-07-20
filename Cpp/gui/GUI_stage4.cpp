#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>


namespace zebrafish {

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 4: Optimization

void GUI::DrawStage4() {

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::Text("Optimization");
    if (ImGui::Button("Start Optimization")) {
        
        Optimization();
    }

    ImGui::Text("Stage 4");
}


////////////////////////////////////////////////////////////////////////////////////////
// Optimization


void GUI::Optimization() {

    logger().info("Not implemented");
}

}  // namespace zebrafish
