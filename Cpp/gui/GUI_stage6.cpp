#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>


namespace zebrafish {

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 6: 

void GUI::DrawStage6() {

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Cylinder Filter", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::Text("TBD");
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Cluster Filter", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::Text("TBD");
    }

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::Text("Stage 6: ");
}


////////////////////////////////////////////////////////////////////////////////////////
// 



}  // namespace zebrafish
