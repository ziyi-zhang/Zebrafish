#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>


namespace zebrafish {

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 6: Iterative Closest Point

void GUI::DrawStage6() {

    if (stage5to6Flag) {

        FinalizeClusterLoc();
        stage5to6Flag = false;
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Iterative Closest Point", ImGuiTreeNodeFlags_DefaultOpen)) {

        
        ImGui::Text("TBD");
    }

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::Text("Stage 6: ICP");
}


////////////////////////////////////////////////////////////////////////////////////////
// 

}  // namespace zebrafish
