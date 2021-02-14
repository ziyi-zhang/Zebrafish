#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>
#include <igl/project.h>


namespace zebrafish {

namespace {

}  // anonymous namespace

void GUI::DrawText(
    Eigen::Vector3d pos, 
    const std::string &text,
    const Eigen::Vector4f color) {
// DO NOT USE THIS

    Eigen::Vector3f coord = igl::project(Eigen::Vector3f(pos.cast<float>()),
    viewer.core().view, viewer.core().proj, viewer.core().viewport);

    // Draw text labels slightly bigger than normal text
    double pixel_ratio_ = 2;
    ImDrawList* drawList = ImGui::GetWindowDrawList();
    drawList->AddText(ImGui::GetFont(), ImGui::GetFontSize() * 1.2,
        ImVec2(coord[0]/pixel_ratio_, (viewer.core().viewport[3] - coord[1])/pixel_ratio_),
        ImGui::GetColorU32(ImVec4(
            color(0),
            color(1),
            color(2),
            color(3))),
        &text[0], &text[0] + text.size());
}

////////////////////////////////////////////////////////////////////////////////////////
// Special module about visualizing some effects

void GUI::DrawReferenceDots() {

    if (!show_refPoints) return;
    static Eigen::Matrix3d loc = (Eigen::Matrix3d() << 0, 0, 1, imgCols, imgRows, layerPerImg, imgCols-1, imgRows-1, layerPerImg-1).finished();
    static Eigen::MatrixXd referencePointColor = [] {
        Eigen::MatrixXd tmp(1, 3);
        tmp << 0.33, 0.83, 0.33;
        return tmp;
    } ();
    viewer.data().add_points(loc, referencePointColor);
}


void GUI::ShowAllMarkerIndex() {

    if (markerArray.empty()) return;
    const markerRecord_t &markerRecord = markerArray[frameToShow];
    const int N = markerRecord.num;
    if (N == 0) return;

    viewer.data().show_labels = true;
    Eigen::Vector3d p;
    for (int i=0; i<N; i++) {
        p << 0.5+markerRecord.loc(i, 1), (imgRows-0.5)-markerRecord.loc(i, 0), markerRecord.loc(i, 2) + 0.003;
        viewer.data().add_label(p, std::to_string(i));
    }
}

}  // namespace zebrafish
