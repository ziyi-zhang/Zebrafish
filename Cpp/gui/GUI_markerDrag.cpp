#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>

#include <igl/unproject_onto_mesh.h>


namespace zebrafish {

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Special module about mouse draging a marker

void GUI::RenderMarkerDragGUI() {

    if (ImGui::TreeNode("Manually move markers")) {

        const float inputWidth = ImGui::GetWindowWidth() / 3.0;
        ImGui::PushItemWidth(inputWidth);
        ImGui::Checkbox("[Mouse] Move markers", &markerDragActive);
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Manually change the XY location of the markers using mouse.\nCheck this box and move close to a marker. Click and hold to select one marker. Drag it to a new location and release button.");
        }
        if (ImGui::Button("Optimize this frame")) {
            MarkerDragReset();
            MarkerDepthCorrection(frameToShow, 0, 0, false);
            // update visualization variable
            UpdateMarkerPointLocArray();
        }
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Optimize the XY location of markers in this frame again. This is optional.\nIf any marker has been manually moved, it is highly suggested to optimize this frame again.\nThis button does not optimize the depth.");
        }
        ImGui::PopItemWidth();

        ImGui::TreePop();
        ImGui::Separator();
    }
}


void GUI::MouseSelectMarker(const Eigen::Vector2f &mouse) {
/// called by "MouseMoveCallback" when "markerDragActive" is true and "markerDragFocused" is false

    static int fid;
    static Eigen::Vector3f bc;
    static double x, y;
    static bool hit;

    hit = igl::unproject_onto_mesh(
        Eigen::Vector2f(mouse(0), viewer.core().viewport(3)-mouse(1)), 
        viewer.core().view, 
        viewer.core().proj,
        viewer.core().viewport, 
        V_texture, 
        F_texture, 
        fid, 
        bc);

    if (hit) {
        // has hit

        y =                 (V_texture(F_texture(fid, 0), 0) * bc(0) + V_texture(F_texture(fid, 1), 0) * bc(1) + V_texture(F_texture(fid, 2), 0) * bc(2)) - 0.5;
        x = imgRows - 0.5 - (V_texture(F_texture(fid, 0), 1) * bc(0) + V_texture(F_texture(fid, 1), 1) * bc(1) + V_texture(F_texture(fid, 2), 1) * bc(2));

        markerDragLoc << x, y;

        if (!markerDragFocused) {
            // search for the nearest marker
            double dist_square, minDistSquare;
            const double distSquareThres = 3.0 * 3.0;
            minDistSquare = distSquareThres;  // reset
            markerRecord_t &markerRecord = markerArray[frameToShow];
            for (int i=0; i<markerRecord.num; i++) {

                dist_square = (markerRecord.loc(i, 0)-x)*(markerRecord.loc(i, 0)-x) + (markerRecord.loc(i, 1)-y)*(markerRecord.loc(i, 1)-y);
                if (dist_square >= mousePickDistSquareThres) continue;
                if (dist_square < minDistSquare) {
                    // find a new nearest marker

                    minDistSquare = dist_square;
                    markerDragHitIndex = i;
                }
            }
            // update markerDragHit
            markerDragHit = minDistSquare < distSquareThres;
        }
    } else {

        // update markerDragHit
        markerDragHit = false;
    }
}


void GUI::MarkerDragSelect() {

    markerDragFocused = true;
    logger().debug("Manual move marker: select index {}", markerDragHitIndex);
}


void GUI::MarkerDragSetNewLoc() {

    double x = markerDragLoc(0);
    double y = markerDragLoc(1);

    // override the location
    markerRecord_t &markerRecord = markerArray[frameToShow];
    if (markerDragHitIndex < 0 || markerDragHitIndex > markerRecord.num-1) {
        logger().error("[Fatal error] markerDragHitIndex is invalid: {}", markerDragHitIndex);
    } else {
        markerRecord.loc(markerDragHitIndex, 0) = x;
        markerRecord.loc(markerDragHitIndex, 1) = y;
        logger().debug("Manually move marker: index {} moved to location {}x{}", markerDragHitIndex, x, y);
    }

    // update markerDrag flags
    markerDragHit = false;
    markerDragFocused = false;
    markerDragHitIndex = -1;

    // update visualization
    UpdateMarkerPointLocArray();
}


void GUI::MarkerDragReset() {

    markerDragActive = false;
    markerDragHit = false;
    markerDragFocused = false;
    markerDragHitIndex = -1;
}


////////////////////////////////////////////////////////////////////////////////////////
// Visualization


void GUI::MarkerDragVisualization() {

    if (!markerDragActive) return;

    // change color to the chosen marker
    if (markerDragHit) {

        Eigen::MatrixXd selectedLoc(1, 3);
        static const double deltaZ = 0.3;
        selectedLoc(0) = markerPointLocArray[frameToShow](markerDragHitIndex, 0);
        selectedLoc(1) = markerPointLocArray[frameToShow](markerDragHitIndex, 1);
        selectedLoc(2) = markerPointLocArray[frameToShow](markerDragHitIndex, 2) + deltaZ;

        Eigen::MatrixXd pointColor(1, 3);
        pointColor << 0.0, 1.0, 1.0;
        viewer.data().add_points(selectedLoc, pointColor);
    }

    // add a new visualization point at tip of mouse
    if (markerDragFocused) {

        Eigen::MatrixXd destLoc(1, 3);
        static const double deltaZ = 0.3;
        destLoc(0) = markerDragLoc(1) + 0.5;
        destLoc(1) = (imgRows-0.5) - markerDragLoc(0);
        destLoc(2) = markerPointLocArray[frameToShow](markerDragHitIndex, 2);

        viewer.data().add_points(destLoc, markerPointColor);
    }
}

}  // namespace zebrafish
