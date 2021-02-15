// cylinder test
#include <zebrafish/Common.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/Padding.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include <imgui/imgui.h>
#include <catch.hpp>

#include <math.h>
#include <stdlib.h>
#include <random>
#include <unordered_set>

using namespace std;
using namespace Eigen;
using namespace zebrafish;


const auto Initialize = [](Eigen::MatrixXd &V, Eigen::MatrixXi &F, RCMap_t &RCMap) {
    F.resize(8, 3);
    F << 8, 6, 5, 
         6, 7, 3, 
         4, 5, 2, 
         5, 3, 1, 
         4, 8, 5, 
         5, 6, 3, 
         2, 5, 1, 
         1, 3, 0;
    double s3 = std::sqrt(3);
    double s32 = s3 * 2.0;
    V.resize(9, 3);
    V << s32, 5, 1, 
         s32, 3, 1, 
         s32, 1, 1, 
         s3, 4, 1, 
         s3, 0, 1, 
         s3, 2, 1, 
         0, 3, 1, 
         0, 5, 1, 
         0, 1, 1;
    RCMap.insert({0, {2, 3}});
    RCMap.insert({1, {2, 2}});
    RCMap.insert({2, {2, 1}});
    RCMap.insert({3, {1, 2}});
    RCMap.insert({4, {1, 0}});
    RCMap.insert({5, {1, 1}});
    RCMap.insert({6, {0, 2}});
    RCMap.insert({7, {0, 3}});
    RCMap.insert({8, {0, 1}});
};


////////////////////////////////////////////////////////////////////////////

TEST_CASE("Padding_9", "[PaddingTest]") {

    Eigen::MatrixXd V, appendV;
    Eigen::MatrixXi F, appendF;
    RCMap_t RCMap;

    Initialize(V, F, RCMap);

    padding::ComputeOneRing(V, F, RCMap, appendV, appendF);
    padding::AddOneRing<Eigen::MatrixXd>(appendV, appendF, V, F);


    cout << "appendV" << endl;
    cout << appendV << endl;
    cout << "appendF" << endl;
    cout << appendF << endl;

    std::unordered_set<int> vids;
    for (int i=0; i<F.rows(); i++) {
        vids.insert(F(i, 0));
        vids.insert(F(i, 1));
        vids.insert(F(i, 2));
    }
    igl::opengl::glfw::Viewer viewer;
    Eigen::Vector3d p;
    p << 0.0, 0.0, 0.005;
    viewer.data().set_mesh(V, F);
    for (int vid : vids) {
        viewer.data().add_label(V.row(vid) + p.transpose(), std::to_string(vid));
    }
    printf("#vids = %lu\n", vids.size());

    // activate label rendering
    viewer.data().show_labels = true;
    // Rendering of text labels is handled by ImGui, so we need to enable the ImGui
    // plugin to show text labels.
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = [](){};
    viewer.plugins.push_back(&menu);
    viewer.launch();
}
