// point_mesh_dist test
#include <zebrafish/Common.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/Cage.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/point_mesh_squared_distance.h>

#include <imgui/imgui.h>
#include <catch.hpp>

#include <math.h>
#include <stdlib.h>
#include <random>
#include <unordered_set>

using namespace std;
using namespace Eigen;
using namespace zebrafish;


////////////////////////////////////////////////////////////////////////////

const auto Init = [](Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
    F.resize(4, 3);
    F << 5, 4, 3, 
         2, 3, 0, 
         2, 5, 3, 
         3, 4, 1;
    V.resize(6, 3);
    V << 5, 13, 1, 
         14, 8, 1, 
         3, 8, 1, 
         9, 8, 1, 
         11, 3, 1,
         5, 3, 1;
};


TEST_CASE("pmd", "[PointMeshDistTest]") {

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    Init(V, F);

    Eigen::MatrixXd p(4, 3);
    p << 5, 13, 1, 
         5, 13, 2.2, 
         0, 0, 0, 
         0, 8, 1;

    Eigen::VectorXd squaredDist;
    Eigen::MatrixXd I, C;
    igl::point_mesh_squared_distance(p, V, F, squaredDist, I, C);

    cout << squaredDist << endl;

    /*
    ///////////////////////////////////
    // Viewer
    igl::opengl::glfw::Viewer viewer;
    Eigen::Vector3d p;
    p << 0.0, 0.0, 0.005;
    viewer.data().set_mesh(V, F);
    for (int i=0; i<V.rows(); i++) {
        viewer.data().add_label(V.row(i) + p.transpose(), std::to_string(i));
    }

    // activate label rendering
    viewer.data().show_labels = true;
    // Rendering of text labels is handled by ImGui, so we need to enable the ImGui
    // plugin to show text labels.
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = [](){};
    viewer.plugins.push_back(&menu);
    viewer.launch();
    */
}
