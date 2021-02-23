// point_mesh_dist test
#include <zebrafish/Common.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/Cage.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>

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


const auto Init1 = [](Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
    F.resize(22, 3);
    F << 0, 1, 3, 
        3, 1, 2, 
        4, 5, 7, 
        7, 5, 6, 
        9, 8, 11, 
        9, 11, 10, 
        0, 1, 4, 
        4, 1, 5, 
        2, 3, 7, 
        2, 7, 6, 
        8, 9, 0, 
        0, 9, 1, 
        10, 11, 3, 
        10, 3, 2, 
        1, 2, 5, 
        5, 2, 6, 
        3, 0, 4, 
        3, 4, 7, 
        9, 10, 1, 
        1, 10, 2, 
        11, 8, 0, 
        11, 0, 3;

    V.resize(12, 3);
    V << 0, 0, 0, 
        1, 0, 0, 
        1, 1, 0, 
        0, 1, 0, 
        0, 0, 1, 
        1, 0, 1, 
        1, 1, 1, 
        0, 1, 1, 
        0, 0, -1, 
        1, 0, -1, 
        1, 1, -1, 
        0, 1, -1;
};


TEST_CASE("pmd", "[PointMeshDistTest]") {

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    Init(V, F);

    Eigen::MatrixXd p(4, 3);
    p << 3.2, 8, 1, 
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


TEST_CASE("pmd2", "[PointMeshDistTest]") {

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    Init1(V, F);

    const std::string switches = "zpq1.4VYY";
    Eigen::MatrixXd TV;
    Eigen::MatrixXi TT, TF;
    igl::copyleft::tetgen::tetrahedralize(V, F, switches, TV, TT, TF);

    cout << "TV" << endl << TV << endl;
    cout << "TT" << endl << TT << endl;
    cout << "TF" << endl << TF << endl;

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
