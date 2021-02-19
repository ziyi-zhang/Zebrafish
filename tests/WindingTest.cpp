// in-out filter test
#include <zebrafish/Common.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/Padding.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/winding_number.h>

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

TEST_CASE("winding", "[WindingTest]") {

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    const double t = std::sqrt(3) / 2.0;
    V.resize(6, 3);
    V << 0, 0, 0, 
         1, 0, 0, 
         2, 0, 0, 
         0.5, t, 0, 
         1.5, t, 0, 
         1, t*2, 0;
    F.resize(4, 3);  // normal pointing positive z
    F << 0, 1, 3, 
         1, 2, 4, 
         3, 4, 5, 
         3, 2, 4;

    Eigen::MatrixXd sample(8, 3);
    sample << 1, 1, -1, 
              1, 1, 1, 
              1, 1, 0, 
              0, 2, 0,
              0, 0, -0.1, 
         1, 0, -0.1, 
         2, 0, -0.1, 
         0.5, t, -0.1;

    Eigen::VectorXd res;
    igl::winding_number(V, F, sample, res);

    std::cout << res << std::endl;
}
