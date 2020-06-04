#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <cassert>
#include <Eigen/Dense>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

// types
typedef std::vector<Eigen::MatrixXd> image_t;  // a 3D double matrix

// global variables
extern double minRadius;  // minimal radius of a cylinder
extern int heightLayers;  // desired height layers used in quadrature

}  // namespace zebrafish
