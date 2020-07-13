#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/autodiff.h>

#include <cassert>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

///////////
// types //
///////////
typedef std::vector<Eigen::MatrixXd> image_t;  // a 3D double matrix
typedef std::vector<image_t> imageData_t;  // a sequence of 3D image
typedef Eigen::Vector3d gradient_t;  // gradient
typedef DScalar1<double, gradient_t> DScalar;

}  // namespace zebrafish
