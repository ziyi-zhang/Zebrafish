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
typedef Eigen::Vector3d gradient_t;  // gradient
typedef DScalar1<double, gradient_t> DScalar;

//////////////////////
// global variables //
//////////////////////


// pre-process
extern double pixelQuantile;  // pixels with value larger than this quantile will be trimmed

}  // namespace zebrafish
