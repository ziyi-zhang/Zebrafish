#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/common.h>

#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

void Interp3D(const image_t &image, const Eigen::MatrixX3d &sample, Eigen::VectorXd &res);
/// Calculate the interpolated B-spline result at "sample" points.
/// Note: this function does not check for input validity
///
/// @param[in]   image     { 3D matrix of the image }
/// @param[in]   sample    { #sample x 3 input sample points }
/// @param[out]  res       { #sample x 1 output interpolated results}
///

}  // namespace zebrafish
