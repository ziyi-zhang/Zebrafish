#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/autodiff.h>

#include <cassert>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

// types
typedef std::vector<Eigen::MatrixXd> image_t;  // a 3D double matrix
typedef Eigen::Vector3d gradient_t;  // gradient
typedef DScalar1<double, gradient_t> DScalar;

// global variables
extern int degree;            // B-spline degree
                              // "3" for cubic basis function
                              // "2" for quadratic basis function

extern double resolutionX;    // The distance between two pixels in X-direction (in micrometers)
extern double resolutionY;
extern double resolutionZ;

extern double pixelQuantile;  // pixels with value larger than this quantile will be trimmed

extern double minRadius;      // minimal radius of a cylinder
extern int heightLayers;      // desired height layers used in quadrature

}  // namespace zebrafish
