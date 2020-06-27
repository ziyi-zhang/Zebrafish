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


// hypre solver
extern int solverMaxIt;       // Hypre solver max iterations
extern double solverConvTol;  // Hypre solver convergence tolerance
extern double solverTol;      // Hypre solver tolerance

// pre-process
extern double pixelQuantile;  // pixels with value larger than this quantile will be trimmed

// cylinder / enerygy function
extern int diskQuadMethod;    // index of the disk quadrature method
extern double minRadius;      // minimal radius of a cylinder
extern int heightLayers;      // desired height layers used in quadrature

}  // namespace zebrafish
