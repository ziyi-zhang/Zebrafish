#include <Common.h>

namespace zebrafish {

// Bspline
int degree = 2;

// image info
double resolutionX;
double resolutionY;
double resolutionZ;

// hypre solver
int solverMaxIt = 20000;  // 1000
double solverConvTol = 1e-15;  // 1e-10
double solverTol = 1e-15;  // 1e-10

//
double pixelQuantile = 0.995;

double minRadius = 2.0;
int heightLayers = 5;

}  // namespace zebrafish
