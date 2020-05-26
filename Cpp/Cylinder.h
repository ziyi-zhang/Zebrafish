#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

#define heightLayers 5  // desired height layers used in quadrature

////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

typedef struct cylinder_t {
// uniquely defines a cylinder
    double x, y, z, r, h;  // bottom x(row), y(col), z + radius + height
} cylinder_t;


typedef struct samplePoints {
//
    Eigen::MatrixX3d points;  // [#points x 3] point positions
    Eigen::VectorXd weights;  // [#points] (Gaussian quadrature weight)*(energy function constants)*(subtraction sigma function)
} samplePoints_t;

///////////////////////////////////////

class cylinder {

private:
    cylinder_t cyl;
    double energy;  // evaluated energy function value
    samplePoints_t samplePoints;

public:
    bool SampleCylinder(const Eigen::MatrixXd &image);
    /// Calculate the sample points (both interior and exterior point) and 
    /// save them to corresponding private variables
    /// 'image' is only used to do boundary check
    ///
    /// @param[in]  image     { 3D matrix of the image }
    /// @return      {whether the cylinder is inside the boundary of the image}
    ///

    void EvaluateCylinder(const Eigen::MatrixXd &image);
    /// Calculate the energy for this cylinder and save the value to a private variable "energy"

    // maintenance methods
    cylinder();
    ~cylinder();
};

}  // namespace zebrafish
