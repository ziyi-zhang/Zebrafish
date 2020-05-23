#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

typedef struct cylinder {
// uniquely defines a cylinder
    double x, y, z, r, h;  // bottom x, y, z + radius + height
} cylinder_t;


class cylinder {

private:
    cylinder_t cyl;

public:
    void SampleCylinder();
    /// Calculate the sample points (both interior and exterior) and save them to private variables
    
    void EvaluateCylinder();
    /// Calculate 

    // maintenance methods
    cylinder();
    ~cylinder();
};

}  // namespace zebrafish
