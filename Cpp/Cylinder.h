#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {



class cylinder {

private:
    double x, y, z, r, h; // bottom x, y, z + radius + height

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
