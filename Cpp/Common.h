#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/autodiff.h>

#include <cassert>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

enum MOUSE_TYPE {
    MOUSEDOWN = 1, 
    MOUSEUP   = 2,
    MOUSEMOVE = 3
};

enum REJECT_MODE {
    REJECT_SINGLE = 0,
    REJECT_AREA   = 1
};

///////////
// types //
///////////
typedef std::vector<Eigen::MatrixXd> image_t;  // a 3D double matrix
typedef std::vector<image_t> imageData_t;  // a sequence of 3D image
typedef Eigen::Vector3d gradient_t;  // gradient
typedef DScalar1<double, gradient_t> DScalar;

}  // namespace zebrafish
