#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/Common.h>

#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

double QuantileImage(const zebrafish::image_t &image, double q);
/// Calculate quantile(image(:), q) where 0 < q < 1
/// Note: for performance consideration, q should be larger than 0.95

}  // namespace zebrafish
