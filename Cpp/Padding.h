#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/Common.h>

#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

typedef std::map<int, std::array<int, 2> > RCMap_t;

////////////////////////////////////////////////
// padding [static class]

class padding {

private:

public:
    static void AddOneRing(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const RCMap_t &RCMap, Eigen::MatrixXd &appendV, Eigen::MatrixXi &appendF);
};

}  // namespace zebrafish
