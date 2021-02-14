#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/Common.h>

#include <Eigen/Dense>
#include <map>
#include <array>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

typedef std::map<int, std::array<int, 2> > RCMap_t;

////////////////////////////////////////////////
// padding [static class]

class padding {

private:

public:
    static void ComputeOneRing(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const RCMap_t &RCMap, Eigen::MatrixXd &appendV, Eigen::MatrixXi &appendF);
    static void AddOneRing(const Eigen::MatrixXd &appendV, const Eigen::MatrixXi &appendF, Eigen::MatrixXd &V, Eigen::MatrixXi &F);
};

}  // namespace zebrafish
