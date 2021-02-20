#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/Common.h>

#include <Eigen/Dense>
#include <array>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

////////////////////////////////////////////////
// cage [static class]

class cage {

private:

public:
    static void ComputeCage(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXd &sideV, Eigen::MatrixXi &sideF, bool above = true);
    static void AddCage(const Eigen::MatrixXd &sideV, const Eigen::MatrixXi &sideF, Eigen::MatrixXd &V, Eigen::MatrixXi &F);
};

}  // namespace zebrafish
