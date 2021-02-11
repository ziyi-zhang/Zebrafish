#pragma once

#include <Eigen/Core>

#include <vector>

namespace zebrafish
{
    void compute_analysis(const std::vector<Eigen::MatrixXd> &V, const Eigen::MatrixXi &F,
                          const double resolutionX, const double resolutionY, const double resolutionZ,
                          const std::string &path,
                          const double E, const double nu,
                          const double offset, const double min_area,
                          const double discr_order, const bool is_linear, const int n_refs, const double vismesh_rel_area);
}
