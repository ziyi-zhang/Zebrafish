#pragma once

#include <Eigen/Core>

#include <vector>

namespace zebrafish {
    void compute_analysis(std::vector<Eigen::MatrixXd> &V, Eigen::MatrixXi &F,
    const std::string &path,
    const double E, const double nu,
     const double offset, const double min_area,
     const int discr_order, const bool is_linear, const int n_refs, const double vismesh_rel_area,
     const bool saveinput);
}
