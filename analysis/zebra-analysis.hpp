#pragma once

#include <Eigen/Core>

#include <vector>
#include <map>
#include <array>

namespace zebrafish
{
    void compute_analysis(const std::vector<Eigen::MatrixXd> &V, const Eigen::MatrixXi &F,
                          const std::string &path,
                          const double E, const double nu,
                          const double offset, const double radius_edge_ratio, const double min_area,
                          const int discr_order, const bool is_linear, const int n_refs, const double vismesh_rel_area, const int upsample,
                          const std::map<int, std::array<int, 2> > &markerRCMap,
                          const bool saveinput);
}
