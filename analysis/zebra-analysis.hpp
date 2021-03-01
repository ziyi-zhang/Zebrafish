#pragma once

#include <Eigen/Core>

#include <vector>
#include <map>
#include <array>

namespace zebrafish
{
    bool compute_analysis(const std::vector<Eigen::MatrixXd> &V, const Eigen::MatrixXi &F,
                          const int bm_offset_v, const int bm_offset_f,
                          const std::string &path,
                          const double E, const double nu,
                          const double min_area,
                          const int discr_order, const bool is_linear, const int n_refs, const double vismesh_rel_area, const int upsample,
                          const std::map<int, std::array<int, 2>> &markerRCMap, const int imgRows, const int imgCols, const int layerPerImg,
                          const double resolutionX, const double resolutionY, const double resolutionZ,
                          const bool saveinput, const bool aboveCage = false, const bool belowCage = false, 
                          const int frameStart=-1, const int frameEnd=-1);
}
