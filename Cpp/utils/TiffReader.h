#pragma once


#include <Eigen/Dense>

#include <vector>


namespace zebrafish
{
    bool read_tif_image(const std::string &path, std::vector<Eigen::MatrixXd> &img);
}
