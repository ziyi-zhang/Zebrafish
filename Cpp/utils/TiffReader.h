#pragma once


#include <Eigen/Dense>


namespace zebrafish
{
    bool read_tif_image(const std::string &path, Eigen::MatrixXd &img);
}
