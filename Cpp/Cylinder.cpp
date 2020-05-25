#include <Cylinder.h>
#include <cmath>

namespace zebrafish {

bool cylinder::SampleCylinder(const Eigen::MatrixXd &image) {

    const double &x = this->cyl.x;
    const double &y = this->cyl.y;
    const double &z = this->cyl.z;
    const double &r = this->cyl.r;
    const double &h = this->cyl.h;
    const int &xmax = image.rows();
    const int &ymax = image.cols();

    // boundary check
    if (x < -0.5 + r*1.5 || x > xmax-0.5-r*1.5 ||
        y < -0.5 + r*1.5 || y > ymax-0.5-r*1.5 ||
        z < 0 || z > ) {
        return false;
    }
}

}  // namespace zebrafish
