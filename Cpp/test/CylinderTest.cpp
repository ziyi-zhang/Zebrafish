#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>

using namespace std;
using namespace Eigen;
using namespace zebrafish;

double func(double x, double y, double z) {

    return (x-14.5)*(x-14.5) + (y-14.5)*(y-14.5);
}

int main() {

    image_t image;  // 30 * 30 * 10
    int sizeX, sizeY, sizeZ;
    int x, y, z;

    sizeX = 30;  // 0, 1, ..., 29
    sizeY = 30;
    sizeZ = 10;

    // generate sample grid (3D)
    for (z=0; z<sizeZ; z++) {
        
        MatrixXd layer(sizeX, sizeY);
        for (x=0; x<sizeX; x++)
            for (y=0; y<sizeY; y++) {
                layer(x, y) = func(x, y, z);
            }

        image.push_back(layer);
    }

    // prepare B-spline
    bspline bsplineSolver;
    bsplineSolver.CalcControlPts(image, 1, 1, 1);

    // prepare cylinder
    cylinder cylinder;
    if (!cylinder.SampleCylinder(image, bsplineSolver, 14.5, 14.5, 3, 5, 4)) {
        cerr << "Invalid cylinder";
    }
    cout << "Evaluated result: " << cylinder.EvaluateCylinder(image, bsplineSolver) << endl;
}
