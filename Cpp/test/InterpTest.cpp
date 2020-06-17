// Interp TEST
#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <math.h>
#include <stdlib.h>
#include <random>
#include <queue>
#include <vector>

using namespace std;
using namespace Eigen;
using namespace zebrafish;

double func(double x, double y, double z) {

    // return x + y;

    // x^2 + y^2
    return (x-14.5)*(x-14.5) + (y-14.5)*(y-14.5);

    // (x^2 + y^2)^(3/2)
    // return pow((x-14.5)*(x-14.5) + (y-14.5)*(y-14.5), 1.5);
    // return (x-14.5)*(x-14.5)*(x-14.5) + (y-14.5)*(y-14.5)*(y-14.5);

    // x^4 + y^4 + 2 * x^2 * y^2
    // return (x-14.5)*(x-14.5)*(x-14.5)*(x-14.5) + (y-14.5)*(y-14.5)*(y-14.5)*(y-14.5) +
    //     2 * (y-14.5)*(y-14.5)* (x-14.5)*(x-14.5);

    // x^5 + y^5
    // return (x-14.5)*(x-14.5)*(x-14.5)*(x-14.5)*(x-14.5) + (y-14.5)*(y-14.5)*(y-14.5)*(y-14.5)*(y-14.5);

    // (x^2 + y^2)^(5/2)
    // return pow((x-14.5)*(x-14.5) + (y-14.5)*(y-14.5), 2.5);
}

DECLARE_DIFFSCALAR_BASE();

int main() {

    image_t image;  // 30 * 30 * 10
    int sizeX, sizeY, sizeZ;
    double x, y, z;

    sizeX = 25;  // 0, 1, ..., 29
    sizeY = 25;
    sizeZ = 25;

    // user input
    resolutionX = 0.325;
    resolutionY = 0.325;
    resolutionZ = 0.5;

    // generate sample grid (3D)
    double maxPixel = 0;
    for (z=0; z<sizeZ; z++) {
        
        MatrixXd layer(sizeX, sizeY);
        for (x=0; x<sizeX; x++)
            for (y=0; y<sizeY; y++) {
                layer(x, y) = func(x, y, z);
            }

        image.push_back(layer);
        if (layer.maxCoeff() > maxPixel) maxPixel = layer.maxCoeff();
    }

    // normalize it
    /*
    for (z=0; z<sizeZ; z++) {
        
        MatrixXd &layer = image[z];
        layer.array() /= maxPixel;
    }
    */

    // prepare B-spline
    bspline bsplineSolver;
    bsplineSolver.CalcControlPts(image, 0.7, 0.7, 1);

    // random
    srand(time(NULL));
    uniform_real_distribution<double> unif(0, sizeX-1);
    default_random_engine re(0);

    // moving median
    double tt, l_;
    priority_queue<double, vector<double>, greater<double> > u;
    priority_queue<double, vector<double>, less<double> > l;

    // interp test
    double err, sumerr = 0, minerr = 1.0, maxerr = 0.0;
    int trialNum = 500;
    for (int i = 0; i<trialNum; i++) {
        
        x = unif(re);
        y = unif(re);
        z = unif(re);
            // x = rand() % (sizeX-1)+0.2;
            // y = rand() % (sizeY-1)+0.3;
            // z = rand() % (sizeZ-1)+0.1;
        // x = 0; y = 0; z = 5;
        // cout << x << " " << y << " " << z << ":      ";

        Eigen::Matrix<DScalar, Eigen::Dynamic, 3> in;
        in.resize(1, 3);
        in << DScalar(x), DScalar(y), DScalar(z);
        Eigen::Matrix<DScalar, Eigen::Dynamic, 1> out;

        bsplineSolver.Interp3D(in, out);
        err = func(x, y, z) - out(0).getValue();
        // cout << err << endl;
        err = fabs(err);

        sumerr += err;
        if (err > maxerr) maxerr = err;
        if (err < minerr) minerr = err;

        // moving median
        if (l.empty()) {l.push(err); continue;}
        l_ = l.top();
        if (err >= l_) {u.push(err);} else {l.push(err);}
        if (u.size() > l.size()) {
            tt = u.top();
            u.pop();
            l.push(tt);
        }
        if (l.size() > u.size() + 1) {
            tt = l.top();
            l.pop();
            u.push(tt);
        }
    }

    cout << "Degree = " << degree << endl;
    cout << "Mean error = " << sumerr / double(trialNum) << endl;
    cout << "Median error = " << l.top() << endl;
    cout << "Min  error = " << minerr << endl;
    cout << "Max  error = " << maxerr << endl;
}
