#include <zebrafish/GUI.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/FileDialog.h>
#include <zebrafish/Logger.hpp>


namespace zebrafish {

namespace {

void GridSearch() {

    // prepare sampleInput & sampleOutput for grid search
    MatrixXd sampleInput, sampleOutput;
    const int Nx = bsplineSolver.Get_Nx();
    const int Ny = bsplineSolver.Get_Ny();
    const int Nz = bsplineSolver.Get_Nz();
    vector<double> rArray = {3, 4, 5, 6};
    const int Nr = rArray.size();
    const double gapX = 1.0;
    const double gapY = 1.0;
    const double gapZ = 1.0;
    double xx, yy, zz, rr;
    int sampleCount = 0;
    sampleInput.resize(Nx*Ny*Nz*Nr, 4);
    for (int ix=0; ix<floor(Nx/gapX); ix++)
        for (int iy=0; iy<floor(Ny/gapY); iy++)
            for (int iz=0; iz<floor(Nz/gapZ); iz++)
                for (int ir=0; ir<Nr; ir++) {

                    xx = ix * gapX;
                    yy = iy * gapY;
                    zz = iz * gapZ;
                    rr = rArray[ir];
                    if (!ValidStartingPoint(image, bsplineSolver, xx, yy, zz, rr)) continue;

                    sampleInput(sampleCount, 0) = xx;
                    sampleInput(sampleCount, 1) = yy;
                    sampleInput(sampleCount, 2) = zz;
                    sampleInput(sampleCount, 3) = rr;
                    sampleCount++;
                }
    sampleOutput.resize(sampleCount, 1);
    logger().info("Grid search #starting points = {}", sampleCount);
    logger().info("Grid search samples prepared...");

    // Search
    logger().info("Before search");
    tbb::parallel_for( tbb::blocked_range<int>(0, sampleCount),
        [&sampleInput, &sampleOutput, &bsplineSolver](const tbb::blocked_range<int> &r) {

            for (int ii = r.begin(); ii != r.end(); ++ii) {
                cylinder::EvaluateCylinder(bsplineSolver, sampleInput(ii, 0), sampleInput(ii, 1), sampleInput(ii, 2), sampleInput(ii, 3), 3, sampleOutput(ii));
            }
        });
    logger().info("After search");
}

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 1: 

void GUI::DrawStage1() {

    if (ImGui::Button("Compute Bspline")) {

        logger().info("Computing Bspine for the first frame");
        const int bsplineDegree = 2;
        bsplineSolver.SetResolution(resolutionX, resolutionY, resolutionZ);
        bsplineSolver.CalcControlPts(img, 0.7, 0.7, 0.7, bsplineDegree);
    }

    ImGui::Separator();

    ImGui::Text("Grid Search");

    ImGui::Separator();

    ImGui::Text("Optimization");

    ImGui::Text("Stage 1");
}

}  // namespace zebrafish
