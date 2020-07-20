#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>


namespace zebrafish {

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 4: Optimization

void GUI::DrawStage4() {

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::Text("Optimization");
    if (ImGui::Button("Start Optimization")) {
        
        Optimization();
    }

    ImGui::Text("Stage 4");
}


////////////////////////////////////////////////////////////////////////////////////////
// Optimization


void GUI::Optimization() {

    logger().info("Not implemented");
}


void GUI::UpdateSampleNewton(const Eigen::MatrixXd &gridSampleInput, const Eigen::MatrixXd &gridSampleOutput) {

    assert(gridSampleInput.rows() == gridSampleOutput.rows());
    const int N = gridSampleOutput.rows();
    int M = 0;
    int i, count = 0;

    for (i=0; i<N; i++)
        if (gridSampleOutput(i) < gridEnergyThres) M++;  // count to reserve space

    pointRecord.num = M;
    pointRecord.alive.resize(M, Eigen::NoChange);
    pointRecord.grid_search.resize(M, Eigen::NoChange);
    pointRecord.optimization.resize(M, Eigen::NoChange);

    for (i=0; i<N; i++) {
        if (gridSampleOutput(i) > gridEnergyThres) continue;

        // alive
        pointRecord.alive(count) = true;

        // grid search
        pointRecord.grid_search(count, 0) = gridSampleInput(i, 0);
        pointRecord.grid_search(count, 1) = gridSampleInput(i, 1);
        pointRecord.grid_search(count, 2) = gridSampleInput(i, 2);
        pointRecord.grid_search(count, 3) = gridSampleInput(i, 3);
        pointRecord.grid_search(count, 4) = gridSampleOutput(i);

        // optimization
        pointRecord.optimization(count, 0) = 0;
        pointRecord.optimization(count, 1) = 0;
        pointRecord.optimization(count, 2) = 0;
        pointRecord.optimization(count, 3) = 0;
        pointRecord.optimization(count, 4) = 0;
        pointRecord.optimization(count, 5) = 0;

        count++;
    }
}

}  // namespace zebrafish
