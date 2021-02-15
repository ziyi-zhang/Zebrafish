#include <CLI/CLI.hpp>
#include <zebrafish/zebra-analysis.hpp>

#include <igl/read_triangle_mesh.h>

int main(int argc, char **argv)
{
    std::string path = "";
    std::string file_name = "";
    std::string out = "sim/out";
    int n_meshes = 2;

    double offset = 1;
    double radius_edge_ratio = 1.2;
    double min_area = 500;

    double E = 1e-3;
    double nu = 0.45;

    bool is_linear = true;
    int discr_order = 1;
    int n_refs = 0;
    double vismesh_rel_area = 0.00001;
    int upsample = 3;
    std::map<int, std::array<int, 2> > markerRCMap;  // dummy

    CLI::App command_line{"ZebraFishAnalysis"};
    command_line.add_option("-p,--path", path, "Input folder with objs to process")->check(CLI::ExistingDirectory);
    command_line.add_option("-f,--file_name", file_name, "Prefix file name for obj, each file is file_name<int>.obj")->required();
    command_line.add_option("-n,--n_meshes", n_meshes, "Input number of meshes in folder")->required();
    command_line.add_option("-o,--out", out, "Output file prefix, name will be out<int>.vtu")->required();

    command_line.add_option("--offset", offset, "Diagonal multiplier for box mesh");
    command_line.add_option("-a,--area", min_area, "Minimum tet area used by tetgen");

    command_line.add_option("--E", E, "Young's modulus")->required();
    command_line.add_option("--nu", nu, "Poisson's ratio")->required();

    command_line.add_flag("--non_linear", is_linear, "Use non-linear material");
    command_line.add_option("-d,--discr_order", discr_order, "Analysis discretization order");
    command_line.add_option("-r,--n_refs", n_refs, "Number of mesh uniform refinements");
    command_line.add_option("-v,--output_density", vismesh_rel_area, "Desnsity of the output visualization");

    try
    {
        command_line.parse(argc, argv);
    }
    catch (const CLI::ParseError &e)
    {
        return command_line.exit(e);
    }

    std::vector<Eigen::MatrixXd> V(n_meshes);
    Eigen::MatrixXi F;

    for (int i = 0; i < n_meshes; ++i)
    {
        igl::read_triangle_mesh(path + file_name + std::to_string(i) + ".obj", V[i], F);
    }

    zebrafish::compute_analysis(V, F, out, E, nu, offset, radius_edge_ratio, min_area, discr_order, is_linear, n_refs, vismesh_rel_area, upsample, markerRCMap, false);

    return EXIT_SUCCESS;
}
