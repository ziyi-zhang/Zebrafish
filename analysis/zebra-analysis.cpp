#include <zebrafish/zebra-analysis.hpp>

#include <polyfem/State.hpp>
#include <polyfem/GenericProblem.hpp>

#include <igl/barycenter.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/writeMESH.h>
#include <igl/write_triangle_mesh.h>

namespace zebrafish {
    void compute_analysis(std::vector<Eigen::MatrixXd> &V, Eigen::MatrixXi &F,
    const std::string &path,
    const double E, const double nu,
     const double offset, const double min_area,
     const double discr_order, const bool is_linear, const int n_refs, const double vismesh_rel_area)
    {

        const Eigen::MatrixXd &V0 = V.front();
        Eigen::MatrixXd barys;
        igl::barycenter(V0, F, barys);

        const auto set_bc = [&barys](const Eigen::MatrixXd &v, const bool boundary)
        {
            if(boundary)
                return 1;
            double min = std::numeric_limits<double>::max();
            int min_index = 0;
            for(int i = 0; i < barys.rows(); ++i)
            {
                const double norm = (v - barys.row(i)).squaredNorm();
                if(norm < min)
                {
                    min_index = i;
                    min = norm;
                }
            }

            if(min < 1)
                return 2;
            return 0;
        };

        const auto bc = [&barys, &V0, &F](const Eigen::MatrixXd &currentv, const double x, const double y, const double z)
        {
            const Eigen::Vector3d v(x, y, z);
            Eigen::MatrixXd vv(1, 3); vv << x, y, z;

            double min = std::numeric_limits<double>::max();
            int min_index = -1;
            for(int i = 0; i < barys.rows(); ++i)
            {
                const double norm = (vv - barys.row(i)).squaredNorm();
                if(norm < min)
                {
                    min_index = i;
                    min = norm;
                }
            }

            const Eigen::MatrixXi tri = F.row(min_index);
            const Eigen::Vector3d a = V0.row(tri(0));
            const Eigen::Vector3d b = V0.row(tri(1));
            const Eigen::Vector3d c = V0.row(tri(2));

            const Eigen::Vector3d e1 = a-b;
            const Eigen::Vector3d e2 = b-c;
            const Eigen::Vector3d e3 = c-a;

            const double area = e1.cross(e2).norm();

            const double lambda1 = e2.cross(v-c).norm()/area;
            const double lambda2 = e3.cross(v-a).norm()/area;
            const double lambda3 = e1.cross(v-b).norm()/area;

            const Eigen::MatrixXd v1 = currentv.row(tri(0));
            const Eigen::MatrixXd v2 = currentv.row(tri(1));
            const Eigen::MatrixXd v3 = currentv.row(tri(2));

            Eigen::MatrixXd res(1, 3);
            res = lambda1*v1+lambda2*v2+lambda3*v3;
            return res;
        };

        //////////////////////////////////////////////////////////////
        // get the bounding box for the first V
        Eigen::MatrixXd box_min = V0.colwise().minCoeff();
        Eigen::MatrixXd box_max = V0.colwise().maxCoeff();
        assert(box_min.size() == 3);
        assert(box_max.size() == 3);

        // get the bounding box for all V's
        for(int i = 1; i < V.size(); ++i)
        {
            const Eigen::MatrixXd cbox_min = V[i].colwise().minCoeff();
            const Eigen::MatrixXd cbox_max = V[i].colwise().maxCoeff();

            for(int d = 0; d < 3; ++d)
            {
                box_min(d) = std::min(box_min(d), cbox_min(d));
                box_max(d) = std::max(box_max(d), cbox_max(d));
            }
        }

        // offset -> larger box
        const double diag = (box_max - box_min).norm();
        box_min.array() -= offset * diag;
        box_max.array() += offset * diag;

        //////////////////////////////////////////////////////////////
        // generate bounding box mesh
        Eigen::MatrixXi cube_f(12, 3);
        cube_f << 4, 2, 0,
                    2, 7, 3,
                    6, 5, 7,
                    1, 7, 5,
                    0, 3, 1,
                    4, 1, 5,
                    4, 6, 2,
                    2, 6, 7,
                    6, 4, 5,
                    1, 3, 7,
                    0, 2, 3,
                    4, 0, 1;
        Eigen::MatrixXd cube_v(8, 3);
        cube_v << box_max(0), box_max(1), box_max(2),
                    box_max(0), box_max(1), box_min(2),
                    box_max(0), box_min(1), box_max(2),
                    box_max(0), box_min(1), box_min(2),
                    box_min(0), box_max(1), box_max(2),
                    box_min(0), box_max(1), box_min(2),
                    box_min(0), box_min(1), box_max(2),
                    box_min(0), box_min(1), box_min(2);

        Eigen::MatrixXd mesh_v(cube_v.rows()+V0.rows(), 3);
        Eigen::MatrixXi mesh_f(cube_f.rows()+F.rows(), 3);

        mesh_v.block(0, 0, cube_v.rows(), 3) = cube_v;
        mesh_v.block(cube_v.rows(), 0, V0.rows(), 3) = V0;

        mesh_f.block(0, 0, cube_f.rows(), 3) = cube_f;
        mesh_f.block(cube_f.rows(), 0, F.rows(), 3) = F;
        mesh_f.block(cube_f.rows(), 0, F.rows(), 3).array() += cube_v.rows();

        //////////////////////////////////////////////////////////////
        // Generate tet mesh
        Eigen::MatrixXd nodes;
        Eigen::MatrixXi elem, faces;
        const std::string switches = "zpq1.414a"+std::to_string(min_area)+"Q";
        // const std::string switches = "zpq1.414a10000Q";
        igl::copyleft::tetgen::tetrahedralize(mesh_v, mesh_f, switches, nodes, elem, faces);

        //////////////////////////////////////////////////////////////
        // polyfem
        polyfem::State state;
        state.mesh = std::make_unique<polyfem::Mesh3D>();
        state.mesh->build_from_matrices(nodes, elem);
		state.args["normalize_mesh"] = false;
		state.args["n_refs"] = n_refs;
		state.args["vismesh_rel_area"] = vismesh_rel_area;
        state.load_mesh();
        state.mesh->compute_boundary_ids(set_bc);

        // {
        //     Eigen::MatrixXd points;
        //     Eigen::MatrixXi faces;
        //     Eigen::MatrixXd sidesets;

        //     state.get_sidesets(points, faces, sidesets);
        //     igl::write_triangle_mesh("text.obj", points, faces);
        // }

        Eigen::RowVector3d zero(0, 0, 0);
        for(int sim = 1; sim < V.size(); ++sim)
        {
            const json args = {
                {"tensor_formulation", is_linear ? "LinearElasticity" : "NeoHookean"},
                {"discr_order", discr_order},

                {"params",{
                    {"E", E},
                    {"nu", nu}
                }},

                {"problem", "GenericTensor"}
            };

            const Eigen::MatrixXd currentv = V[sim] - V0;

            state.init(args);
            polyfem::GenericTensorProblem &tproblem = *dynamic_cast<polyfem::GenericTensorProblem *>(state.problem.get());
            tproblem.add_dirichlet_boundary(1, zero, true, true, true);
            tproblem.add_dirichlet_boundary(2, [&currentv, &bc](const double x, const double y, const double z){return bc(currentv, x, y, z); }, true, true, true);

            state.solve();

            const std::string out_path = path + std::to_string(sim) + ".vtu";
            std::cerr << out_path << std::endl;
			state.save_vtu(out_path, 0);
        }
    }
}
