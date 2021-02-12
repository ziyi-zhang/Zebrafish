#include <zebrafish/zebra-analysis.hpp>

#include <polyfem/State.hpp>
#include <polyfem/GenericProblem.hpp>
#include <polyfem/VTUWriter.hpp>
#include <polyfem/PointBasedProblem.hpp>

#include <igl/barycenter.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/writeMESH.h>
#include <igl/upsample.h>
#include <igl/doublearea.h>
#include <igl/write_triangle_mesh.h>
#include <igl/remove_unreferenced.h>
#include <highfive/H5Easy.hpp>

namespace zebrafish
{
    void compute_analysis(const std::vector<Eigen::MatrixXd> &V, const Eigen::MatrixXi &FF,
                          const std::string &path,
                          const double E, const double nu,
                          const double offset, const double min_area,
                          const int discr_order, const bool is_linear, const int n_refs, const double vismesh_rel_area,
                          const bool saveinput)
    {

        //const std::string rbf_function = "gaussian";
        const std::string rbf_function = "thin-plate";
        // const std::string rbf_function = "cubic";
        const double eps = 10;

        Eigen::MatrixXi F;
        Eigen::MatrixXd V0;
        const Eigen::MatrixXd &ref_v = V.front();
        igl::upsample(ref_v, FF, V0, F, 2);
        Eigen::MatrixXd barys;
        igl::barycenter(V0, F, barys);

        const auto set_bc = [&barys](const Eigen::MatrixXd &v, const bool boundary) {
            if (boundary)
                return 1;
            double min = std::numeric_limits<double>::max();
            int min_index = 0;
            for (int i = 0; i < barys.rows(); ++i)
            {
                const double norm = (v - barys.row(i)).squaredNorm();
                if (norm < min)
                {
                    min_index = i;
                    min = norm;
                }
            }

            if (min < 1)
                return 2;
            return 0;
        };

        const auto WriteInputToFile = [&V, &FF, &path, &E, &nu, &offset, &min_area, &discr_order, &is_linear, &n_refs, &vismesh_rel_area]() {
            const int frames = V.size();
            const int Nverts = V[0].rows();

            Eigen::MatrixXd V_concatenated(frames * Nverts, 3);
            for (int i = 0; i < frames; i++)
            {
                V_concatenated.block(Nverts * i, 0, Nverts, 3) = V[i];
            }

            std::string fileName = path + "-AnalysisInput" + ".h5";
            H5Easy::File file(fileName, H5Easy::File::ReadWrite | H5Easy::File::Create);
            H5Easy::dump(file, "E", E);
            H5Easy::dump(file, "nu", nu);
            H5Easy::dump(file, "offset", offset);
            H5Easy::dump(file, "min_area", min_area);
            H5Easy::dump(file, "discr_order", discr_order);
            H5Easy::dump(file, "is_linear", is_linear);
            H5Easy::dump(file, "n_refs", n_refs);
            H5Easy::dump(file, "vismesh_rel_area", vismesh_rel_area);

            H5Easy::dump(file, "frames", frames);
            H5Easy::dump(file, "Nverts", Nverts);
            H5Easy::dump(file, "V", V_concatenated);
            H5Easy::dump(file, "F", FF);
        };

        if (saveinput)
            WriteInputToFile();

        //////////////////////////////////////////////////////////////
        // get the bounding box for the first V
        Eigen::MatrixXd box_min = ref_v.colwise().minCoeff();
        Eigen::MatrixXd box_max = ref_v.colwise().maxCoeff();
        assert(box_min.size() == 3);
        assert(box_max.size() == 3);

        // get the bounding box for all V's
        for (int i = 1; i < V.size(); ++i)
        {
            const Eigen::MatrixXd cbox_min = V[i].colwise().minCoeff();
            const Eigen::MatrixXd cbox_max = V[i].colwise().maxCoeff();

            for (int d = 0; d < 3; ++d)
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

        igl::upsample(cube_v, cube_f, cube_v, cube_f, 3);

        Eigen::MatrixXd mesh_v(cube_v.rows() + V0.rows(), 3);
        Eigen::MatrixXi mesh_f(cube_f.rows() + F.rows(), 3);

        mesh_v.block(0, 0, cube_v.rows(), 3) = cube_v;
        mesh_v.block(cube_v.rows(), 0, V0.rows(), 3) = V0;

        mesh_f.block(0, 0, cube_f.rows(), 3) = cube_f;
        mesh_f.block(cube_f.rows(), 0, F.rows(), 3) = F;
        mesh_f.block(cube_f.rows(), 0, F.rows(), 3).array() += cube_v.rows();

        //////////////////////////////////////////////////////////////
        // Generate tet mesh
        Eigen::MatrixXd nodes;
        Eigen::MatrixXi elem, faces;
        const std::string switches = "zpq1.414a" + std::to_string(min_area) + "QYY";
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
        //     igl::write_triangle_mesh("text1.obj", V0, F);
        // }

        Eigen::RowVector3d zero(0, 0, 0);
        for (int sim = 1; sim < V.size(); ++sim)
        {
            const json args = {
                {"tensor_formulation", is_linear ? "LinearElasticity" : "NeoHookean"},
                {"discr_order", discr_order},

                {"params", {{"E", E}, {"nu", nu}}},

                {"problem", "PointBasedTensor"}};

            const Eigen::MatrixXd currentv = V[sim] - ref_v;
            // Eigen::MatrixXd currentv = V[sim] - ref_v;
            // currentv.setConstant(.5);

            state.init(args);
            polyfem::PointBasedTensorProblem &tproblem = *dynamic_cast<polyfem::PointBasedTensorProblem *>(state.problem.get());
            tproblem.add_constant(1, zero);
            Eigen::Matrix<bool, 3, 1> dirichet_dims;
            dirichet_dims.setConstant(true);
            const Eigen::MatrixXd twod = ref_v.block(0, 0, ref_v.rows(), 2);
            tproblem.add_function(2, currentv, twod, rbf_function, eps, 2, dirichet_dims);
            state.solve();

            Eigen::MatrixXd vals;
            Eigen::MatrixXd traction_forces, traction_forces_flip;
            Eigen::MatrixXd stress, stress_flip;
            Eigen::MatrixXd mises, mises_flip;
            Eigen::MatrixXi F_flip = F;
            F_flip.col(0) = F.col(1);
            F_flip.col(1) = F.col(0);

            state.interpolate_boundary_function_at_vertices(V0, F, state.sol, vals);
            state.interpolate_boundary_tensor_function(V0, F, state.sol, vals, true, traction_forces, stress, mises);
            state.interpolate_boundary_tensor_function(V0, F_flip, state.sol, vals, true, traction_forces_flip, stress_flip, mises_flip);

            Eigen::MatrixXd vertex_traction_forces(V0.rows(), 3);
            vertex_traction_forces.setZero();

            Eigen::MatrixXd vertex_stress(V0.rows(), 9);
            vertex_stress.setZero();

            Eigen::MatrixXd vertex_mises(V0.rows(), 1);
            vertex_mises.setZero();

            Eigen::MatrixXd vertex_traction_forces_flip(V0.rows(), 3);
            vertex_traction_forces_flip.setZero();

            Eigen::MatrixXd vertex_stress_flip(V0.rows(), 9);
            vertex_stress_flip.setZero();

            Eigen::MatrixXd vertex_mises_flip(V0.rows(), 1);
            vertex_mises_flip.setZero();

            Eigen::MatrixXd area, vertex_area(V0.rows(), 1);
            vertex_area.setZero();
            igl::doublearea(V0, F, area);

            for (int f = 0; f < F.rows(); ++f)
            {
                for (int d = 0; d < 3; ++d)
                {
                    vertex_traction_forces.row(F(f, d)) += traction_forces.row(f) * area(f);
                    vertex_traction_forces_flip.row(F(f, d)) += traction_forces_flip.row(f) * area(f);

                    vertex_stress.row(F(f, d)) += stress.row(f) * area(f);
                    vertex_stress_flip.row(F(f, d)) += stress_flip.row(f) * area(f);

                    vertex_mises(F(f, d)) += mises(f) * area(f);
                    vertex_mises_flip(F(f, d)) += mises_flip(f) * area(f);

                    vertex_area(F(f, d)) += area(f);
                }
            }

            for (int d = 0; d < 3; ++d)
            {
                vertex_traction_forces.col(d).array() /= vertex_area.array();
                vertex_traction_forces_flip.col(d).array() /= vertex_area.array();
            }
            for (int d = 0; d < 9; ++d)
            {
                vertex_stress.col(d).array() /= vertex_area.array();
                vertex_stress_flip.col(d).array() /= vertex_area.array();
            }
            vertex_mises.array() /= vertex_area.array();
            vertex_mises_flip.array() /= vertex_area.array();

            const std::string out_path = path + "-frame" + std::to_string(sim) + ".vtu";
            std::cerr << out_path << std::endl;
            state.save_vtu(out_path, 0);
            polyfem::VTUWriter writer;
            writer.add_field("displacement", vals);
            writer.add_field("traction_forces", vertex_traction_forces);
            writer.add_field("traction_forces_other", vertex_traction_forces_flip);

            for(int i = 0; i < 3; ++i){
                for(int j = 0; j < 3; ++j){
                    const std::string index  = "_" + std::to_string(i) + std::to_string(j);
                    writer.add_field("stress" + index, vertex_stress.col(i*3+j));
                    writer.add_field("stress_other" + index, vertex_stress_flip.col(i*3+j));
                }
            }

            writer.add_field("von_mises", vertex_mises);
            writer.add_field("von_mises_other", vertex_mises_flip);


            writer.write_mesh(out_path, V0, F);

            // state.save_vtu(out_path + ".all.vtu", 0);
        }
    }
} // namespace zebrafish
