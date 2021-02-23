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
#include <igl/writeMESH.h>
#include <igl/remove_unreferenced.h>
#include <igl/winding_number.h>
#include <igl/point_mesh_squared_distance.h>
#include <highfive/H5Easy.hpp>

namespace zebrafish
{
    void compute_analysis(const std::vector<Eigen::MatrixXd> &VV, const Eigen::MatrixXi &F,
                          const int bm_offset_v, const int bm_offset_f,
                          const std::string &path,
                          const double E, const double nu,
                          const double offset, const double radius_edge_ratio, const double max_tet_vol,
                          const int discr_order, const bool is_linear, const int n_refs, const double vismesh_rel_area, const int upsample,
                          const std::map<int, std::array<int, 2>> &markerRCMap, const int imgRows, const int imgCols, const int layerPerImg,
                          const double resolutionX, const double resolutionY, const double resolutionZ,
                          const bool saveinput, bool useWindingNumber, bool windingNumberOtherSide)
    {
        // Necesary for the scaling!
        std::vector<Eigen::MatrixXd> V = VV;

        // RBF setup
        //const std::string rbf_function = "gaussian";
        const std::string rbf_function = "thin-plate";
        // const std::string rbf_function = "cubic";
        const double eps = 10;

        // save to analysis file for future re-analysis
        const auto WriteInputToFile = [&V, &F, &path, &E, &nu, &offset, &radius_edge_ratio, &max_tet_vol, &discr_order, &is_linear, &n_refs, &vismesh_rel_area, &upsample,
                                       &markerRCMap, &imgRows, &imgCols, &layerPerImg, &resolutionX, &resolutionY, &resolutionZ]() {
            const int frames = V.size();
            const int Nverts = V[0].rows();

            Eigen::MatrixXd V_concatenated(frames * Nverts, 3);
            for (int i = 0; i < frames; i++)
            {
                V_concatenated.block(Nverts * i, 0, Nverts, 3) = V[i];
            }

            Eigen::MatrixXi markerRCMap_mat(markerRCMap.size(), 3);
            int count = 0;
            for (auto it : markerRCMap)
            {
                markerRCMap_mat(count, 0) = it.first;
                markerRCMap_mat(count, 1) = it.second[0];
                markerRCMap_mat(count, 2) = it.second[1];
                count++;
            }

            std::string fileName = path + "-AnalysisInput" + ".h5";
            H5Easy::File file(fileName, H5Easy::File::ReadWrite | H5Easy::File::Create);
            H5Easy::dump(file, "E", E);
            H5Easy::dump(file, "nu", nu);
            H5Easy::dump(file, "offset", offset);
            H5Easy::dump(file, "radius_edge_ratio", radius_edge_ratio);
            H5Easy::dump(file, "max_tet_vol", max_tet_vol);
            H5Easy::dump(file, "discr_order", discr_order);
            H5Easy::dump(file, "is_linear", is_linear);
            H5Easy::dump(file, "n_refs", n_refs);
            H5Easy::dump(file, "vismesh_rel_area", vismesh_rel_area);
            H5Easy::dump(file, "upsample", upsample);

            // for re-analysis visualization
            H5Easy::dump(file, "imgRows", imgRows);
            H5Easy::dump(file, "imgCols", imgCols);
            H5Easy::dump(file, "layerPerImg", layerPerImg);
            H5Easy::dump(file, "resolutionX", resolutionX);
            H5Easy::dump(file, "resolutionY", resolutionY);
            H5Easy::dump(file, "resolutionZ", resolutionZ);

            // V, F
            // NOTE: V is in visualization format, with mean displacement removed, but before converting to physical unit
            H5Easy::dump(file, "frames", frames);
            H5Easy::dump(file, "Nverts", Nverts);
            H5Easy::dump(file, "V", V_concatenated);
            H5Easy::dump(file, "F", F);
            // for ring padding
            H5Easy::dump(file, "markerRCMap_mat", markerRCMap_mat);
        };

        if (saveinput)
            WriteInputToFile();

        // To physical unit
        if (resolutionX > 0 && resolutionY > 0 && resolutionZ > 0)
        {
            for (int i = 0; i < V.size(); i++)
            {
                V[i].col(0) *= resolutionY;
                V[i].col(1) *= resolutionX;
                V[i].col(2) *= resolutionZ;
            }
        }
        else
        {
            std::cout << "Analysis using pixel unit" << std::endl;
        }

        const Eigen::MatrixXd &V0 = V.front();
        const Eigen::MatrixXd bm_v = V0.block(0, 0, bm_offset_v, 3);
        const Eigen::MatrixXi bm_f = F.block(0, 0, bm_offset_f, 3);

        //////////////////////////////////////////////////////////////
        // Generate tet mesh
        Eigen::MatrixXd nodes;
        Eigen::MatrixXi elem, faces;
        // allow exterior steiner points
        const std::string switches = "zpq" + std::to_string(radius_edge_ratio) + "a" + std::to_string(max_tet_vol) + "V";
        igl::copyleft::tetgen::tetrahedralize(V0, F, switches, nodes, elem, faces);
        const auto AboveBM = [&bm_v, &bm_f, &nodes](int p3) -> bool {
            Eigen::VectorXd wnumber;
            igl::winding_number(bm_v, bm_f, nodes.row(p3), wnumber);
            return wnumber(0)>0;
        };

        // DEBUG
        std::vector<int> ss;

        Eigen::MatrixXd result_v = nodes;
        Eigen::MatrixXi result_f(elem.rows()*4, 3);
        const auto SaveMsh = [&bm_v, &bm_f, &ss, &result_f](const char* fileName, Eigen::MatrixXd &V, Eigen::MatrixXi &T) {
            H5Easy::File file("./" + std::string(fileName), H5Easy::File::Overwrite);
            H5Easy::dump(file, "V", V);
            H5Easy::dump(file, "T", T);

            H5Easy::dump(file, "bm_v", bm_v);
            H5Easy::dump(file, "bm_f", bm_f);
            H5Easy::dump(file, "Tid", ss);
            H5Easy::dump(file, "result_f", result_f);
        };

        // after TetGen {bm_v, bm_f} -> {result_v, result_f}
        std::vector<bool> on_bm_surface(nodes.rows(), false);
        Eigen::VectorXd squaredDist;
        Eigen::MatrixXd I, C;
        igl::point_mesh_squared_distance(nodes, bm_v, bm_f, squaredDist, I, C);
        const double thres = 1e-14;
        for (int i=0; i<nodes.rows(); i++) {
            if (squaredDist(i) < thres) on_bm_surface[i] = true;
        }
        int cntFaces = 0, cntInvert = 0;
        for (int i=0; i<elem.rows(); i++) {
            for (int j=0; j<4; j++) {
                int f0 = elem(i, (j+1)%4);
                int f1 = elem(i, (j+2)%4);
                int f2 = elem(i, (j+3)%4);
                if (on_bm_surface[f0] && on_bm_surface[f1] && on_bm_surface[f2]) {
                    // make sure this face is on bm
                    Eigen::MatrixXi f(1, 3);
                    f << f0, f1, f2;
                    Eigen::MatrixXd bc;
                    igl::barycenter(nodes, f, bc);
                    igl::point_mesh_squared_distance(bc, bm_v, bm_f, squaredDist, I, C);
                    if (squaredDist(0) > thres) continue;  // a pseudo on-bm face, a bridge!
                    // orientation matters here
                    if (AboveBM(elem(i, j))) {
                        result_f.row(cntFaces) << f0, f1, f2;
                        cntFaces++;
                        ss.push_back(i);  // may have duplicated tet indices
                    } else {
                        cntInvert++;
                    }
                }
            }
        }
        result_f.conservativeResize(cntFaces, 3);
        // should clean mesh here TODO
        // log
        int cntV = 0;
        for (int i=0; i<on_bm_surface.size(); i++)
            if (on_bm_surface[i]) cntV++;
        std::cout << "bm_v.size = " << bm_v.rows() << " cntV = " << cntV << std::endl;
        std::cout << "bm_f.size = " << bm_f.rows() << " result_f.size = " << cntFaces << " | invert = " << cntInvert << std::endl;

        // DEBUG
        SaveMsh("test.h5", nodes, elem);

        // assert bm surface area
        double bm_area = 0, result_area = 0;
        Eigen::VectorXd areaV;
        igl::doublearea(bm_v, bm_f, areaV);
        bm_area = areaV.sum();
        igl::doublearea(result_v, result_f, areaV);
        result_area = areaV.sum();
        if (bm_area != result_area) {
            std::cerr << "bm_area = " << bm_area << " result_area = " << result_area << std::endl;
        }
        

        // lamda for bc
        Eigen::MatrixXd barys;
        igl::barycenter(result_v, result_f, barys);

        // get the bounding box for the first V
        const Eigen::MatrixXd box_min = V0.colwise().minCoeff();
        const Eigen::MatrixXd box_max = V0.colwise().maxCoeff();
        assert(box_min.size() == 3);
        assert(box_max.size() == 3);

        const auto set_bc = [&barys, &box_min, &box_max](const Eigen::MatrixXd &v, const bool boundary) {
            if (boundary)
            {
                if (fabs(v(2) - box_min(2)) < 0.1)
                    return 1;
                if (fabs(v(2) - box_max(2)) < 0.1)
                    return 1;  // change me to 0 if free to move

                return 0;
            }
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

            if (min < 0.1)
                return 2;

            return 0;
        };

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

        {
            Eigen::MatrixXd points;
            Eigen::MatrixXi faces;
            Eigen::MatrixXd sidesets;

            state.get_sidesets(points, faces, sidesets);
            igl::write_triangle_mesh("text.obj", points, faces);
            igl::write_triangle_mesh("text1.obj", V0, F);
        }

        const Eigen::RowVector3d zero(0, 0, 0);
        Eigen::MatrixXd vals;
        Eigen::MatrixXd traction_forces, traction_forces_flip;
        Eigen::MatrixXd stress, stress_flip;
        Eigen::MatrixXd mises, mises_flip;
        Eigen::MatrixXi result_f_flip = result_f;
        result_f_flip.col(0) = result_f.col(1);
        result_f_flip.col(1) = result_f.col(0);

        for (int sim = 1; sim < V.size(); ++sim)
        {
            const json args = {
                {"tensor_formulation", is_linear ? "LinearElasticity" : "NeoHookean"},
                {"discr_order", discr_order},

                {"params", {{"E", E / double(1e6)}, {"nu", nu}}}, // E has unit [Pascal], convert that to [um] unit

                {"problem", "PointBasedTensor"}};

            const Eigen::MatrixXd currentv = V[sim].block(0, 0, bm_offset_v, 3) - bm_v;
            // Eigen::MatrixXd currentv = V[sim] - ref_v;
            // currentv.setConstant(.5);

            state.init(args);
            polyfem::PointBasedTensorProblem &tproblem = *dynamic_cast<polyfem::PointBasedTensorProblem *>(state.problem.get());
            tproblem.add_constant(1, zero);
            Eigen::Matrix<bool, 3, 1> dirichet_dims;
            dirichet_dims.setConstant(true);
            const Eigen::MatrixXd twod = bm_v.block(0, 0, bm_v.rows(), 2);
            tproblem.add_function(2, currentv, twod, rbf_function, eps, 2, dirichet_dims);
            state.solve();

            //output
            state.interpolate_boundary_function_at_vertices(result_v, result_f, state.sol, vals);
            state.interpolate_boundary_tensor_function(result_v, result_f, state.sol, vals, true, traction_forces, stress, mises);
            state.interpolate_boundary_tensor_function(result_v, result_f_flip, state.sol, vals, true, traction_forces_flip, stress_flip, mises_flip);

            Eigen::MatrixXd vertex_traction_forces(result_v.rows(), 3);
            vertex_traction_forces.setZero();

            Eigen::MatrixXd vertex_stress(result_v.rows(), 9);
            vertex_stress.setZero();

            Eigen::MatrixXd vertex_mises(result_v.rows(), 1);
            vertex_mises.setZero();

            Eigen::MatrixXd vertex_traction_forces_flip(result_v.rows(), 3);
            vertex_traction_forces_flip.setZero();

            Eigen::MatrixXd vertex_stress_flip(result_v.rows(), 9);
            vertex_stress_flip.setZero();

            Eigen::MatrixXd vertex_mises_flip(result_v.rows(), 1);
            vertex_mises_flip.setZero();

            Eigen::MatrixXd area, vertex_area(result_v.rows(), 1);
            vertex_area.setZero();
            igl::doublearea(result_v, result_f, area);

            for (int f = 0; f < result_f.rows(); ++f)
            {
                for (int d = 0; d < 3; ++d)
                {
                    vertex_traction_forces.row(result_f(f, d)) += traction_forces.row(f) * area(f);
                    vertex_traction_forces_flip.row(result_f(f, d)) += traction_forces_flip.row(f) * area(f);

                    vertex_stress.row(result_f(f, d)) += stress.row(f) * area(f);
                    vertex_stress_flip.row(result_f(f, d)) += stress_flip.row(f) * area(f);

                    vertex_mises(result_f(f, d)) += mises(f) * area(f);
                    vertex_mises_flip(result_f(f, d)) += mises_flip(f) * area(f);

                    vertex_area(result_f(f, d)) += area(f);
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

            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    const std::string index = "_" + std::to_string(i) + std::to_string(j);
                    writer.add_field("stress" + index, vertex_stress.col(i * 3 + j));
                    writer.add_field("stress_other" + index, vertex_stress_flip.col(i * 3 + j));
                }
            }

            writer.add_field("von_mises", vertex_mises);
            writer.add_field("von_mises_other", vertex_mises_flip);

            writer.write_mesh(out_path, result_v, result_f);

            state.save_vtu(out_path + ".all.vtu", 0);
        }
    }
} // namespace zebrafish
