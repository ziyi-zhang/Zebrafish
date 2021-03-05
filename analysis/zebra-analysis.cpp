#include <zebrafish/zebra-analysis.hpp>

#include <polyfem/GenericProblem.hpp>
#include <polyfem/PointBasedProblem.hpp>
#include <polyfem/State.hpp>
#include <polyfem/VTUWriter.hpp>

#include <highfive/H5Easy.hpp>
#include <igl/barycenter.h>
#include <igl/boundary_loop.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/doublearea.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/remove_unreferenced.h>
#include <igl/upsample.h>
#include <igl/winding_number.h>
#include <igl/writeMESH.h>
#include <igl/write_triangle_mesh.h>

#include <algorithm>

namespace zebrafish
{
    bool compute_analysis(const std::vector<Eigen::MatrixXd> &VV, const Eigen::MatrixXi &F,
                          const int bm_offset_v, const int bm_offset_f,
                          const std::string &path,
                          const double E, const double nu,
                          const double max_tet_vol,
                          const int discr_order, const bool is_linear, const int n_refs, const double vismesh_rel_area, const int upsample,
                          const std::map<int, std::array<int, 2>> &markerRCMap, const int imgRows, const int imgCols, const int layerPerImg,
                          const double resolutionX, const double resolutionY, const double resolutionZ,
                          const bool saveinput, const bool aboveCage, const bool belowCage,
                          const int frameStart, const int frameEnd)
    {
        // Necesary for the scaling!
        std::vector<Eigen::MatrixXd> V = VV;

        // RBF setup
        //const std::string rbf_function = "gaussian";
        const std::string rbf_function = "thin-plate";
        // const std::string rbf_function = "cubic";
        const double eps = 10;

        // save to analysis file for future re-analysis
        const auto WriteInputToFile = [&V, &F, &bm_offset_v, &bm_offset_f, &path, &E, &nu, &max_tet_vol, &discr_order, &is_linear, &n_refs, &vismesh_rel_area, &upsample,
                                       &markerRCMap, &imgRows, &imgCols, &layerPerImg, &resolutionX, &resolutionY, &resolutionZ]() {
            // NOTE: only store the part of {V, F} without the cage

            const int frames = V.size();
            const int Nverts = bm_offset_v;

            Eigen::MatrixXd V_concatenated(frames * Nverts, 3);
            for (int i = 0; i < frames; i++)
            {
                V_concatenated.block(Nverts * i, 0, Nverts, 3) = V[i].topRows(Nverts);
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
            H5Easy::dump(file, "F", F.topRows(bm_offset_f));
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

        /////////////////////////////////////////////////////////////////////////////////////////////
        // Generate tet mesh
        Eigen::MatrixXd nodes;
        Eigen::MatrixXi elem, faces;
        // allow exterior steiner points
        // const std::string switches = "zpq" + std::to_string(radius_edge_ratio) + "a" + std::to_string(max_tet_vol) + "V";
        std::stringstream buf;
        buf.precision(100);
        buf.setf(std::ios::fixed, std::ios::floatfield);
        buf << "Qpq1.414a" << max_tet_vol;
        igl::copyleft::tetgen::tetrahedralize(V0, F, buf.str(), nodes, elem, faces);

        // igl::writeMESH("out.mesh", nodes, elem, faces);

        Eigen::MatrixXd tetMeshBC;
        igl::barycenter(nodes, elem, tetMeshBC);
        const auto AboveBM = [&bm_v, &bm_f, &tetMeshBC, &aboveCage](int i) -> bool {
            Eigen::VectorXd wnumber;
            igl::winding_number(bm_v, bm_f, tetMeshBC.row(i), wnumber);
            if (aboveCage)
                return wnumber(0) < 0;
            else
                return wnumber(0) > 0;
        };
        const auto NormalPointingUp = [&nodes](int p0, int p1, int p2) -> bool {
            Eigen::Vector3d v01, v02, v_cross;
            v01 = nodes.row(p1) - nodes.row(p0);
            v02 = nodes.row(p2) - nodes.row(p0);
            v_cross = v01.cross(v02);
            return v_cross(2) > 0; // cannot == 0
        };

        // DEBUG
        std::vector<int> ss;

        // after TetGen {bm_v, bm_f} -> {result_v, result_f}
        Eigen::MatrixXd result_v = nodes;
        Eigen::MatrixXi result_f(elem.rows() * 4, 3);
        std::vector<bool> on_bm_surface(nodes.rows(), false);
        Eigen::VectorXd squaredDist;
        Eigen::MatrixXd I, C;
        igl::point_mesh_squared_distance(nodes, bm_v, bm_f, squaredDist, I, C);
        const double thres = 1e-14;
        for (int i = 0; i < nodes.rows(); i++)
        {
            if (squaredDist(i) < thres)
            {
                on_bm_surface[i] = true;
            }
        }

        int cntFaces = 0, cntInvert = 0;
        for (int i = 0; i < elem.rows(); i++)
        {
            for (int j = 0; j < 4; j++)
            {
                int f0 = elem(i, (j + 1) % 4);
                int f1 = elem(i, (j + 2) % 4);
                int f2 = elem(i, (j + 3) % 4);
                if (on_bm_surface[f0] && on_bm_surface[f1] && on_bm_surface[f2])
                {
                    // // make sure this face is on bm
                    Eigen::MatrixXi f(1, 3);
                    f << f0, f1, f2;
                    Eigen::MatrixXd bc;
                    igl::barycenter(nodes, f, bc);
                    igl::point_mesh_squared_distance(bc, bm_v, bm_f, squaredDist, I, C);
                    if (squaredDist(0) > thres)
                        continue; // a pseudo on-bm face, a bridge!
                    // orientation matters here
                    if (AboveBM(i))
                    {
                        if (NormalPointingUp(f0, f1, f2))
                            result_f.row(cntFaces) << f0, f1, f2;
                        else
                            result_f.row(cntFaces) << f1, f0, f2;

                        cntFaces++;
                        ss.push_back(i); // may have duplicated tet indices
                    }
                    else
                    {
                        cntInvert++;
                    }
                }
            }
        }
        result_f.conservativeResize(cntFaces, 3);

        // remove unreferenced V
        Eigen::MatrixXi Im, result_f_ = result_f;
        Eigen::MatrixXd result_v_ = result_v;
        igl::remove_unreferenced(result_v_, result_f_, result_v, result_f, Im);

        // log
        std::cout << "bm_v.size = " << bm_v.rows() << " result_v.size = " << result_v.rows() << std::endl;
        std::cout << "bm_f.size = " << bm_f.rows() << " result_f.size = " << cntFaces << " | invert.size = " << cntInvert << std::endl;

        // DEBUG
        const auto SaveMsh = [&bm_v, &bm_f, &ss, &result_f](const char *fileName, Eigen::MatrixXd &V, Eigen::MatrixXi &T) {
            H5Easy::File file("./" + std::string(fileName), H5Easy::File::Overwrite);
            H5Easy::dump(file, "V", V);
            H5Easy::dump(file, "T", T);

            H5Easy::dump(file, "bm_v", bm_v);
            H5Easy::dump(file, "bm_f", bm_f);
            H5Easy::dump(file, "Tid", ss);
            H5Easy::dump(file, "result_f", result_f);
        };
        // SaveMsh("test.h5", nodes, elem);
        // igl::write_triangle_mesh("extracted.obj", result_v, result_f);

        // assert bm surface area
        double bm_area = 0, result_area = 0;
        Eigen::VectorXd areaV;
        igl::doublearea(bm_v, bm_f, areaV);
        bm_area = areaV.sum();
        igl::doublearea(result_v, result_f, areaV);
        result_area = areaV.sum();
        if (std::fabs(bm_area - result_area) > 1e-10)
        {
            std::cerr << "[WARNING] bm_area = " << bm_area << " result_area = " << result_area << std::endl;
            return false;
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

            if (min < 0.01)
                return 2;

            if (boundary)
            {
                if (fabs(v(2) - box_min(2)) < 0.01)
                    return 1;
                if (fabs(v(2) - box_max(2)) < 0.01)
                    return 1; // change me to 0 if free to move

                return 0;
            }

            return 0;
        };

        /////////////////////////////////////////////////////////////////////////////////////////////
        // polyfem
        polyfem::State state;
        state.mesh = std::make_unique<polyfem::Mesh3D>();
        state.mesh->build_from_matrices(nodes, elem);
        //igl::writeMESH("out.mesh", nodes, elem, faces);
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

        const Eigen::RowVector3d zero(0, 0, 0);

        int fStart = frameStart;
        int fEnd = frameEnd;
        if (frameStart == -1)
            fStart = 1;
        if (frameEnd == -1)
            fEnd = V.size() - 1;
        for (int sim = fStart; sim <= fEnd; ++sim)
        {
            const json args = {
                {"tensor_formulation", is_linear ? "LinearElasticity" : "NeoHookean"},
                {"discr_order", discr_order},
                {"params", {{"E", E / double(1e6)}, {"nu", nu}}}, // E has unit [Pascal], convert that to [um] unit
                {"problem", "PointBasedTensor"}};

            const Eigen::MatrixXd currentv = V[sim].block(0, 0, bm_offset_v, 3) - bm_v;

            state.init(args);
            polyfem::PointBasedTensorProblem &tproblem = *dynamic_cast<polyfem::PointBasedTensorProblem *>(state.problem.get());
            tproblem.add_constant(1, zero);
            Eigen::Matrix<bool, 3, 1> dirichet_dims;
            dirichet_dims.setConstant(true);
            tproblem.add_function(2, currentv, bm_v, rbf_function, eps, -1, dirichet_dims);
            state.solve();

            // output
            const auto OutputHelper = [&state, &belowCage, &aboveCage](const std::string &out_path, const Eigen::MatrixXd &mesh_v, const MatrixXi &mesh_f) {
                Eigen::MatrixXd displacement_vec;
                Eigen::MatrixXd traction_forces, traction_forces_flip;
                Eigen::MatrixXd stress, stress_flip;
                Eigen::MatrixXd mises, mises_flip;
                Eigen::MatrixXi mesh_f_flip = mesh_f;
                mesh_f_flip.col(0) = mesh_f.col(1);
                mesh_f_flip.col(1) = mesh_f.col(0);

                state.interpolate_boundary_function_at_vertices(mesh_v, mesh_f, state.sol, displacement_vec);
                // guaranteed that (belowCage || aboveCage) == True
                if (belowCage)
                    state.interpolate_boundary_tensor_function(mesh_v, mesh_f, state.sol, displacement_vec, true, traction_forces, stress, mises);
                if (aboveCage)
                    state.interpolate_boundary_tensor_function(mesh_v, mesh_f_flip, state.sol, displacement_vec, true, traction_forces_flip, stress_flip, mises_flip);

                // per-triangle data -> per-vertex data by taking average
                const auto ToVertexData = [&mesh_v, &mesh_f](
                                              const Eigen::MatrixXd &traction_forces_tri,
                                              const Eigen::MatrixXd &stress_tri,
                                              const Eigen::MatrixXd &mises_tri,
                                              Eigen::MatrixXd &traction_forces_ver,
                                              Eigen::MatrixXd &stress_ver,
                                              Eigen::MatrixXd &mises_ver) {
                    traction_forces_ver = Eigen::MatrixXd::Zero(mesh_v.rows(), 3);
                    stress_ver = Eigen::MatrixXd::Zero(mesh_v.rows(), 9); // 3x3 matrix
                    mises_ver = Eigen::MatrixXd::Zero(mesh_v.rows(), 1);

                    Eigen::VectorXd area, vertex_area(mesh_v.rows());
                    vertex_area.setZero();
                    igl::doublearea(mesh_v, mesh_f, area);

                    for (int f = 0; f < mesh_f.rows(); ++f)
                    {
                        for (int d = 0; d < 3; ++d)
                        {
                            int vid = mesh_f(f, d);
                            traction_forces_ver.row(vid) += traction_forces_tri.row(f) * area(f);
                            stress_ver.row(vid) += stress_tri.row(f) * area(f);
                            mises_ver(vid) += mises_tri(f) * area(f);
                            vertex_area(vid) += area(f);
                        }
                    }

                    for (int d = 0; d < 3; ++d)
                    {
                        traction_forces_ver.col(d).array() /= vertex_area.array();
                    }
                    for (int d = 0; d < 9; ++d)
                    {
                        stress_ver.col(d).array() /= vertex_area.array();
                    }
                    mises_ver.array() /= vertex_area.array();
                };

                // triangle to vertex
                Eigen::MatrixXd traction_force_v, traction_force_v_flip, traction_force_res;
                Eigen::MatrixXd stress_v, stress_v_flip, stress_res;
                Eigen::MatrixXd mises_v, mises_v_flip, mises_res;
                if (belowCage)
                    ToVertexData(traction_forces, stress, mises, traction_force_v, stress_v, mises_v);
                if (aboveCage)
                    ToVertexData(traction_forces_flip, stress_flip, mises_flip, traction_force_v_flip, stress_v_flip, mises_v_flip);
                std::string suffix = "";
                if (aboveCage && !belowCage)
                {
                    traction_force_res = traction_force_v_flip;
                    stress_res = stress_v_flip;
                    mises_res = mises_v_flip;
                    suffix = "_above";
                }
                else if (!aboveCage && belowCage)
                {
                    traction_force_res = traction_force_v;
                    stress_res = stress_v;
                    mises_res = mises_v;
                    suffix = "_below";
                }
                else if (aboveCage && belowCage)
                {
                    traction_force_res = (traction_force_v_flip + traction_force_v).array() / 2.0;
                    stress_res = (stress_v_flip + stress_v).array() / 2.0;
                    mises_res = (mises_v_flip + mises_v).array() / 2.0;
                    suffix = "_avg";
                }
                else
                {
                    std::cerr << "neither aboveCage or belowCage? impossible!" << std::endl;
                    return;
                }

                // write to VTU
                polyfem::VTUWriter VTUwriter;
                VTUwriter.add_field("displacement", displacement_vec);
                VTUwriter.add_field("traction_forces" + suffix, traction_force_res);

                for (int i = 0; i < 3; ++i)
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        const std::string index = "_" + std::to_string(i) + std::to_string(j);
                        VTUwriter.add_field("stress" + suffix + index, stress_res.col(i * 3 + j));
                    }
                }

                VTUwriter.add_field("von_mises" + suffix, mises_res);

                VTUwriter.write_mesh(out_path, mesh_v, mesh_f);
            }; // OutputHelper lambda function

            const auto RemoveOneRing = [&result_f](Eigen::MatrixXi &new_f) {
                std::vector<int> boundary_v;
                igl::boundary_loop(result_f, boundary_v);
                new_f.resizeLike(result_f);
                int cnt = 0;
                for (int i = 0; i < result_f.rows(); i++)
                {
                    if (std::find(boundary_v.begin(), boundary_v.end(), result_f(i, 0)) != boundary_v.end() ||
                        std::find(boundary_v.begin(), boundary_v.end(), result_f(i, 1)) != boundary_v.end() ||
                        std::find(boundary_v.begin(), boundary_v.end(), result_f(i, 2)) != boundary_v.end())
                    {
                        // if this triangle is at boundary
                        continue;
                    }
                    else
                    {
                        new_f.row(cnt) = result_f.row(i);
                        cnt++;
                    }
                }
                new_f.conservativeResize(cnt, 3);
            }; // KillOneRing lambda function

            const std::string out_path = path + "-frame" + std::to_string(sim);
            Eigen::MatrixXi result_f_without_boundary;
            RemoveOneRing(result_f_without_boundary);
            OutputHelper(out_path + ".withBoundary.vtu", result_v, result_f);
            OutputHelper(out_path + ".vtu", result_v, result_f_without_boundary);
            state.save_vtu(out_path + ".all.vtu", 0);
            state.save_surface(out_path + ".surf.vtu");

            std::cout << out_path << "  analysis and export done!" << std::endl;
        }

        return true;
    }
} // namespace zebrafish
