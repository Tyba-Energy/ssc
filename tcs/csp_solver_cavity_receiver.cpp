/**
BSD-3-Clause
Copyright 2019 Alliance for Sustainable Energy, LLC
Redistribution and use in source and binary forms, with or without modification, are permitted provided
that the following conditions are met :
1.	Redistributions of source code must retain the above copyright notice, this list of conditions
and the following disclaimer.
2.	Redistributions in binary form must reproduce the above copyright notice, this list of conditions
and the following disclaimer in the documentation and/or other materials provided with the distribution.
3.	Neither the name of the copyright holder nor the names of its contributors may be used to endorse
or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT HOLDER, CONTRIBUTORS, UNITED STATES GOVERNMENT OR UNITED STATES
DEPARTMENT OF ENERGY, NOR ANY OF THEIR EMPLOYEES, BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
OR CONSEQUENTIAL DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "csp_solver_cavity_receiver.h"
#include "csp_solver_core.h"
#include "csp_solver_util.h"
#include "sam_csp_util.h"

#include "Ambient.h"
#include "definitions.h"
#include <math.h>

#include "../splinter/Core"
#include "../splinter/LU"
#include "../splinter/Cholesky"
#include "../splinter/QR"
#include "../splinter/SVD"
#include "../splinter/Geometry"
#include "../splinter/Eigenvalues"

bool sort_pair_ascending(pair<double,double> i, pair<double, double> j)
{
    if (i.first > j.first) {
        return true;
    }
    else if (i.first == j.first && i.second > j.second) {
        return true;
    }
    else {
        return false;
    }
}

C_cavity_receiver::C_cavity_receiver(double dni_des /*W/m2*/, double hel_stow_deploy /*-*/,
    int field_fl /*-*/, util::matrix_t<double> field_fl_props,
    double od_rec_tube /*m*/, double th_rec_tube /*m*/, int tube_mat_code /*-*/,
    double rec_height /*m*/, double rec_width /*m*/, double toplip_height /*m*/, double botlip_height /*m*/,
    double eps_active_sol /*-*/, double eps_passive_sol /*-*/, double eps_active_therm /*-*/, double eps_passive_therm /*-*/,
    double pipe_loss_per_m /*Wt/m*/, double pipe_length_add /*m*/, double pipe_length_mult /*-*/)
{
    m_dni_des = dni_des;                    //[W/m2]
    m_hel_stow_deploy = hel_stow_deploy;    //[deg]
    m_field_fl = field_fl;                  //[-]
    m_field_fl_props = field_fl_props;      //[-]

    m_od_rec_tube = od_rec_tube;            //[m]
    m_th_rec_tube = th_rec_tube;            //[m]
    m_tube_mat_code = tube_mat_code;        //[-]

    m_receiverHeight = rec_height;          //[m]
    m_receiverWidth = rec_width;            //[m]
    m_topLipHeight = toplip_height;         //[m]
    m_botLipHeight = botlip_height;         //[m]
    m_e_act_sol = eps_active_sol;           //[-]
    m_e_pass_sol = eps_passive_sol;         //[-]
    m_e_act_therm = eps_active_therm;       //[-]
    m_e_pass_therm = eps_passive_therm;     //[-]

    m_pipe_loss_per_m = pipe_loss_per_m;    //[Wt/m]
    m_pipe_length_add = pipe_length_add;    //[m]
    m_pipe_length_mult = pipe_length_mult;  //[-]

    m_area_active_total = std::numeric_limits<double>::quiet_NaN();
    m_d_in_rec_tube = std::numeric_limits<double>::quiet_NaN();
    m_A_cs_tube = std::numeric_limits<double>::quiet_NaN();
    m_Ntubes = 0;
    m_Q_dot_piping_loss = std::numeric_limits<double>::quiet_NaN();
    m_rel_roughness = std::numeric_limits<double>::quiet_NaN();
}

void C_cavity_receiver::genOctCavity()
{
    // Matlab function:
    //function[PANEL1, PANEL2, PANEL3, PANEL4, FLOOR, COVER, TOPLIP, BOTLIP, APERTURE] = genOctCavity(height, width, lipTop, lipBot, varargin) % #codegen
    
    // Creates geometry (i.e. defines vertices) for a 4-panel half-octagonal cavity receiver
    // of specified width, height, and with or without upper and lower lips
    mv_rec_surfs.resize(m_nPanels + 5);

    mv_rec_surfs[PANEL1].type = 0; // 1;
    mv_rec_surfs[PANEL2].type = 0; // 1;
    mv_rec_surfs[PANEL3].type = 0; // 1;
    mv_rec_surfs[PANEL4].type = 0; // 1;
    mv_rec_surfs[FLOOR].type = 3; // 2;
    mv_rec_surfs[COVER].type = 3; // 2;
    mv_rec_surfs[TOPLIP].type = 0; // 1;
    mv_rec_surfs[BOTLIP].type = 0; // 1;
    mv_rec_surfs[APERTURE].type = 3; // 2;

    mv_rec_surfs[PANEL1].is_active_surf = true;
    mv_rec_surfs[PANEL2].is_active_surf = true;
    mv_rec_surfs[PANEL3].is_active_surf = true;
    mv_rec_surfs[PANEL4].is_active_surf = true;
    mv_rec_surfs[FLOOR].is_active_surf = false;
    mv_rec_surfs[COVER].is_active_surf = false;
    mv_rec_surfs[TOPLIP].is_active_surf = false;
    mv_rec_surfs[BOTLIP].is_active_surf = false;
    mv_rec_surfs[APERTURE].is_active_surf = false;

    mv_rec_surfs[PANEL1].is_flipRoute = true;
    mv_rec_surfs[PANEL2].is_flipRoute = false;
    mv_rec_surfs[PANEL3].is_flipRoute = true;
    mv_rec_surfs[PANEL4].is_flipRoute = false;

    size_t n_surfs = mv_rec_surfs.size();
    for (size_t i = 0; i < n_surfs; i++) {
        if (i == APERTURE) {
            mv_rec_surfs[i].eps_sol = 1.0;
            mv_rec_surfs[i].eps_therm = 1.0;
        }
        else {
            if (mv_rec_surfs[i].is_active_surf) {
                mv_rec_surfs[i].eps_sol = m_e_act_sol;        //[-]
                mv_rec_surfs[i].eps_therm = m_e_act_therm;    //[-]
            }
            else {
                mv_rec_surfs[i].eps_sol = m_e_pass_sol;       //[-]
                mv_rec_surfs[i].eps_therm = m_e_pass_therm;   //[-]
            }
        }
    }

    // Polygon interior angle
    // Assumes center of polygon is center of aperture
    double theta = CSP::pi / (double)m_nPanels;   //[rad]

    // Calculate element sizes for mesh based on dimensions and model resolution parameter
    double panelWidth = 2.0*(m_receiverWidth/2.0)*sin(theta/2.0);    //[m]
    mv_rec_surfs[PANEL1].surf_elem_size = panelWidth / (double)m_pipeWindings / (double)m_modelRes;
    mv_rec_surfs[PANEL4].surf_elem_size = mv_rec_surfs[PANEL3].surf_elem_size = mv_rec_surfs[PANEL2].surf_elem_size = mv_rec_surfs[PANEL1].surf_elem_size;
    mv_rec_surfs[FLOOR].surf_elem_size = mv_rec_surfs[COVER].surf_elem_size = m_receiverWidth / 4.0 / m_modelRes;
    mv_rec_surfs[TOPLIP].surf_elem_size = mv_rec_surfs[BOTLIP].surf_elem_size = m_receiverWidth / 6.0 / m_modelRes;

    // matrix_t(nr, nc, val)
    // (x, y, z)
    mv_rec_surfs[FLOOR].vertices.resize_fill(m_nPanels + 1, 3, 0.0);
    mv_rec_surfs[COVER].vertices.resize_fill(m_nPanels + 1, 3, 0.0);

    for (size_t i = 0; i < m_nPanels+1; i++) {
        mv_rec_surfs[FLOOR].vertices(i, 0) = mv_rec_surfs[COVER].vertices(i, 0) = m_receiverWidth * cos(i * theta)/2.0;
        mv_rec_surfs[FLOOR].vertices(i, 1) = mv_rec_surfs[COVER].vertices(i, 1) = m_receiverWidth * sin(i * theta)/2.0;
        mv_rec_surfs[FLOOR].vertices(i, 2) = -m_receiverHeight / 2.0;
        mv_rec_surfs[COVER].vertices(i, 2) = m_receiverHeight / 2.0;
    }

    util::matrix_t<double> temp_total(4, 3, std::numeric_limits<double>::quiet_NaN());
    util::matrix_t<double> temp_p1(4, 3, std::numeric_limits<double>::quiet_NaN());

    for (size_t n = 0; n < m_nPanels; n++) {
        temp_p1 = mv_rec_surfs[FLOOR].vertices.row(n);
        for (size_t i = 0; i < temp_p1.length(); i++) {
            temp_total(0, i) = temp_p1(0, i);
        }
        temp_p1 = mv_rec_surfs[COVER].vertices.row(n);
        for (size_t i = 0; i < temp_p1.length(); i++) {
            temp_total(1, i) = temp_p1(0, i);
        }
        temp_p1 = mv_rec_surfs[COVER].vertices.row(n + 1);
        for (size_t i = 0; i < temp_p1.length(); i++) {
            temp_total(2, i) = temp_p1(0, i);
        }
        temp_p1 = mv_rec_surfs[FLOOR].vertices.row(n + 1);
        for (size_t i = 0; i < temp_p1.length(); i++) {
            temp_total(3, i) = temp_p1(0, i);
        }
        if (n == 0) {
            mv_rec_surfs[PANEL1].vertices = temp_total;
        }
        else if (n == 1) {
            mv_rec_surfs[PANEL2].vertices = temp_total;
        }
        else if (n == 2) {
            mv_rec_surfs[PANEL3].vertices = temp_total;
        }
        else if (n == 3) {
            mv_rec_surfs[PANEL4].vertices = temp_total;
        }
    }

    mv_rec_surfs[APERTURE].vertices.resize_fill(4, 3, std::numeric_limits<double>::quiet_NaN());
    mv_rec_surfs[TOPLIP].vertices.resize_fill(4, 3, std::numeric_limits<double>::quiet_NaN());
    mv_rec_surfs[BOTLIP].vertices.resize_fill(4, 3, std::numeric_limits<double>::quiet_NaN());

    if (m_topLipHeight <= 0.0) {
        temp_p1 = mv_rec_surfs[COVER].vertices.row(0);
        for (size_t i = 0; i < temp_p1.length(); i++) {
            mv_rec_surfs[APERTURE].vertices(0, i) = temp_p1(0, i);
        }
        temp_p1 = mv_rec_surfs[COVER].vertices.row(m_nPanels);
        for (size_t i = 0; i < temp_p1.length(); i++) {
            mv_rec_surfs[APERTURE].vertices(1, i) = temp_p1(0, i);
        }
    }
    else {
        util::matrix_t<double> temp_a(2, 3, std::numeric_limits<double>::quiet_NaN());
        for (size_t j = 0; j < 3; j++) {
            mv_rec_surfs[TOPLIP].vertices(0,j) = mv_rec_surfs[COVER].vertices(0,j);
            mv_rec_surfs[TOPLIP].vertices(1,j) = mv_rec_surfs[COVER].vertices(mv_rec_surfs[COVER].vertices.nrows()-1,j);
        }
        for (size_t i = 0; i < 2; i++) {
            temp_a(i, 0) = mv_rec_surfs[TOPLIP].vertices(i, 0);
            temp_a(i, 1) = mv_rec_surfs[TOPLIP].vertices(i, 1);
            temp_a(i, 2) = mv_rec_surfs[TOPLIP].vertices(i, 2) - m_topLipHeight;
        }
        util::matrix_t<double> temp_b;
        flipup(temp_a, temp_b);
        for (size_t j = 0; j < 3; j++) {
            mv_rec_surfs[TOPLIP].vertices(2,j) = temp_b(0,j);
            mv_rec_surfs[TOPLIP].vertices(3,j) = temp_b(1,j);
            mv_rec_surfs[APERTURE].vertices(0,j) = mv_rec_surfs[TOPLIP].vertices(3,j);
            mv_rec_surfs[APERTURE].vertices(1,j) = mv_rec_surfs[TOPLIP].vertices(2,j);
        }
    }

    if (m_botLipHeight <= 0.0) {
        temp_p1 = mv_rec_surfs[FLOOR].vertices.row(m_nPanels);
        for (size_t i = 0; i < temp_p1.length(); i++) {
            mv_rec_surfs[APERTURE].vertices(2, i) = temp_p1(0, i);
        }
        temp_p1 = mv_rec_surfs[FLOOR].vertices.row(0);
        for (size_t i = 0; i < temp_p1.length(); i++) {
            mv_rec_surfs[APERTURE].vertices(3, i) = temp_p1(0, i);
        }
    }
    else {
        util::matrix_t<double> temp_a(2,3,std::numeric_limits<double>::quiet_NaN());
        for (size_t j = 0; j < 3; j++) {
            mv_rec_surfs[BOTLIP].vertices(0,j) = mv_rec_surfs[FLOOR].vertices(0,j);
            mv_rec_surfs[BOTLIP].vertices(1,j) = mv_rec_surfs[FLOOR].vertices(mv_rec_surfs[FLOOR].vertices.nrows()-1,j);
        }
        for (size_t i = 0; i < 2; i++) {
            temp_a(i,0) = mv_rec_surfs[BOTLIP].vertices(i,0);
            temp_a(i,1) = mv_rec_surfs[BOTLIP].vertices(i,1);
            temp_a(i,2) = mv_rec_surfs[BOTLIP].vertices(i,2) + m_botLipHeight;
        }
        util::matrix_t<double> temp_b;
        flipup(temp_a, temp_b);
        for (size_t j = 0; j < 3; j++) {
            mv_rec_surfs[BOTLIP].vertices(2,j) = temp_b(0,j);
            mv_rec_surfs[BOTLIP].vertices(3,j) = temp_b(1,j);
            mv_rec_surfs[APERTURE].vertices(2,j) = mv_rec_surfs[BOTLIP].vertices(2,j);
            mv_rec_surfs[APERTURE].vertices(3,j) = mv_rec_surfs[BOTLIP].vertices(3,j);
        }
    }

    return;
}


void C_cavity_receiver::meshGeometry()
{
    /*% Mesh Surfaces
        [nodesGlobal, elemsPanel1, elemsPanel2, elemsPanel3, elemsPanel4, ...
        elemsFloor, elemsCover, elemsTopLip, elemsBotLip, elemsApert] = meshGeometry(elemSizes, meshTypes, ...
            PANEL1, PANEL2, PANEL3, PANEL4, FLOOR, COVER, TOPLIP, BOTLIP, APERTURE);*/

    // Meshes each surface, depending on each surface "type"
    // Also uses the same constant 'm_elemSize' for each surface
    // Defines member data:
    // -- m_v_elems: each vector index is a surface; each row lists the nodes that define a mesh element
    // -- m_nodesGlobal: each row lists x,y,z coordinates of each node

    vector<util::matrix_t<double>> v_nodes;

    int lastNodeID = 0;

    for (size_t i = 0; i < mv_rec_surfs.size(); i++)
    {
        util::matrix_t<double> surf = mv_rec_surfs[i].vertices;
        size_t type = mv_rec_surfs[i].type;

        util::matrix_t<double> nodes;
        util::matrix_t<int> elems;
        if (std::isnan(surf(0, 0))) {

            nodes.resize_fill(1,1,std::numeric_limits<double>::quiet_NaN());
            elems.resize_fill(1,1,-1);
            v_nodes.push_back(nodes);
        }
        else {
            if (type == 0) {
                // Mesh with quads
                meshMapped(surf, mv_rec_surfs[i].surf_elem_size, nodes, elems);
                v_nodes.push_back(nodes);
            }
            else if (type == 1) {
                // Mesh with triangles: meshHalfNgon method
                throw(C_csp_exception("Triangle mesh method meshHalfNgon not currently supported"));
            }
            else if (type == 2) {
                // Mesh with triangles: meshPolygon method
                throw(C_csp_exception("Triangle mesh method meshPolygon not currently supported. Need to add qhull project"));
                meshPolygon(surf, mv_rec_surfs[i].surf_elem_size);
            }
            else {
                // Mesh as a single element
                v_nodes.push_back(surf);
                elems.resize(1, surf.nrows());
                for (size_t j = 0; j < surf.nrows(); j++) {
                    elems(0,j) = j;
                }
            }

            // Shift node IDs to account for previous element sets
            add_constant_to_each_element(lastNodeID, elems);
            m_v_elems.push_back(elems);
            lastNodeID = lastNodeID + v_nodes[i].nrows();
        }
    }

    m_nodesGlobal.resize_fill(lastNodeID, 3, 0.0);
    lastNodeID = 0;
    int nodeCount = 0;
    for (size_t k = 0; k < mv_rec_surfs.size(); k++) {
        if ( std::isfinite(v_nodes[k](0,0) )) { 
            nodeCount = v_nodes[k].nrows();
            for (size_t i = 0; i < nodeCount; i++) {
                for (size_t j = 0; j < 3; j++) {
                    m_nodesGlobal(lastNodeID + i,j) = v_nodes[k](i,j);
                }
            }
            lastNodeID += nodeCount;
        }
    }

    return;
}

void C_cavity_receiver::makeGlobalElems()
{
    size_t n_surfs = mv_rec_surfs.size();

    util::matrix_t<int> type(n_surfs, 1, 0);
    util::matrix_t<int> count(n_surfs, 1, 0);

    util::matrix_t<int> int_neg(1, 1, -1);

    m_nElems = 0;

    util::matrix_t<int> temp;
    for (size_t i = 0; i < n_surfs; i++) {
        if (std::isfinite(mv_rec_surfs[i].vertices(0,0))) {

            count(i,0) = m_v_elems[i].nrows();  //[-] number of elements for surface
            type(i,0) = m_v_elems[i].ncols();   //[-] number of nodes defining the element

            // adjust to global index
            temp.resize(count(i,0),1);
            for (size_t j = 0; j < count(i, 0); j++) {
                temp(j,0) = m_nElems + j;
            }
            m_surfIDs.push_back(temp);

            m_nElems += count(i,0);     //[-] adjust global element count
        }
        else {
            m_surfIDs.push_back(int_neg);
        }
    }

    // initialize output arrays
    int typemax = round(max_int_first_column(type));
    m_elements.resize_fill(m_nElems, typemax, 0);
    m_areas.resize_fill(m_nElems, 1, std::numeric_limits<double>::quiet_NaN());
    m_centroids.resize_fill(m_nElems, 3, std::numeric_limits<double>::quiet_NaN());

    for (size_t i_surf = 0; i_surf < n_surfs; i_surf++) {
        if (std::isfinite(mv_rec_surfs[i_surf].vertices(0, 0))) {

            if (type(i_surf, 0) == 3) {
                throw(C_csp_exception("makeGlobalElems: Triangle meshes not currently supported"));
            }
            else if (type(i_surf, 0) == 4) {

                for (size_t i = 0; i < m_v_elems[i_surf].nrows(); i++) {
                    for (size_t j = 0; j < m_v_elems[i_surf].ncols(); j++) {
                        m_elements(m_surfIDs[i_surf](i,0),j) = m_v_elems[i_surf](i,j);
                    }
                }

                util::matrix_t<double> diff1, diff2, cross;
                double vecnorm;
                for (size_t i = 0; i < m_v_elems[i_surf].nrows(); i++) {
                    diffrows(m_nodesGlobal.row(m_v_elems[i_surf](i,2)), m_nodesGlobal.row(m_v_elems[i_surf](i, 0)), diff1);
                    diffrows(m_nodesGlobal.row(m_v_elems[i_surf](i,3)), m_nodesGlobal.row(m_v_elems[i_surf](i, 1)), diff2);
                    crossproduct(diff1, diff2, cross);
                    m_areas(m_surfIDs[i_surf](i,0),0) = mag_vect(cross) / 2.0;
                }

                for (size_t i = 0; i < m_v_elems[i_surf].nrows(); i++) {
                    for (size_t j = 0; j < 3; j++) {
                        m_centroids(m_surfIDs[i_surf](i, 0), j) = 0.0;
                        for (size_t k = 0; k < 4; k++) {
                            m_centroids(m_surfIDs[i_surf](i, 0), j) += m_nodesGlobal(m_v_elems[i_surf](i, k), j);
                        }
                        m_centroids(m_surfIDs[i_surf](i, 0), j) *= 0.25;
                    }
                }
            }
            else if (type(i_surf, 0) > 4) {

                for (size_t i = 0; i < m_v_elems[i_surf].nrows(); i++) {
                    for (size_t j = 0; j < m_v_elems[i_surf].ncols(); j++) {
                        m_elements(m_surfIDs[i_surf](i, 0), j) = m_v_elems[i_surf](i, j);
                    }
                }

                // Area of arbitrary polygon
                util::matrix_t<double> theseAreas(count(i_surf, 0), 1, 0.0);
                util::matrix_t<double> poly_b(type(i_surf,0), 3, std::numeric_limits<double>::quiet_NaN());
                util::matrix_t<double> poly_looped, cross, diff1, diff2, cross2, nhat;
                for (size_t i = 0; i < count(i_surf, 0); i++) {
                    // http://geomalgorithms.com/a01-_area.html
                    for (size_t j = 0; j < type(i_surf, 0); j++) {
                        for (size_t k = 0; k < 3; k++) {
                            poly_b(j,k) = m_nodesGlobal(m_v_elems[i_surf](i,j),k);
                        }
                    }

                    poly_looped = poly_b;
                    poly_looped.resize_preserve(poly_b.nrows() + 1, 3, std::numeric_limits<double>::quiet_NaN());
                    for (size_t k = 0; k < 3; k++) {
                        poly_looped(poly_looped.nrows()-1, k) = poly_b(0,k);
                    }

                    util::matrix_t<double> toSum(1,3,0.0);
                    for (size_t j = 0; j < type(i_surf, 0); j++) {
                        crossproduct(poly_looped.row(j), poly_looped.row(j+1), cross);
                        for (size_t k = 0; k < 3; k++) {
                            toSum(0,k) += cross(0,k);
                        }
                    }

                    diffrows(poly_b.row(1), poly_b.row(0), diff1);
                    diffrows(poly_b.row(3), poly_b.row(0), diff2);
                    crossproduct(diff1, diff2, cross2);
                    norm3Dvect(cross2, nhat);

                    double theseAreas = dotprod3D(nhat, toSum) / 2.0;

                    m_areas(m_surfIDs[i_surf](i, 0), 0) = theseAreas;
                }

                // Centroid of arbitrary polygon
                for (size_t i = 0; i < m_v_elems[i_surf].nrows(); i++) {
                    for (size_t j = 0; j < 3; j++) {
                        m_centroids(m_surfIDs[i_surf](i, 0), j) = 0.0;
                        for (size_t k = 0; k < type(i_surf, 0); k++) {
                            m_centroids(m_surfIDs[i_surf](i, 0), j) += m_nodesGlobal(m_v_elems[i_surf](i, k), j);
                        }
                        m_centroids(m_surfIDs[i_surf](i, 0), j) /= (double)type(i_surf, 0);
                    }
                }
            }
            else {
                throw(C_csp_exception("makeGlobalElems: Incorrect dimensions for element sets"));
            }

        }
    }

    // Calculate total active area and configure global surface index tracker
    m_global_to_surf_index.resize(m_nElems);
    m_area_active_total = 0.0;      //[m2]
    for (size_t i_surf = 0; i_surf < n_surfs; i_surf++) {
        if (std::isfinite(mv_rec_surfs[i_surf].vertices(0, 0)) && mv_rec_surfs[i_surf].is_active_surf) {
            for (size_t i = 0; i < m_v_elems[i_surf].nrows(); i++) {
                int i_global = m_surfIDs[i_surf](i, 0);
                m_area_active_total += m_areas(m_surfIDs[i_surf](i, 0), 0);
                m_global_to_surf_index[i_global] = i_surf;
            }
        }
    }

    // Create Eigen version of centroids
    matrixt_to_eigen(m_centroids, mE_centroids);

    return;
}

void C_cavity_receiver::surfValuesToElems()
{
    // CREATE EPSILON LOCAL VECTORS

    int lastElemID = 0;
    size_t n_surfs = mv_rec_surfs.size();
    std::vector<util::matrix_t<bool>> is_active_value(mv_rec_surfs.size());

    util::matrix_t<bool> temp;
    for (size_t i_surf = 0; i_surf < n_surfs; i_surf++) {

        // skip i_surf if surface not defined (e.g. no lip)
        if (std::isfinite(mv_rec_surfs[i_surf].vertices(0, 0))) {

            // elements
            int nElems = m_v_elems[i_surf].nrows();
            temp.resize_fill(nElems, 1, mv_rec_surfs[i_surf].is_active_surf);
            is_active_value[i_surf] = temp;
            lastElemID += nElems;
        }
    }

    // Combine local vectors into global vector
    m_epsilonSol.resize_fill(lastElemID, 1, std::numeric_limits<double>::quiet_NaN());
    m_epsilonTherm.resize_fill(lastElemID, 1, std::numeric_limits<double>::quiet_NaN());
    lastElemID = 0;
    for (size_t i_surf = 0; i_surf < n_surfs; i_surf++) {

        // skip i_surf if surface not defined (e.g. no lip)
        if (std::isfinite(mv_rec_surfs[i_surf].vertices(0, 0))) {
            int nElems = is_active_value[i_surf].nrows();
            for (size_t i = 0; i < nElems; i++) {
                m_epsilonSol(i+lastElemID,0) = mv_rec_surfs[i_surf].eps_sol;
                m_epsilonTherm(i + lastElemID, 0) = mv_rec_surfs[i_surf].eps_therm;
            }
            lastElemID += nElems;
        }
    }

}

void C_cavity_receiver::zigzagRouting()
{
    double tol = 0.05;      //[-] fraction of receiver height
    int maxDim = m_centroids.nrows();

    size_t maxRow = 0;
    size_t maxCol = 0;

    size_t panels_per_path = m_nPanels / m_nPaths;
    util::matrix_t<int> flow_route(m_nPaths, panels_per_path, 0);

    // resize Fluid Connectivity Array
    m_FCA.resize(m_nPaths);

    for (size_t j = 0; j < panels_per_path; j++){
        if (m_is_centerOutFlow) {
            flow_route(0,j) = panels_per_path - j -1;
            flow_route(1,j) = panels_per_path + j;
        }
        else {
            flow_route(0,j) = j;
            flow_route(1,j) = m_nPanels - j - 1;
        }
    }
        
    for (size_t h = 0; h < m_nPaths; h++) {

        util::matrix_t<int> FCM(maxDim, maxDim, -1);
        size_t count = -1;

        // For now, assume flow always in-to-out
        Eigen::VectorXi E_zagOrder;
        bool is_flow_out_to_in = false;
        if (is_flow_out_to_in) {
            E_zagOrder = Eigen::VectorXi::LinSpaced(m_pipeWindings, 0, m_pipeWindings - 1);
        }
        else {
            E_zagOrder = Eigen::VectorXi::LinSpaced(m_pipeWindings, -(m_pipeWindings - 1), 0).cwiseAbs();
        }

        // For now, assume flow always top-to-bot
        // (this implies always even number of pipeWindings)
        bool is_flow_top_to_bot = false;
        bool lastFlip = true;
        if (is_flow_top_to_bot)
            lastFlip = false;       // pipe inlet at top of receiver
        else
            lastFlip = true;        // pipe inlet at bottom of receiver

        for (size_t i_path = 0; i_path < flow_route.ncols(); i_path++) {
            size_t i = flow_route(h, i_path);

            util::matrix_t<int> elemIDs = m_surfIDs[i];
            size_t nElems = elemIDs.nrows();
            Eigen::VectorXi E_elemIDs = Eigen::VectorXi::Constant(nElems, -1);
            for (size_t j = 0; j < nElems; j++) {
                E_elemIDs(j) = elemIDs(j, 0);
            }

            util::matrix_t<double> cents(nElems, 3);
            util::matrix_t<double> cent_local;
            for (size_t j = 0; j < nElems; j++) {
                cent_local = m_centroids.row(elemIDs[j]);
                for (size_t k = 0; k < 3; k++) {
                    cents(j, k) = cent_local(0, k);
                }
            }

            Eigen::MatrixXd E_cents;
            matrixt_to_eigen(cents, E_cents);

            // translate elements to 2D domain in consistent coordinates
            Eigen::MatrixXd E_aimpoint(1, 3);
            E_aimpoint.row(0) << 10.0*0, 10.0*m_receiverHeight, 10.0*m_receiverHeight;

            Eigen::MatrixXd E_panelOrigin = furthest(E_cents, E_aimpoint);

            E_aimpoint.row(0) << 10.0, -10.0 * m_receiverHeight, 10.0 * m_receiverHeight;

            Eigen::MatrixXd E_panelXaxis = furthest(E_cents, E_aimpoint).array() - E_panelOrigin.array();

            Eigen::MatrixXd E_panelYaxis = nearest(E_cents, E_aimpoint).array() - E_panelOrigin.array();

            Eigen::Vector3d V_panelXaxis(3);
            V_panelXaxis << E_panelXaxis(0,0), E_panelXaxis(0,1), E_panelXaxis(0,2);
            Eigen::Vector3d V_panelYaxis(3);
            V_panelYaxis << E_panelYaxis(0, 0), E_panelYaxis(0, 1), E_panelYaxis(0, 2);
            Eigen::MatrixXd E_panelNorm = V_panelXaxis.cross(V_panelYaxis).transpose();

            util::matrix_t<double> panelOrigin, panelNorm, panelXaxis;
            eigen_to_matrixt(E_panelOrigin, panelOrigin);
            eigen_to_matrixt(E_panelNorm, panelNorm);
            eigen_to_matrixt(E_panelXaxis, panelXaxis);
            util::matrix_t<double> cents_2D, poly_rt;
            to2D(cents, panelOrigin, panelNorm, panelXaxis, cents_2D, poly_rt);

            Eigen::MatrixXd E_cents_2D;
            matrixt_to_eigen(cents_2D, E_cents_2D);

            // determine step heights - works best if the number of relevant
            //    elements is a multiple of nSteps
            Eigen::VectorXd A_steps = Eigen::VectorXd::LinSpaced(m_pipeWindings,
               E_cents_2D.col(0).minCoeff(), E_cents_2D.col(0).maxCoeff());

            // sort elements by the step to which they are closest
            //util::matrix_t<int> elemStep(nElems, 1, -1);
            Eigen::VectorXi E_elemStep = Eigen::VectorXi::Constant(nElems, -1);
            Eigen::VectorXd jdiff;
            double k_v_min;
            size_t k_min;
            for (size_t j = 0; j < nElems; j++) {
                jdiff = (-A_steps.array() + E_cents_2D(j,0)).abs();
                k_v_min = 1.E6;
                k_min = 0;
                for (size_t k = 0; k < jdiff.size(); k++) {
                    if (jdiff(k) < k_v_min) {
                        k_v_min = jdiff(k);
                        k_min = k;
                    }
                }
                E_elemStep(j) = k_min;
            }

            // loop through each step
            for (size_t j_index = 0; j_index < E_zagOrder.size(); j_index++) {
                int j = E_zagOrder(j_index);


                Eigen::VectorXi is_selected = (E_elemStep.array() == j).cast<int>();
                Eigen::MatrixXd E_centsStep(is_selected.sum(), E_cents_2D.cols());
                Eigen::VectorXi E_elemIDsStep(is_selected.sum());
                int rownew = 0;
                for (int ii = 0; ii < E_elemStep.rows(); ii++) {
                    if (is_selected[ii]) {
                        E_centsStep.row(rownew) = E_cents_2D.row(ii);
                        E_elemIDsStep(rownew) = E_elemIDs(ii);
                        rownew++;
                    }
                }

                // Determine number of columns and their locations
                std::vector<double> v_width(E_centsStep.rows());
                for (size_t ii = 0; ii < E_centsStep.rows(); ii++) {
                    v_width[ii] = E_centsStep(ii, 1);
                }

                std::vector<double>::iterator last = std::unique(v_width.begin(), v_width.end());
                v_width.erase(last, v_width.end());
                std::sort(v_width.begin(), v_width.end());
                last = std::unique(v_width.begin(), v_width.end());
                v_width.erase(last, v_width.end());

                int columns = v_width.size();

                if (lastFlip) {
                    lastFlip = false;
                }
                else {
                    std::reverse(v_width.begin(), v_width.end());

                    lastFlip = true;
                }

                for (size_t k = 0; k < columns; k++) {

                    Eigen::VectorXi E_test = (E_centsStep.col(1).array() <= v_width[k] + tol &&
                                        E_centsStep.col(1).array() >= v_width[k] -tol).cast<int>();

                    if (E_test.any()) {
                        count++;


                        Eigen::MatrixXi E_elemGroup(E_test.sum(),1);
                        rownew = 0;
                        for (int p = 0; p < E_test.rows(); p++) {
                            if (E_test(p)) {
                                E_elemGroup(rownew,0) = E_elemIDsStep(p);
                            }
                        }

                        E_elemGroup.transpose();

                        int elemGroupSize = E_elemGroup.rows();

                        if (elemGroupSize > maxCol) {
                            maxCol = elemGroupSize;
                        }

                        for (size_t j = 0; j < elemGroupSize; j++) {
                            FCM(count, j) = E_elemGroup(0,j);
                        }

                    }
                }
            }

            if (count > maxRow) {
                maxRow = count;
            }

        }

        m_FCA[h] = FCM;

        size_t j_col_nonzero = 0;
        size_t n_rows = m_FCA[h].nrows();

        for (size_t j = 1; j < m_FCA[h].ncols(); j++) {
            bool is_non_zero = false;
            for (size_t i = 0; i < n_rows; i++) {
                if (m_FCA[h](i, j) > -1) {
                    is_non_zero = true;
                    j_col_nonzero = j;
                    break;
                }
            }
            if (!is_non_zero) {
                break;
            }
        }

        size_t i_row_nonzero = 0;
        size_t n_cols = m_FCA[h].ncols();

        for (size_t i = 1; i < n_rows; i++) {
            bool is_non_zero = false;
            for (size_t j = 0; j < n_cols; j++) {
                if (m_FCA[h](i, j) > -1) {
                    is_non_zero = true;
                    i_row_nonzero = i;
                    break;
                }
            }
            if (!is_non_zero) {
                break;
            }
        }

        m_FCA[h].resize_preserve(i_row_nonzero + 1, j_col_nonzero + 1, -1);
    }

    // trim padding zeros
}

Eigen::MatrixXd C_cavity_receiver::furthest(const Eigen::MatrixXd cents, const Eigen::MatrixXd aimpoint)
{
    size_t nrows = cents.rows();

    double max = 0.0;
    double norm = 0.0;
    size_t i_max = 0;

    Eigen::MatrixXd diff;
    for (size_t i = 0; i < nrows; i++) {
        diff = cents.row(i).array() - aimpoint.row(0).array();
        norm = diff.norm();
        if (norm > max) {
            max = norm;
            i_max = i;
        }
    }

    return cents.row(i_max);
}

Eigen::MatrixXd C_cavity_receiver::nearest(const Eigen::MatrixXd cents, const Eigen::MatrixXd aimpoint)
{
    size_t nrows = cents.rows();

    double min = 1.E6;
    double norm = 0.0;
    size_t i_min = 0;

    Eigen::MatrixXd diff;
    for (size_t i = 0; i < nrows; i++) {
        diff = cents.row(i).array() - aimpoint.row(0).array();
        norm = diff.norm();
        if (norm < min) {
            min = norm;
            i_min = i;
        }
    }

    return cents.row(i_min);
}

void C_cavity_receiver::zigzagRouting(size_t n_steps)
{
    double tol = 0.05;      //[-] fraction of receiver height
    int maxDim = m_centroids.nrows();

    int n_surf_all = mv_rec_surfs.size();
    int n_active = 0;
    for (size_t i = 0; i < n_surf_all; i++) {
        if (mv_rec_surfs[i].is_active_surf) {
            n_active++;
        }
    }

    // resize Fluid Connectivity Array
    m_FCA.resize(n_active);

    size_t maxRow = 0;
    size_t maxCol = 0;

    for (size_t i_surf = 0; i_surf < n_active; i_surf++) {

        size_t nElems = m_surfIDs[i_surf].nrows();

        util::matrix_t<int> FCM(nElems, nElems, -1);
        util::matrix_t<double> cents(nElems, 3, std::numeric_limits<double>::quiet_NaN());

        for (size_t i = 0; i < nElems; i++) {
            for (size_t j = 0; j < 3; j++) {
                cents(i,j) = m_centroids(m_surfIDs[i_surf][i],j);
            }
        }

        int count = -1;

        // translate elements to 2D domain
        util::matrix_t<double> panelXaxis;
        diffrows(cents.row(1), cents.row(0), panelXaxis);
        util::matrix_t<double> diff1;
        diffrows(cents.row(cents.nrows()-1), cents.row(0), diff1);
        util::matrix_t<double> panelNorm;
        crossproduct(panelXaxis, diff1, panelNorm);

        util::matrix_t<double> cents2D, temp;
        to2D(cents, cents.row(0), panelNorm, panelXaxis, cents2D, temp);


        // determine step heights - this could be done better
        double min_cent = min_column_val(cents2D, 1);
        double max_cent = max_column_val(cents2D, 1);
        double delta_cent = (max_cent - min_cent) / (double)(n_steps-1);
        util::matrix_t<double> steps(n_steps, 1, std::numeric_limits<double>::quiet_NaN());
        for (size_t i = 0; i < n_steps; i++) {
            steps(i,0) = min_cent + delta_cent*i;
        }

        // Sort elements by the step to which they are closest
        util::matrix_t<int> elemStep(nElems, 1, 0);
        int i_min;
        double minval, ival;
        for (size_t j = 0; j < nElems; j++) {
            i_min = 0;
            minval = abs(cents2D(j,1) - steps(i_min,0));
            for (size_t i = 1; i < n_steps; i++) {
                ival = abs(cents2D(j, 1) - steps(i, 0));
                if (ival < minval) {
                    minval = ival;
                    i_min = i;
                }
            }
            elemStep(j,0) = i_min;
        }

        // Loop through each step
        for (int j_step = n_steps - 1; j_step >= 0; j_step--) {

            std::vector<int> elem_rows;
            for (size_t i = 0; i < nElems; i++) {
                if (elemStep(i, 0) == j_step) {
                    elem_rows.push_back(i);
                }
            }

            int n_er = elem_rows.size();
            if (n_er > 0) {
                util::matrix_t<double> centsStep(n_er, 2, std::numeric_limits<double>::quiet_NaN());
                util::matrix_t<int> elemIDsStep(n_er, 1, -1);

                for (size_t i = 0; i < n_er; i++) {
                    for (size_t j = 0; j < 2; j++) {
                        centsStep(i,j) = cents2D(elem_rows[i],j);
                        elemIDsStep(i,0) = m_surfIDs[i_surf](elem_rows[i],0);
                    }
                }

                // Determine number of columns and their locations
                std::vector<double> v_width(n_er);
                for (size_t i = 0; i < n_er; i++) {
                    v_width[i] = centsStep(i,0);
                }

                std::vector<double>::iterator last = std::unique(v_width.begin(), v_width.end());
                v_width.erase(last, v_width.end());
                std::sort(v_width.begin(), v_width.end());
                last = std::unique(v_width.begin(), v_width.end());
                v_width.erase(last, v_width.end());

                int columns = v_width.size();
                util::matrix_t<double> width(columns, 1);
                for (size_t i = 0; i < columns; i++) {
                    width(i,0) = v_width[i];
                }

                util::matrix_t<double> width_copy = width;
                if (mv_rec_surfs[i_surf].is_flipRoute) {
                    if (j_step % 2 == 1) { // alternates flow direction for each step
                        flipup(width_copy, width);
                    }
                }
                else {
                    if ((j_step + 1) % 2 == 1) {
                        flipup(width_copy, width);
                    }
                }
                
                for (size_t k = 0; k < columns; k++) {
                    std::vector<bool> test(n_er);
                    bool is_any_true = false;
                    for (size_t i = 0; i < n_er; i++) {
                        test[i] = abs(centsStep(i,0)-width(k,0)) < tol;
                        if (test[i]) {
                            is_any_true = true;
                        }
                    }

                    if (is_any_true) {
                        count++;
                    }

                    int n_test = test.size();

                    if (n_test > 0) {
                        std::vector<int> elemGroup;
                        for (size_t i = 0; i < n_test; i++) {
                            if (test[i]) {
                                elemGroup.push_back(elemIDsStep[i]);
                            }
                        }

                        if (elemGroup.size() > maxCol) {
                            maxCol = elemGroup.size();
                        }
                        
                        for (size_t j = 0; j < elemGroup.size(); j++) {
                            FCM(count,j) = elemGroup[j];
                        }
                    }
                }
            }
        }

        if (count > maxRow) {
            maxRow = count;
        }

        m_FCA[i_surf] = FCM;

        size_t j_col_nonzero = 0;
        size_t n_rows = m_FCA[i_surf].nrows();

        for (size_t j = 1; j < m_FCA[i_surf].ncols(); j++) {
            bool is_non_zero = false;
            for (size_t i = 0; i < n_rows; i++) {
                if (m_FCA[i_surf](i, j) > -1) {
                    is_non_zero = true;
                    j_col_nonzero = j;
                    break;
                }
            }
            if (!is_non_zero) {
                break;
            }
        }

        size_t i_row_nonzero = 0;
        size_t n_cols = m_FCA[i_surf].ncols();

        for (size_t i = 1; i < n_rows; i++) {
            bool is_non_zero = false;
            for (size_t j = 0; j < n_cols; j++) {
                if (m_FCA[i_surf](i, j) > -1) {
                    is_non_zero = true;
                    i_row_nonzero = i;
                    break;
                }
            }
            if (!is_non_zero) {
                break;
            }
        }

        m_FCA[i_surf].resize_preserve(i_row_nonzero + 1, j_col_nonzero + 1, -1);
    }

    return;
}

void C_cavity_receiver::VFMatrix()
{
    size_t n_surfs = m_v_elems.size();
    int nElems = 0;
    for (size_t i = 0; i < n_surfs; i++) {
        nElems += m_v_elems[i].nrows();
    }

    m_F.resize_fill(nElems, nElems, 0.0);

    int iLast = 0; // last used row index in F matrix

    util::matrix_t<double> ELEM_I;
    util::matrix_t<double> ELEM_J;

    for (size_t g_surf = 0; g_surf < n_surfs - 1; g_surf++) {   // loop though surfaces
        if (std::isfinite(mv_rec_surfs[g_surf].vertices(0, 0))) {   // check for empty
            int gElems = m_v_elems[g_surf].nrows();
            int gType = m_v_elems[g_surf].ncols();

            int jLast = iLast + gElems; // last used column index in F matrix
            for (size_t h_surf = g_surf + 1; h_surf < n_surfs; h_surf++) {
                if (std::isfinite(mv_rec_surfs[h_surf].vertices(0, 0))) {   // check for empty
                    int hElems = m_v_elems[h_surf].nrows();
                    int hType = m_v_elems[h_surf].ncols();

                    ELEM_I.resize_fill(gType, 3, std::numeric_limits<double>::quiet_NaN());
                    for (size_t i_gElems = 0; i_gElems < gElems; i_gElems++) {
                        
                        for (size_t i = 0; i < gType; i++) {
                            for (size_t j = 0; j < 3; j++) {
                                ELEM_I(i,j) = m_nodesGlobal(m_v_elems[g_surf](i_gElems,i),j);
                            }
                        }

                        ELEM_J.resize_fill(hType, 3, std::numeric_limits<double>::quiet_NaN());
                        for (size_t j_hElems = 0; j_hElems < hElems; j_hElems++) {
                            for (size_t i = 0; i < hType; i++) {
                                for (size_t j = 0; j < 3; j++) {
                                    ELEM_J(i,j) = m_nodesGlobal(m_v_elems[h_surf](j_hElems,i),j);
                                }
                            }
                            // calculate and assign view factors to F matrix
                            // use reciprocity
                            viewFactor(ELEM_I, ELEM_J, m_F(iLast+i_gElems,jLast+j_hElems), m_F(jLast+j_hElems,iLast+i_gElems));
                        }
                    }

                    jLast = jLast + hElems; // increment last used column index by the number of elements in the last set 'h'
                }
            }

            iLast = iLast + gElems; // increment last used row index by the number of elements in the last set 'g'
        }
    }

    return;
}

void C_cavity_receiver::FHatMatrix(const util::matrix_t<double>& eps,
    util::matrix_t<double>& F_hat, util::matrix_t<double>& rho,
    Eigen::MatrixXd& E_F_hat, Eigen::MatrixXd& E_rho)
{
    rho.resize_fill(m_nElems, 1, std::numeric_limits<double>::quiet_NaN());
    F_hat.resize_fill(m_nElems, m_nElems, std::numeric_limits<double>::quiet_NaN());

    E_rho.resize(m_nElems, 1);

    for (size_t i = 0; i < m_nElems; i++) {
        rho(i,0) = 1.0 - eps(i,0);
        E_rho(i,0) = rho(i,0);
    }

    util::matrix_t<double> A(m_nElems, m_nElems, 0.0);
    for (size_t i = 0; i < m_nElems; i++) {
        A(i,i) = 1.0;
        for (size_t j = 0; j < m_nElems; j++) {
            A(i,j) -= m_F(i,j)*rho(j,0);
        }
    }

    Eigen::MatrixXd Ae(m_nElems, m_nElems);
    for (size_t i = 0; i < m_nElems; i++) {
        for (size_t j = 0; j < m_nElems; j++) {
            Ae(i,j) = A(i,j);
        }
    }

    Eigen::MatrixXd Fe(m_F.nrows(), m_F.ncols());
    for (size_t i = 0; i < m_F.nrows(); i++) {
        for (size_t j = 0; j < m_F.ncols(); j++) {
            Fe(i,j) = m_F(i,j);
        }
    }

    E_F_hat = Ae.colPivHouseholderQr().solve(Fe);

    for (size_t i = 0; i < m_nElems; i++) {
        for (size_t j = 0; j < m_nElems; j++) {
            F_hat(i, j) = E_F_hat(i, j);
        }
    }

    return;
}

void C_cavity_receiver::matrixt_to_eigen(const util::matrix_t<double>& matrixt,
    Eigen::MatrixXd& eigenx)
{
    eigenx.resize(matrixt.nrows(), matrixt.ncols());
    for (size_t i = 0; i < matrixt.nrows(); i++) {
        for (size_t j = 0; j < matrixt.ncols(); j++) {
            eigenx(i,j) = matrixt(i,j);
        }
    }
}

void C_cavity_receiver::eigen_to_matrixt(const Eigen::MatrixXd& eigenx,
    util::matrix_t<double>& matrixt)
{
    size_t nrows = eigenx.rows();
    size_t ncols = eigenx.cols();
    matrixt.resize(eigenx.rows(), eigenx.cols());
    for (size_t i = 0; i < nrows; i++) {
        for (size_t j = 0; j < ncols; j++) {
            matrixt(i,j) = eigenx(i,j);
        }
    }
}

void C_cavity_receiver::hbarCorrelation(const Eigen::MatrixXd& T, double T_amb /*K*/, Eigen::MatrixXd& h /*W/m2-K*/)
{
    double A_total = mE_areas.sum() - mE_areas(mE_areas.rows()-1,0);

    double T_avg = 0.0;
    for (size_t i = 0; i < mE_areas.rows() - 1; i++) {
        T_avg += T(i,0)*(mE_areas(i,0)/A_total);
    }

    // Siebers and Kraabel
    double beta = 1.0 / T_amb;
    double nu = 1.03450643178104E-17*pow(T_amb,4) - 4.85019754418772E-14*pow(T_amb,3) + 1.3580075963433E-10*pow(T_amb,2)
                        + 2.27985665430374E-8*T_amb - 2.0313337298359E-6;
    double k = -1.24607229972985E-16*pow(T_amb,4) + 5.01096786429384E-12*pow(T_amb,3) - 2.940474355754410E-8*pow(T_amb,2)
                        + 9.05978900277077E-5*T_amb + 9.82003734668099E-4;
    double Gr = (CSP::grav * beta * (T_avg - T_amb) * pow(m_receiverHeight,3)) / pow(nu,2);
    double Nuss = 0.088*pow(Gr,(1./3.)) * pow((T_avg / T_amb),0.18);

    h.setConstant(mE_areas.rows() - 1, 1, (Nuss*k) / m_receiverHeight * 1.0);       //[W/m2-K]
}

void C_cavity_receiver::interpSolarFlux(const util::matrix_t<double>& fluxDist)
{
    for (size_t i_panel = 0; i_panel < 4; i_panel++) {

        util::matrix_t<double> panel = mv_rec_surfs[i_panel].vertices;
        util::matrix_t<double> origin = panel.row(0);
        util::matrix_t<double> xAxis, last_less_first, normal;
        diffrows(panel.row(1), origin, xAxis);
        diffrows(panel.row(panel.nrows()-1), origin, last_less_first);
        crossproduct(xAxis, last_less_first, normal);
    }
}

void C_cavity_receiver::viewFactor(const util::matrix_t<double>& poly_a, const util::matrix_t<double>& poly_b, double& F_AB, double& F_BA)
{
    double almostzero = 1.E-9;

    util::matrix_t<double> n_A;
    double areaA;
    int N;

    polygon_normal_and_area(poly_a, n_A, areaA, N);

    util::matrix_t<double> n_B;
    double areaB;
    int M;

    polygon_normal_and_area(poly_b, n_B, areaB, M);

    util::matrix_t<double> sumTerms(N, M, std::numeric_limits<double>::quiet_NaN());     // terms to sum to yield conductance
    util::matrix_t<int> skewPairs(N, M, 0);    //tracks which terms come from parallel edges

    for (size_t p = 0; p < M; p++) {    // loop through vertices of polygon B
        for (size_t i = 0; i < N; i++) {    // loop through vertices of polygon A
            util::matrix_t<double> r_i = poly_a.row(i); 
            util::matrix_t<double> r_p = poly_b.row(p);

            // loop pairings of vertices to cycle through edges
            util::matrix_t<double> r_j;
            if (i < N-1) {
                r_j = poly_a.row(i+1);
            }
            else {
                r_j = poly_a.row(0);
            }

            util::matrix_t<double> r_q;
            if (p < M-1) {
                r_q = poly_b.row(p+1);
            }
            else {
                r_q = poly_b.row(0);
            }

            // check for coincident vertices - nudge polygon B vertices if found
            //if (are_rows_equal(r_i, r_p, 0) || are_rows_equal(r_j, r_p, 0)) {
                for (size_t jj = 0; jj < r_p.ncols(); jj++) {
                    r_p(0, jj) += almostzero;
                }
            //}
            //else if (are_rows_equal(r_i, r_q, 0) || are_rows_equal(r_j, r_q, 0)) {
                for (size_t jj = 0; jj < r_q.ncols(); jj++) {
                    r_q(0, jj) += almostzero;
                }
            //}

            // determine parameterized coordinates for each edge, and minimum
            // distance between edge rays(edges extended infinitely into space)
            double dMin;
            util::matrix_t<double> sOrigin, sHat, lHat, lOrigin;
            bool skew;
            edgePairParameters(r_i, r_j, r_p, r_q, dMin, sOrigin, sHat, lHat, lOrigin, skew);

            if (skew) {
                // locate each vertex in the parameterized coordinate system
                util::matrix_t<double> ri_sO, rj_sO, rp_lO, rq_lO;
                diffrows(r_i, sOrigin, ri_sO);
                diffrows(r_j, sOrigin, rj_sO);
                diffrows(r_p, lOrigin, rp_lO);
                diffrows(r_q, lOrigin, rq_lO);

                double s_i = dotprod3D(ri_sO, sHat);
                double s_j = dotprod3D(rj_sO, sHat);
                double l_p = dotprod3D(rp_lO, lHat);
                double l_q = dotprod3D(rq_lO, lHat);

                skewPairs(i, p) = 1;
                double cosAlpha = dotprod3D(sHat, lHat);
                double alpha = acos(cosAlpha);
                double sinAlpha = sin(alpha);

                // Eq.(22a) from paper - calculate final terms that yield the
                // view factor when summed and divided by(4 * pi * area)

                sumTerms(i, p) = cosAlpha * (f_skew(s_j, l_q, alpha, cosAlpha, sinAlpha, dMin)
                    - f_skew(s_i, l_q, alpha, cosAlpha, sinAlpha, dMin)
                    - f_skew(s_j, l_p, alpha, cosAlpha, sinAlpha, dMin)
                    + f_skew(s_i, l_p, alpha, cosAlpha, sinAlpha, dMin));

                if (!isfinite(sumTerms(i, p))) {
                    throw(C_csp_exception("viewFactor skewTrue returned nan"));
                }
            }
            else {  // alternate expression for when alpha approaches zero
                lHat = sHat;
                // locate each vertex in the parameterized coordinate system
                util::matrix_t<double> ri_sO, rj_sO, rp_lO, rq_lO;
                diffrows(r_i, sOrigin, ri_sO);
                diffrows(r_j, sOrigin, rj_sO);
                diffrows(r_p, lOrigin, rp_lO);
                diffrows(r_q, lOrigin, rq_lO);

                double s_i = dotprod3D(ri_sO, sHat);
                double s_j = dotprod3D(rj_sO, sHat);
                double l_p = dotprod3D(rp_lO, lHat);
                double l_q = dotprod3D(rq_lO, lHat);

                skewPairs(i,p) = 0;
                sumTerms(i,p) = dotprod3D(sHat,lHat)*(fParallel(s_j, l_q, dMin) -
                    fParallel(s_i, l_q, dMin) - fParallel(s_j, l_p, dMin) + fParallel(s_i, l_p, dMin));

                if (!isfinite(sumTerms(i, p))) {
                    throw(C_csp_exception("viewFactor skewFalse returned nan"));
                }
            }
        }
    }

    double radUA = 0.0;
    for (size_t p = 0; p < M; p++) {    // loop through vertices of polygon B
        for (size_t i = 0; i < N; i++) {    // loop through vertices of polygon A
            radUA += sumTerms(i,p);
        }
    }

    radUA = abs(radUA)/(4.0*CSP::pi);

    F_AB = radUA / areaA;
    F_BA = radUA / areaB;

    if (!isfinite(F_AB) || !isfinite(F_BA)) {
        throw(C_csp_exception("viewFactor calculated view factor is nan"));
    }

    return;
}

double C_cavity_receiver::fParallel(double s, double l, double d)
{
    // Eq. 23 from paper
    if (d == 0.0) {
        d = 1.E-9;
    }

    double sMinusl = s - l;
    double sMinusl2 = sMinusl*sMinusl;
    double s2 = s*s;
    double l2 = l*l;
    double d2 = d*d;

    double acos_arg = min(1.0, max(-1.0, sMinusl / sqrt(s2 + l2 - 2. * s * l + d2)));
    double F = 0.5*(sMinusl2-d2)*log(sMinusl2+d2)-2.*sMinusl*d*acos(acos_arg)+s*l;
    return F;
}

double C_cavity_receiver::f_skew(double s, double l, double alpha, double cosAlpha, double sinAlpha, double d)
{
    // Eq. 22b from paper
    double s2 = s*s;
    double l2 = l*l;
    double d2 = d*d;
    double sinAlpha2 = sinAlpha*sinAlpha;

    double wsqrt = sqrt(s2 + d2/sinAlpha2);
    double psqrt = sqrt(l2 + d2/sinAlpha2);
    double wdim = 1.E-9;
    if (abs(s + wsqrt) > 0) {
        wdim = s + wsqrt;
    }
    double pdim = 1.E-9;
    if (abs(l + psqrt) > 0) {
        pdim = l + psqrt;
    }

    double F = (0.5 * cosAlpha * (s2 + l2) - s * l)* log(s2 + l2 - 2 * s * l * cosAlpha + d2)
        + s * sinAlpha * wsqrt * atan2(sqrt(s2 * sinAlpha2 + d2), (l - s * cosAlpha)) 
        + l * sinAlpha * psqrt * atan2(sqrt(l2 * sinAlpha2 + d2), (s - l * cosAlpha)) + s * l 
        + 0.5 * (d2 / sinAlpha) * (imagLi_2((wdim / pdim), alpha) + imagLi_2((pdim / wdim), alpha) - 2 * imagLi_2((wdim - 2 * s) / pdim, (CSP::pi - alpha)));

    return F;
}

double C_cavity_receiver::imagLi_2(double mag, double angle)
{
    double imaginaryPart = std::numeric_limits<double>::quiet_NaN();
    if (mag > 1.E-9) {
        double omega = atan2(mag*sin(angle), (1. - mag*cos(angle)));
        imaginaryPart = 0.5 * Cl(2 * angle) + 0.5 * Cl(2 * omega) - 0.5 * Cl(2 * omega + 2 * angle) + log(mag) * omega;
    }
    else {
        imaginaryPart = mag * sin(angle);
    }

    return imaginaryPart;
}

double C_cavity_receiver::Cl(double theta_in)
{
    double almostZero = 1.E-9;

    double theta = std::fmod(theta_in, 2.0*CSP::pi); // theta % (2.0 * CSP::pi);
    double chebArg = theta / CSP::pi - 1.0;
    double b[] = { 1.865555351433979e-1, 6.269948963579612e-2, 3.139559104552675e-4, 
        3.916780537368088e-6, 6.499672439854756e-8, 1.238143696612060e-9, 
        5.586505893753557e-13 };
    // Chebyshev polynomials of degrees 2 * n + 1 (n = 1:6) found using the sym command :
    // >> chebyshevT((2 * (0:6) + 1), sym(chebArg));
    double T[] = { chebArg, 4 * pow(chebArg,3) - 3 * chebArg, 
        16 * pow(chebArg,5) - 20 * pow(chebArg,3) + 5 * chebArg, 
        64 * pow(chebArg,7) - 112 * pow(chebArg,5) + 56 * pow(chebArg,3) - 7 * chebArg, 
        256 * pow(chebArg,9) - 576 * pow(chebArg,7) + 432 * pow(chebArg,5) - 120 * pow(chebArg,3) + 9 * chebArg,
        1024 * pow(chebArg,11) - 2816 * pow(chebArg,9) + 2816 * pow(chebArg,7) - 1232 * pow(chebArg,5) + 220 * pow(chebArg,3) - 11 * chebArg,
        4096 * pow(chebArg,13) - 13312 * pow(chebArg,11) + 16640 * pow(chebArg,9) - 9984 * pow(chebArg,7) + 2912 * pow(chebArg,5) - 364 * pow(chebArg,3) + 13 * chebArg };

    double sumbT = 0.0;
    for (size_t i = 0; i < 6; i++) {
        sumbT += b[i] * T[i];
    }

    double ClausenIntegral = (theta - CSP::pi) * (2 + log((pow(CSP::pi,2)) / 2)) + (2 * CSP::pi - theta) * log((2 * CSP::pi - theta) * (1 - almostZero) + almostZero)
        - theta * log(theta * (1 - almostZero) + almostZero) + sumbT;

    return ClausenIntegral;
}

void C_cavity_receiver::edgePairParameters(const util::matrix_t<double>& Po, const util::matrix_t<double>& Pf, const util::matrix_t<double>& Qo, const util::matrix_t<double>& Qf,
    double& D, util::matrix_t<double>& sOrigin, util::matrix_t<double>& sHat, util::matrix_t<double>& lHat, util::matrix_t<double>& lOrigin, bool& skew)
{
    // http://geomalgorithms.com/a07-_distance.html
    // find shortest distance D between line Po + s * u and Qo + t * v for initial
    // points Poand Qo, parameters sand t, and vectors uand v

    util::matrix_t<double> u, v, w;
    diffrows(Pf, Po, u);
    diffrows(Qf, Qo, v);
    diffrows(Po, Qo, w);

    util::matrix_t<double> u_copy = u;
    norm3Dvect(u_copy, u);
    util::matrix_t<double> v_copy = v;
    norm3Dvect(v_copy, v);

    double a = 1.0;
    double b = dotprod3D(u, v);
    double c = 1.0;
    double d = dotprod3D(u, w);
    double e = dotprod3D(v, w);

    double den = a*c - b*b;

    // Calculate shortest distance between edge rays
    skew = false;
    double s, l;
    if (den > 1.E-9) {
        skew = true;
        s =  (b*e - c*d)/den;
        l = (a*e - b*d)/den;
        util::matrix_t<double> vecsum(1,3,std::numeric_limits<double>::quiet_NaN());
        for (size_t j = 0; j < w.ncols(); j++) {
            vecsum(0,j) = w(0,j) + s*u(0,j) - l*v(0,j);
        }
        D = mag_vect(vecsum);
    }
    else {
        skew = false;
        s = 0.0;
        l = e / c;
        util::matrix_t<double> vecsum(1, 3, std::numeric_limits<double>::quiet_NaN());
        for (size_t j = 0; j < w.ncols(); j++) {
            vecsum(0, j) = w(0, j) - l*v(0,j);
        }
        D = mag_vect(vecsum);
    }

    // see Fig 5 in this paper:
    // Narayanaswamy, Arvind. "An analytic expression for radiation view 
    // factor between two arbitrarily oriented planar polygons." International
    // Journal of Heat and Mass Transfer 91 (2015) : 841 - 847.
    // for description of why these values are calculated in this way.
    //
    // parameter origin is location on edge ray where distance between edges has
    // its smallest value
    sOrigin.resize_fill(1, 3, std::numeric_limits<double>::quiet_NaN());
    lOrigin.resize_fill(1, 3, std::numeric_limits<double>::quiet_NaN());
    for (size_t j = 0; j < u.ncols(); j++) {
        sOrigin(0,j) = Po(0,j) + u(0,j)*s;
        lOrigin(0,j) = Qo(0,j) + v(0,j)*l;
    }

    util::matrix_t<double> pfdiff, qfdiff, podiff, qodiff;
    diffrows(Pf, sOrigin, pfdiff);
    diffrows(Qf, lOrigin, qfdiff);
    diffrows(Po, sOrigin, podiff);
    diffrows(Qo, lOrigin, qodiff);

    double s_toEnd = mag_vect(pfdiff);
    double l_toEnd = mag_vect(qfdiff);

    // unit vectors point from parameter origin to furthest of the two vertices
    if (abs(s) < s_toEnd) {
        norm3Dvect(pfdiff, sHat);
    }
    else {
        norm3Dvect(podiff, sHat);
    }
    if (abs(l) < l_toEnd) {
        norm3Dvect(qfdiff, lHat);
    }
    else {
        norm3Dvect(qodiff, lHat);
    }

}

void C_cavity_receiver::polygon_normal_and_area(const util::matrix_t<double>& poly_a,
    util::matrix_t<double>& norm_vect, double& area, int& n_rows)
{
    double almostZero = 1.E-9;

    int N = poly_a.nrows();
    n_rows = N;
    int n_verts = poly_a.ncols();

    if (n_verts != 3) {
        throw(C_csp_exception("viewFactor: requires 3 dimensions for each vertex"));
    }

    if (N < 3) {
        throw(C_csp_exception("viewFactor: need at least 3 vertices for a polygon"));
    }

    util::matrix_t<double> diff1, diff2, nHat_A, diff_local;
    diffrows(poly_a.row(1), poly_a.row(0), diff1);
    diffrows(poly_a.row(2), poly_a.row(0), diff2);
    crossproduct(diff1, diff2, norm_vect);
    norm3Dvect(norm_vect, nHat_A);

    // test for coplanar vertices
    if (N == 3) {   // test automatically satisfied
        area = mag_vect(norm_vect) / 2.0;
    }
    else {
        for (size_t i = 0; i < N; i++) {
            // the triple product of any combination of vertices must be zero
            // for the polygon to be planar
            diffrows(poly_a.row(i), poly_a.row(0), diff_local);
            double volume = abs(dotprod3D(norm_vect, diff_local));
            if (volume > almostZero) {
                throw(C_csp_exception("viewFactor: input 1 vertices not coplanar"));
            }
        }
        if (N == 4) {
            diffrows(poly_a.row(3), poly_a.row(1), diff_local);
            util::matrix_t<double> cross_local;
            crossproduct(diff2, diff_local, cross_local);
            area = mag_vect(cross_local) / 2.0;
        }
        else {
            util::matrix_t<double> poly_a_looped = poly_a;
            poly_a_looped.resize_preserve(poly_a.nrows() + 1, poly_a.ncols(), std::numeric_limits<double>::quiet_NaN());
            for (size_t j = 0; j < poly_a_looped.ncols(); j++) {
                poly_a_looped(poly_a_looped.nrows() - 1, j) = poly_a(0, j);
            }
            util::matrix_t<double> toSum(1, 3, 0.0);
            util::matrix_t<double> cross_i;
            for (size_t i = 0; i < N; i++) {
                crossproduct(poly_a_looped.row(i), poly_a_looped.row(i + 1), cross_i);
                for (size_t j = 0; j < 3; j++) {
                    toSum(0, j) += cross_i(0, j);
                }
            }
            area = dotprod3D(nHat_A, toSum) / 2.0;
        }
    }
}

void C_cavity_receiver::meshPolygon(const util::matrix_t<double>& poly, double elemSize)
{
    double almostZero = 1.E-7;

    // Confirm correct inputs
    size_t n_verts = poly.nrows();
    size_t n_dims = poly.ncols();

    util::matrix_t<double> nHat;
    util::matrix_t<double> less1_0(1, 3, std::numeric_limits<double>::quiet_NaN());
    if (n_dims == 3) {

        if (n_verts < 3) {
            throw(C_csp_exception("meshPolygon requires at least 3 vertices"));
        }

        // test for coplanar vertices
        util::matrix_t<double> poly0 = poly.row(0);
        util::matrix_t<double> poly1 = poly.row(1);
        util::matrix_t<double> poly2 = poly.row(2);

        util::matrix_t<double> less2_0(1, 3, std::numeric_limits<double>::quiet_NaN());

        diffrows(poly2, poly1, less2_0);
        diffrows(poly1, poly0, less1_0);

        util::matrix_t<double> n;
        crossproduct(less1_0, less2_0, n);

        norm3Dvect(n, nHat);

        if (n_verts > 3) { // check that points are planar

            double volume = 0.0;
            util::matrix_t<double> diff_local(1, 3, std::numeric_limits<double>::quiet_NaN());
            // Only need to check points after the first 3
            for (size_t i = 3; i < n_verts; i++) {
                // the triple product of any combination of vertices must be zero for the polygon to be planar
                // only need to check after first three points used to calculate the normal
                for (size_t j = 0; j < 3; j++) {
                    diff_local(0, j) = poly(i, j) - poly(0, j);
                }
                volume = std::abs(dotprod3D(n, diff_local));

                if (volume > almostZero) {
                    throw(C_csp_exception("meshPolygon polygon vertices not coplanar"));
                }
            }
        }
    }
    else {
        throw(C_csp_exception("meshMapped requires 3 dimensions for a vortex"));
    }

    util::matrix_t<double> center;
    ave_columns(poly, center);

    util::matrix_t<double> poly_2D;
    util::matrix_t<double> poly_rt;
    to2D(poly, center, nHat, less1_0, poly_2D, poly_rt);

    util::matrix_t<double> max_vect;
    util::matrix_t<double> min_vect;
    min_max_vects_from_columns(poly_2D, max_vect, min_vect);

    size_t n_cols = poly_2D.ncols();
    util::matrix_t<double> bbox(2, n_cols, std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < n_cols; i++) {
        bbox(0, i) = min_vect(0, i);
        bbox(1, i) = max_vect(0, i);
    }

    util::matrix_t<double> max_less_min;
    diffrows(max_vect, min_vect, max_less_min);
    double maxDim = max_row_value(max_less_min);

    if (maxDim / elemSize < 3 || maxDim / elemSize > 30) {
        throw(C_csp_exception("meshPolygon: Element size not within the required range"));
    }

    // evenly distribute mesh points on edges
    util::matrix_t<double> pfix_local = poly_2D;
    util::matrix_t<double> pfix = pfix_local;
    for (size_t j = 0; j < n_verts; j++) {

        pfix_local = pfix;

        util::matrix_t<double> pointA = poly_2D.row(j);
        util::matrix_t<double> pointB;
        if (j < n_verts - 1) {
            pointB = poly_2D.row(j+1);
        }
        else{
            pointB = poly_2D.row(0);
        }

        util::matrix_t<double> BlessA;
        diffrows(pointB, pointA, BlessA);
        double edgeLength = mag_vect(BlessA);

        // Determine number of elements on each edge
        int edgeDivs = max(1, (int)std::round(edgeLength / elemSize));

        // divide up edges into mesh points
        if (edgeDivs > 1) {
            util::matrix_t<double> AtoB;
            scale_vect(BlessA, 1./edgeLength, AtoB);
            double segment = edgeLength / (double)edgeDivs;

            util::matrix_t<double> newPoints(edgeDivs - 1, 2, 0.0);

            for (size_t i = 0; i < edgeDivs - 1; i++) {
                for (size_t k = 0; k < 2; k++) {
                    newPoints(i, k) = pointA(0, k) + (i + 1) * segment * AtoB(0, k);
                }
            }

            size_t n_row_pfix = pfix_local.nrows();
            size_t n_row_newPoints = newPoints.nrows();
            size_t n_row_pfix_new = n_row_pfix + n_row_newPoints;

            pfix.resize_preserve(n_row_pfix_new, 2, std::numeric_limits<double>::quiet_NaN());
            for (size_t i = 0; i < n_row_newPoints; i++) {
                for (size_t k = 0; k < 2; k++) {
                    pfix(n_row_pfix + i,k) = newPoints(i,k);
                }
            }
        }
        else {
            pfix = pfix_local;
        }
    }

    // Call the mesh engine
    triMesh2D(elemSize, bbox, pfix, poly_2D);

    //% call the mesh engine
    //    [nodes2D, triangles] = triMesh2D(fd, @huniform, elemSize, bbox, pfix);


}

void C_cavity_receiver::triMesh2D(double h0, const util::matrix_t<double>& bbox, const util::matrix_t<double>& pfix,
    const util::matrix_t<double>& poly_2D)
{
    // function settings
    double dptol = .001;
    double ttol = .1;
    double Fscale = 1.2;
    double deltat = .2;
    double geps = .001 * h0;
    double deps = sqrt(pow(2., -52)) * h0;


    // 1. Create initial distribution in bounding box(equilateral triangles)
    std::vector<double> x_mg;
    for (double ix = bbox(0, 0); ix < bbox(1, 0); ix = ix + h0) {
        x_mg.push_back(ix);
    }
    std::vector<double> y_mg;
    for (double iy = bbox(0, 1); iy < bbox(1, 1); iy = iy + h0 * sqrt(3) / 2.) {
        y_mg.push_back(iy);
    }

    size_t n_x_mg = x_mg.size();
    size_t n_y_mg = y_mg.size();

    util::matrix_t<double> x(n_y_mg, n_x_mg, std::numeric_limits<double>::quiet_NaN());
    util::matrix_t<double> y(n_y_mg, n_x_mg, std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < n_y_mg; i++) {
        for (size_t j = 0; j < n_x_mg; j++) {
            x(i,j) = x_mg[j];
            y(i,j) = y_mg[i];
        }
    }

    for (size_t i = 0; i < x.nrows(); i++) {
        if (i % 2 == 1) {
            for (size_t j = 0; j < x.ncols(); j++) {
                x(i, j) = x(i, j) + h0 / 2.0;       // Shift odd (even in 1-based indices) rows
            }
        }
    }

    util::matrix_t<double> p(n_y_mg*n_x_mg, 2, std::numeric_limits<double>::quiet_NaN());
    for(size_t j = 0; j < x.ncols(); j++){
        for (size_t i = 0; i < x.nrows(); i++) {
            p(j*x.nrows() + i, 0) = x(i,j);
            p(j*x.nrows() + i, 1) = y(i,j);
        }
    }

    // 2. Remove points outside the region, apply the rejection method
    util::matrix_t<double> d;
    pointToPoly(p, poly_2D, d);

    std::vector<size_t> i_p;
    for (size_t i = 0; i < d.nrows(); i++) {
        if (d(i, 0) + h0 / 2 < geps) {
            i_p.push_back(i);
        }
    }

    util::matrix_t<double> p_temp(i_p.size(), 2, std::numeric_limits<double>::quiet_NaN());
    for (size_t k = 0; k < i_p.size(); k++) {
        size_t i = i_p[k];
        p_temp(k,0) = p(i,0);
        p_temp(k,1) = p(i,1);
    }
    p = p_temp;

    std::vector<size_t> i_p_include;
    bool is_unique = true;
    for (size_t i = 0; i < p.nrows(); i++) {
        is_unique = true;
        for (size_t k = 0; k < pfix.nrows(); k++) {
            if (p(i, 0) == pfix(k, 0)) {
                if (p(i, 1) == pfix(k, 1)) {
                    is_unique = false;
                    break;
                }
            }
        }
        if (is_unique) {
            i_p_include.push_back(i);
        }
    }

    size_t n_row_pfix = pfix.nrows();
    size_t n_row_p = i_p_include.size();
    size_t n_row_p_new = n_row_pfix + n_row_p;
    p_temp = pfix;
    p_temp.resize_preserve(n_row_p_new, pfix.ncols(), std::numeric_limits<double>::quiet_NaN());

    for (size_t i = 0; i < i_p_include.size(); i++) {
        p_temp(n_row_pfix + i, 0) = p(i_p_include[i],0);
        p_temp(n_row_pfix + i, 1) = p(i_p_include[i],1);
    }

    p = p_temp;     // unordered

    vector<pair<double, double>> v_pair_p;
    for (size_t i = 0; i < p.nrows(); i++) {
        v_pair_p.push_back(make_pair(p(i,0),p(i,1)));
    }

    // using modified sort() function to sort
    sort(v_pair_p.rbegin(), v_pair_p.rend(), sort_pair_ascending);

    for (size_t i = 0; i < p.nrows(); i++) {
        p(i,0) = v_pair_p[i].first;
        p(i, 1) = v_pair_p[i].second;
    }

    util::matrix_t<double> pold(p.nrows(), p.ncols(), std::numeric_limits<double>::quiet_NaN());
    int count = 0;

    while (true) {

        // 3. Retriangulation by the Delaunay algorithm
        pold = p;       // Save current positions

    }

    /*while true
        % 3. Retriangulation by the Delaunay algorithm
        if max(sqrt(sum((p - pold). ^ 2, 2)) / h0) > ttol% Any large movement ?
            pold = p;% Save current positions
            t = delaunayn(p);% List of triangles
            pmid = (p(t(:, 1), :) + p(t(:, 2), :) + p(t(:, 3), :)) / 3; % Compute centroids
            t = t(feval(fd, pmid, varargin{ : }) < -geps, :);% Keep interior triangles
            % 4. Describe each bar by a unique pair of nodes
            bars = [t(:, [1, 2]); t(:, [1, 3]); t(:, [2, 3])];% Interior bars duplicated
            bars = unique(sort(bars, 2), 'rows');% Bars as node pairs
            % 5. Graphical output of the current mesh
            if plotFormation
                trimesh(t, p(:, 1), p(:, 2), zeros(N, 1))
                view(2); axis equal; axis off; drawnow
                end
                end*/


    return;
}

double C_cavity_receiver::pointToLine(const util::matrix_t<double>& p, const util::matrix_t<double>& a,
    const util::matrix_t<double>& b)
{
    // find the distance between point pand line segment a - b
    double x = p(0,0);
    double y = p(0,1);
    double x1 = a(0,0);
    double y1 = a(0,1);
    double x2 = b(0,0);
    double y2 = b(0,1);

    double A = x - x1;
    double B = y - y1;
    double C = x2 - x1;
    double D = y2 - y1;

    double dott = A*C + B*D;
    double len_sq = C*C + D*D;

    double param = -1.0;
    if (len_sq != 0) {
        param = dott / len_sq;
    }

    double xx = std::numeric_limits<double>::quiet_NaN();
    double yy = std::numeric_limits<double>::quiet_NaN();

    if (param < 0) {
        xx = x1;
        yy = y1;
    }
    else if (param > 1) {
        xx = x2;
        yy = y2;
    }
    else {
        xx = x1 + param*C;
        yy = y1 + param*D;
    }

    double dx = x - xx;
    double dy = y - yy;

    double val = sqrt(dx*dx + dy*dy);
    return val;
}

void C_cavity_receiver::pointToPoly(const util::matrix_t<double>& p, const util::matrix_t<double>& POLY,
    util::matrix_t<double>& d)
{
    int n = p.nrows();
    int m = p.ncols();
    int N = POLY.nrows();
    int M = POLY.ncols();

    if (m == 2 && M == 2) {

        d.resize_fill(n, 1, 0.0);

        for (size_t i = 0; i < n; i++) {

            util::matrix_t<double> D(N, 1, 0.0);

            for (size_t j = 0; j < N; j++) {
                util::matrix_t<double> a = POLY.row(j);

                util::matrix_t<double> b;
                if (j < N-1) {
                    b = POLY.row(j+1);
                }
                else {
                    b = POLY.row(0);
                }

                D(j,0) = std::abs(pointToLine(p.row(i), a, b));
            }

            d(i,0) = min_val_first_colum(D);
        }

        double abc = 1.23;
        util::matrix_t<double> p_x = p.col(0);
        util::matrix_t<double> p_y = p.col(1);

        util::matrix_t<double> poly_x;
        util::matrix_t<double> poly_y;

        transpose_matrix_t(POLY.col(0), poly_x);
        transpose_matrix_t(POLY.col(1), poly_y);

        util::matrix_t<int> in;
        inpolygon(p_x, p_y, poly_x, poly_y, in);

        for (size_t i = 0; i < n; i++) {
            d(i,0) *= (-1.0*in(i,0) + (double)(!in(i,0)));
        }
    }
    else {
        throw(C_csp_exception("pointToPoly: incorrect dimensions"));
    }
}

void C_cavity_receiver::inpolygon(const util::matrix_t<double>& p_x, const util::matrix_t<double>& p_y,
    const util::matrix_t<double>& poly_x, const util::matrix_t<double>& poly_y,
    util::matrix_t<int>& is_in_polygon)
{
    util::matrix_t<double> x = p_x;
    util::matrix_t<double> y = p_y;

    // Last point in polygon should equal first
    util::matrix_t<double> vx = poly_x;
    util::matrix_t<double> vy = poly_y;
    if (poly_x(poly_x.nrows() - 1, 0) != poly_x(0, 0) || poly_y(poly_y.nrows() - 1, 0)) {
        vx.resize_preserve(poly_x.nrows()+1, 1, std::numeric_limits<double>::quiet_NaN());
        vx(poly_x.nrows(),0) = poly_x(0,0);
        vy.resize_preserve(poly_y.nrows()+1, 1, std::numeric_limits<double>::quiet_NaN());
        vy(poly_y.nrows(),0) = poly_y(0,0);
    }

    size_t n_verts = vx.nrows();
    size_t n_points = x.ncols();

    util::matrix_t<double> xx(n_verts, n_points, std::numeric_limits<double>::quiet_NaN());
    util::matrix_t<double> yy(n_verts, n_points, std::numeric_limits<double>::quiet_NaN());
    util::matrix_t<double> v_vx(n_verts, n_points, std::numeric_limits<double>::quiet_NaN());
    util::matrix_t<double> v_vy(n_verts, n_points, std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < n_verts; i++) {
        for (size_t j = 0; j < n_points; j++) {
            xx(i,j) = x(0,j);
            yy(i,j) = y(0,j);
            v_vx(i,j) = vx(i,0);
            v_vy(i,j) = vy(i,0);
        }
    }

    x = xx;
    y = yy;
    vx = v_vx;
    vy = v_vy;

    util::matrix_t<double> avx(n_verts - 1, 1, std::numeric_limits<double>::quiet_NaN());
    util::matrix_t<double> avy(n_verts - 1, 1, std::numeric_limits<double>::quiet_NaN());
    util::matrix_t<double> ScaleFactor(n_verts - 1, 1, std::numeric_limits<double>::quiet_NaN());

    for (size_t i = 0; i < n_verts - 1; i++) {
        avx(i,0) = std::abs(0.5*(vx(i,0) + vx(i+1,0)));
        avy(i,0) = std::abs(0.5*(vy(i,0) + vy(i+1,0)));
        ScaleFactor(i,0) = max(max(avx(i,0),avy(i,0)), avx(i,0)*avy(i,0));
    }

    for (size_t i = 0; i < n_verts; i++) {
        for (size_t j = 0; j < n_points; j++) {
            vx(i,j) -= x(i,j);
            vy(i,j) -= y(i,j);
        }
    }

    util::matrix_t<int> quad(n_verts, n_points, std::numeric_limits<double>::quiet_NaN());
    bool posX, posY, negX, negY;
    for (size_t i = 0; i < n_verts; i++) {
        for (size_t j = 0; j < n_points; j++) {
            posX = vx(i, j) > 1.E-10;
            posY = vy(i, j) > 1.E-10;
            negX = !posX;
            negY = !posY;
            quad(i,j) = (size_t)(negX && posY) + 2*(size_t)(negX && negY) + 3*(size_t)(posX && negY);
        }
    }

    double scaledEps_base = sqrt(pow(2., -52)) * 3.0;
    double scaledEps = std::numeric_limits<double>::quiet_NaN();


    util::matrix_t<double> theCrossProd(n_verts-1, n_points, std::numeric_limits<double>::quiet_NaN());
    util::matrix_t<double> dotProd(n_verts - 1, n_points, std::numeric_limits<double>::quiet_NaN());
    util::matrix_t<int> signCrossProd(n_verts - 1, n_points, -99);
    util::matrix_t<int> diffQuad(n_verts - 1, n_points, -99);
    for (size_t i = 0; i < n_verts - 1; i++) {
        for (size_t j = 0; j < n_points; j++) {

            // Compute the sign() of the cross productand dot product of adjacent vertices.
            theCrossProd(i,j) = vx(i,j)*vy(i+1,j) - vx(i+1,j)*vy(i,j);

            scaledEps = ScaleFactor(i, 0) * scaledEps_base;
            if (abs(theCrossProd(i, j)) < scaledEps) {
                signCrossProd(i,j) = 0;
            }
            else if (theCrossProd(i, j) > 1.E-10) {
                signCrossProd(i,j) = 1;
            }
            else if (theCrossProd(i, j) < -1.E-10) {
                signCrossProd(i, j) = -1;
            }
            else {
                signCrossProd(i, j) = 0;
            }

            dotProd(i,j) = vx(i,j)*vx(i+1,j) + vy(i,j)*vy(i+1,j);

            diffQuad(i,j) = quad(i+1,j) - quad(i,j);

            if (abs(diffQuad(i, j)) == 3) {
                diffQuad(i,j) /= -3;
            }
            if (abs(diffQuad(i, j)) == 2) {
                diffQuad(i,j) = 2 * signCrossProd(i,j);
            }
        }        
    }

    util::matrix_t<int> in(1, n_points, -1);

    sum_int_columns(diffQuad, in);

    util::matrix_t<int> on(1, n_points, -1);

    size_t is_true = 0;
    for (size_t j = 0; j < n_points; j++) {

        if (in(0, j) != 0) {
            in(0, j) = 1;
        }
        else {
            in(0, j) = 0;
        }

        is_true = 0;
        for (size_t i = 0; i < n_verts - 1; i++) {
            if (signCrossProd(i, j) == 0 && dotProd(i, j) <= 0.0) {
                is_true = 1;
            }
        }
        on(0,j) = is_true;
        if (in(0,j) == 1 || on(0,j) == 1) {
            in(0,j) = 1;
        }
        else {
            in(0,j) = 0;
        }
    }

    transpose_int_matrix_t(in, is_in_polygon);

    return;
}


void C_cavity_receiver::transpose_matrix_t(const util::matrix_t<double>& a, util::matrix_t<double>& b)
{
    size_t n_row = a.nrows();
    size_t n_col = a.ncols();
    b.resize_fill(n_col, n_row, std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < n_row; i++) {
        for (size_t j = 0; j < n_col; j++) {
            b(j,i) = a(i,j);
        }
    }
}

void C_cavity_receiver::transpose_int_matrix_t(const util::matrix_t<int>& a, util::matrix_t<int>& b)
{
    size_t n_row = a.nrows();
    size_t n_col = a.ncols();
    b.resize_fill(n_col, n_row, std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < n_row; i++) {
        for (size_t j = 0; j < n_col; j++) {
            b(j, i) = a(i, j);
        }
    }
}


void C_cavity_receiver::meshMapped(const util::matrix_t<double>& poly, double elemSize,
    util::matrix_t<double>& nodes, util::matrix_t<int>& quads)
{
    double almostZero = 1.E-7;

    // Confirm correct inputs
    size_t n_verts = poly.nrows();
    size_t n_dims = poly.ncols();

    util::matrix_t<double> less1_0(1, 3, std::numeric_limits<double>::quiet_NaN());
    util::matrix_t<double> n_hat;
    if (n_dims == 3) {

        // Test for coplanar vertices
        if (n_verts != 4) {
            throw(C_csp_exception("meshMapped requires 4 vertices"));
        }
        else {
            util::matrix_t<double> less2_0(1, 3, std::numeric_limits<double>::quiet_NaN());
            for (size_t i = 0; i < 3; i++) {
                less1_0(0, i) = poly(1, i) - poly(0, i);
                less2_0(0, i) = poly(2, i) - poly(0, i);
            }
            util::matrix_t<double> n;
            crossproduct(less1_0, less2_0, n);
            norm3Dvect(n, n_hat);
            double volume = 0.0;
            util::matrix_t<double> diff_local(1, 3, std::numeric_limits<double>::quiet_NaN());
            for (size_t i = 3; i < n_verts; i++) {
                // the triple product of any combination of vertices must be zero for the polygon to be planar
                // only need to check after first three points used to calculate the normal
                for (size_t j = 0; j < 3; j++) {
                    diff_local(0, j) = poly(i, j) - poly(0, j);
                }
                volume = std::abs(dotprod3D(n, diff_local));

                if (volume > almostZero) {
                    throw(C_csp_exception("meshMapped polygon vertices not coplanar"));
                }
            }
        }
    }
    else {
        throw(C_csp_exception("meshMapped requires 3 dimensions for a vortex"));
    }

    util::matrix_t<double> center;
    ave_columns(poly, center);

    util::matrix_t<double> poly_2D;
    util::matrix_t<double> poly_rt;
    to2D(poly, center, n_hat, less1_0, poly_2D, poly_rt);

    util::matrix_t<double> nodes2D;
    map(poly_2D, elemSize, nodes2D, quads);

    to3D(nodes2D, center, n_hat, less1_0, nodes);

    return;
}

void C_cavity_receiver::to3D(const util::matrix_t<double>& poly_xy, const util::matrix_t<double>& origin,
    const util::matrix_t<double>& normal, const util::matrix_t<double>& xaxis,
    util::matrix_t<double>& poly3d)
{
    size_t n = poly_xy.nrows();     // number of points to process

    util::matrix_t<double> nHat;
    norm3Dvect(normal, nHat);

    util::matrix_t<double> xHat;
    norm3Dvect(xaxis, xHat);

    util::matrix_t<double> yHat;
    crossproduct(nHat, xHat, yHat);

    poly3d.resize_fill(n, 3, 0.0);
    for (size_t i = 0; i < n; i++){
        for (size_t j = 0; j < 3; j++) {
            poly3d(i,j) = origin(0,j) + xHat(0,j)*poly_xy(i,0) + yHat(0,j)*poly_xy(i,1);
        }
    }

    return;
}

void C_cavity_receiver::map(const util::matrix_t<double>& poly2D, double elemSize,
    util::matrix_t<double>& nodes, util::matrix_t<int>& quads)
{
    util::matrix_t<double> A = poly2D.row(0);
    util::matrix_t<double> B = poly2D.row(1);
    util::matrix_t<double> C = poly2D.row(2);
    util::matrix_t<double> D = poly2D.row(3);

    util::matrix_t<double> AtoB;
    diffrows(B, A, AtoB);
    double lengthAB = mag_vect(AtoB);

    util::matrix_t<double> BtoC;
    diffrows(C, B, BtoC);
    double lengthBC = mag_vect(BtoC);

    util::matrix_t<double> CtoD;
    diffrows(D, C, CtoD);
    double lengthCD = mag_vect(CtoD);

    util::matrix_t<double> AtoD;
    diffrows(D, A, AtoD);
    double lengthDA = mag_vect(AtoD);

    double maxDim = max({ lengthAB, lengthBC, lengthCD, lengthDA });

    // test for a reasonable element size
    if (maxDim / elemSize < 0.5 || maxDim / elemSize > 250) {
        throw(C_csp_exception("Element size not within the required range"));
    }

    // find reasonable # of elements per edge for mapping
    size_t elemsM = std::max(1, (int) std::round(0.5/elemSize*(lengthAB + lengthCD)));
    size_t elemsN = std::max(1, (int) std::round(0.5/elemSize*(lengthBC + lengthDA)));

    // initialize nodeand element arrays
    nodes.resize_fill((elemsM+1)*(elemsN+1), 2, 0.0);
    quads.resize_fill(elemsM*elemsN, 4, -1);

    int nodeID = -1;     // last used node index - set to -1 for C++ 0-based index
    int elemID = -1;     // last used element index - set to -1 for C++ 0-based index

    // Matlab code also starts at index 0
    for (size_t m = 0; m < elemsM + 1; m++) {
        // create bridge EF between AB and CD
        util::matrix_t<double> scaledAtoB;
        scale_vect(AtoB, m/(double)elemsM, scaledAtoB);
        util::matrix_t<double> E;
        add_vect_rows(A, scaledAtoB, E);

        util::matrix_t<double> scaledCtoD;
        scale_vect(CtoD, m/(double)elemsM, scaledCtoD);
        util::matrix_t<double> F;
        diffrows(D, scaledCtoD, F);

        util::matrix_t<double> EtoF;
        diffrows(F, E, EtoF);

        // Matlab code also starts at index 0
        for (size_t n = 0; n < elemsN + 1; n++) {
            // walk along bridge EF defining nodes
            nodeID++;
            util::matrix_t<double> scaledEtoF;
            scale_vect(EtoF, n/(double)elemsN, scaledEtoF);
            util::matrix_t<double> E_plus_scaledEtoF;
            add_vect_rows(E, scaledEtoF, E_plus_scaledEtoF);

            for (size_t i = 0; i < E_plus_scaledEtoF.ncols(); i++) {
                nodes(nodeID, i) = E_plus_scaledEtoF(0, i);
            }

            if (m > 0 && n > 0) { // then define new element
                elemID++;

                // Retrieve nodes around element in CCW direction
                // these MUST go around the element in order(CW or CCW) in
                // order to work correctly with viewFactor(...)
                quads(elemID, 0) = nodeID - elemsN - 2;
                quads(elemID, 1) = nodeID - 1;
                quads(elemID, 2) = nodeID;
                quads(elemID, 3) = nodeID - elemsN - 1;
            }

        }
    }

    return;
}

void C_cavity_receiver::to2D(const util::matrix_t<double>& poly, const util::matrix_t<double>& center,
    const util::matrix_t<double>& normal, const util::matrix_t<double>& xaxis,
    util::matrix_t<double>& poly_xy, util::matrix_t<double>& poly_rt)
{
    size_t n = poly.nrows();

    util::matrix_t<double> nHat;
    norm3Dvect(normal, nHat);

    util::matrix_t<double> xHat;
    norm3Dvect(xaxis, xHat);

    util::matrix_t<double> yHat;
    crossproduct(nHat, xHat, yHat);

    poly_xy.resize_fill(n, 2, 0.0);
    poly_rt.resize_fill(n, 2, 0.0);

    for (size_t i = 0; i < n; i++) {
        util::matrix_t<double> point = poly.row(i);
        util::matrix_t<double> arm(1, 3);
        for (size_t j = 0; j < 3; j++) {
            arm(0,j) = point(0,j) - center(0,j);
        }
        double radius = mag_vect(arm);
        double xComp = dotprod3D(arm, xHat);    // x coordinate in 2D CS
        double yComp = dotprod3D(arm, yHat);    // y coordinate in 2D CS
        double theta = atan2(yComp, xComp);     //
        if (theta < 0.0) {
            theta = theta + 2. * CSP::pi;
        }
        poly_xy(i, 0) = xComp;
        poly_xy(i, 1) = yComp;
        poly_rt(i, 0) = radius;
        poly_rt(i, 1) = theta;
    }

    return;
}

void C_cavity_receiver::add_constant_to_each_element(int val, util::matrix_t<int>& a)
{
    for (size_t i = 0; i < a.nrows(); i++) {
        for (size_t j = 0; j < a.ncols(); j++) {
            a(i,j) += val;
        }
    }
}

double C_cavity_receiver::dotprod3D(const util::matrix_t<double>& a, const util::matrix_t<double>& b)
{
    return a(0,0)*b(0,0) + a(0,1)*b(0,1) + a(0,2)*b(0,2);
}

void C_cavity_receiver::scale_vect(const util::matrix_t<double>& a, double scale, util::matrix_t<double>& out_vect)
{
    out_vect = a;
    for (size_t i = 0; i < a.ncols(); i++) {
        out_vect(0,i) = a(0,i)*scale;
    }
}

void C_cavity_receiver::crossproduct(const util::matrix_t<double>& a_vert, const util::matrix_t<double>& b_vert, util::matrix_t<double>& cross)
{
    cross.resize_fill(1, 3, std::numeric_limits<double>::quiet_NaN());

    cross(0,0) = a_vert(0,1)*b_vert(0,2) - a_vert(0,2)*b_vert(0,1);
    cross(0,1) = a_vert(0,2)*b_vert(0,0) - a_vert(0,0)*b_vert(0,2);
    cross(0,2) = a_vert(0,0)*b_vert(0,1) - a_vert(0,1)*b_vert(0,0);
}

void C_cavity_receiver::norm3Dvect(const util::matrix_t<double>& vector_in, util::matrix_t<double>& norm_vect)
{
    norm_vect.resize_fill(1, 3, std::numeric_limits<double>::quiet_NaN());
    double magnitude = mag_vect(vector_in);
    for (size_t i = 0; i < 3; i++) {
        norm_vect(0, i) = vector_in(0, i) / magnitude;
    }
}

double C_cavity_receiver::mag_vect(const util::matrix_t<double>& vector_in)
{
    double sum_of_sq = 0.0;
    for (size_t i = 0; i < vector_in.ncols(); i++) {
        sum_of_sq += std::pow(vector_in(0, i), 2);
    }

    return std::sqrt(sum_of_sq);

    //return std::sqrt(std::pow(vector_in(0, 0), 2) + std::pow(vector_in(0, 1), 2) + std::pow(vector_in(0, 2), 2));
}

void C_cavity_receiver::flipup(const util::matrix_t<double>& a, util::matrix_t<double>& b)
{
    size_t nrows = a.nrows();
    size_t ncols = a.ncols();
    b.resize_fill(nrows, ncols, std::numeric_limits<double>::quiet_NaN());

    for (size_t i = 0; i < nrows; i++) {
        for (size_t j = 0; j < ncols; j++) {
            b(i,j) = a(nrows-1-i,j);
        }
    }
}

void C_cavity_receiver::sumcolumns(const util::matrix_t<double>& a, util::matrix_t<double>& summed)
{
    size_t ncols = a.ncols();
    summed.resize_fill(1, ncols, 0.0);

    for (size_t i = 0; i < a.nrows(); i++) {
        for (size_t j = 0; j < ncols; j++) {
            summed(0, j) += a(i, j);
        }
    }
}

void C_cavity_receiver::sum_int_columns(const util::matrix_t<int>& a, util::matrix_t<int>& summed)
{
    size_t ncols = a.ncols();
    summed.resize_fill(1, ncols, 0.0);

    for (size_t i = 0; i < a.nrows(); i++) {
        for (size_t j = 0; j < ncols; j++) {
            summed(0, j) += a(i, j);
        }
    }
}

void C_cavity_receiver::diffrows(const util::matrix_t<double>& a, const util::matrix_t<double>& b, util::matrix_t<double>& a_less_b)
{
    a_less_b.resize_fill(1, a.ncols(), std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < a.ncols(); i++) {
        a_less_b(0, i) = a(0, i) - b(0, i);
    }
}

void C_cavity_receiver::add_vect_rows(const util::matrix_t<double>& a, const util::matrix_t<double>& b, util::matrix_t<double>& a_plus_b)
{
    a_plus_b.resize_fill(1, a.ncols(), std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < a.ncols(); i++) {
        a_plus_b(0,i) = a(0,i) + b(0,i);
    }
}

void C_cavity_receiver::ave_columns(const util::matrix_t<double>& a, util::matrix_t<double>& averaged)
{
    double nrows = (double)a.nrows();
    sumcolumns(a, averaged);
    for (size_t i = 0; i < 3; i++) {
        averaged(0, i) /= nrows;
    }
}

double C_cavity_receiver::max_row_value(const util::matrix_t<double>& a)
{
    double maxval = a(0, 0);
    for (size_t i = 1; i < a.ncols(); i++) {
        maxval = max(maxval, a(0,i));
    }
    return maxval;
}

int C_cavity_receiver::max_row_int_value(const util::matrix_t<int>& a)
{
    int maxval = a(0, 0);
    for (size_t i = 1; i < a.ncols(); i++) {
        maxval = max(maxval, a(0, i));
    }
    return maxval;
}

double C_cavity_receiver::max_column_val(const util::matrix_t<double>& a, size_t n_c)
{
    double maxval = a(0, n_c);
    for (size_t i = 1; i < a.nrows(); i++) {
        maxval = max(maxval, a(i, n_c));
    }
    return maxval;
}

double C_cavity_receiver::min_val_first_colum(const util::matrix_t<double>& a)
{
    double minval = a(0, 0);
    for (size_t i = 1; i < a.nrows(); i++) {
        minval = min(minval, a(i, 0));
    }
    return minval;
}

double C_cavity_receiver::min_column_val(const util::matrix_t<double>& a, size_t n_c)
{
    double minval = a(0, n_c);
    for (size_t i = 1; i < a.nrows(); i++) {
        minval = min(minval, a(i, n_c));
    }
    return minval;
}

int C_cavity_receiver::max_int_first_column(const util::matrix_t<int>& a)
{
    int maxval = a(0, 0);
    for (size_t i = 1; i < a.nrows(); i++) {
        maxval = max(maxval, a(i, 0));
    }
    return maxval;
}

bool C_cavity_receiver::are_rows_equal(const util::matrix_t<double>& a, const util::matrix_t<double>& b, int i_row)
{
    size_t n_col_a = a.ncols();
    size_t n_col_b = b.ncols();
    if (n_col_a != n_col_b) {
        return false;
    }
    for (size_t j = 0; j < n_col_a; j++) {
        if (a(i_row, j) != b(i_row, j)) {
            return false;
        }
    }

    return true;
}

void C_cavity_receiver::min_max_vects_from_columns(const util::matrix_t<double>& a, util::matrix_t<double>& max_vect, util::matrix_t<double>& min_vect)
{
    size_t ncols = a.ncols();
    max_vect = a.row(0);
    min_vect = a.row(0);
    for (size_t i = 1; i < a.nrows(); i++) {
        for (size_t j = 0; j < ncols; j++) {
            max_vect(0,j) = max(max_vect(0,j), a(i,j));
            min_vect(0,j) = min(min_vect(0,j), a(i,j));
        }
    }

}

void C_cavity_receiver::tube_UA_and_deltaP(std::vector<double> m_dot_paths /*kg/s*/, const Eigen::MatrixXd E_T_HTF /*K*/,
    Eigen::MatrixXd & UA, double & W_dot_pump /*MWe*/)
{
    UA.resize(m_nElems, 1);
    UA.setConstant(0.0);

    double W_dot_htf_dP = 0.0;    //[We]
    // (using "path" here instead of "pipe" from Matlab code)
    for (size_t i_path = 0; i_path < m_nPaths; i_path++) {

        double m_dot_tube = m_dot_paths[i_path] / (double)m_Ntubes; //[kg/s]

        util::matrix_t<int> FCM = m_FCA[i_path];

        size_t nSteps = FCM.nrows();

        if (FCM.ncols() > 1) {
            throw(C_csp_exception("Cavity code currently does not support parallel nodes within a fluid connectivity matrix",
                "C_cavity_receiver"));
        }

        //double l_last = 0;
        double lTotal = 0;

        double dP_path = 0.0;

        double T_in = E_T_HTF(FCM(0, 0),0);   //[K]
        double T_out = E_T_HTF(FCM(nSteps-1,0),0);    //[K]
        double rho_avg = field_htfProps.dens(0.5*(T_in + T_out),1.0);

        for (size_t i = 0; i < nSteps; i++) {

            size_t stepID = FCM(i,0);

            double T_htf = E_T_HTF(stepID,0);           //[K]

            // Use HTF temperature to calculate tube conductivity
            double ktube = tube_material.cond(T_htf); //[W/m-K]

            double k = field_htfProps.cond(T_htf);          //[W/m-K]
            double rho = field_htfProps.dens(T_htf, 1.0);   //[kg/m3]
            double mu = field_htfProps.visc(T_htf);         //[kg/m-s]
            double cp = field_htfProps.Cp(T_htf)*1.E3;      //[J/kg-K] convert from kJ/kg-K
            double u = m_dot_tube/(rho*m_A_cs_tube);        //[m/s]
            double Re = rho*u*m_d_in_rec_tube/mu;           //[-]
            double Pr = cp*mu/k;                            //[-]

            int i_panel = m_global_to_surf_index[stepID];
            double i_elem_size = mv_rec_surfs[i_panel].surf_elem_size;  //[m]
            double l = mE_areas(stepID,0) / i_elem_size;       //[m]

            lTotal += l;

            // Determine convection coefficient
            double Nu, f;
            CSP::PipeFlow(Re, Pr, lTotal / m_d_in_rec_tube, m_rel_roughness, Nu, f);

            double h = Nu * k / m_d_in_rec_tube;    //[W/m2-K]
            double Rcond = log(m_od_rec_tube/m_d_in_rec_tube)/(CSP::pi*l*ktube*m_Ntubes);
            double Rconv = 2.0/(h*m_Ntubes*l*m_d_in_rec_tube*CSP::pi);
            UA(stepID,0) = 1.0/(Rcond + Rconv);     //[W/K]

            if (i == nSteps - 1) {
                dP_path = rho*f*lTotal*u*u/(2.0*m_d_in_rec_tube);   //[Pa]
            }
        }

        W_dot_htf_dP += dP_path * m_dot_paths[i_path] / rho_avg; //[W]
    }

    double m_dot_htf_tot = std::accumulate(m_dot_paths.begin(), m_dot_paths.end(), 0.0);   //[kg/s]

    double est_load = fmax(0.25, m_dot_htf_tot / m_m_dot_htf_des) * 100;		//[%] Relative pump load. Limit to 25%
    double eta_pump_adj = m_eta_pump * (-2.8825E-9 * pow(est_load, 4) + 6.0231E-7 * pow(est_load, 3) - 1.3867E-4 * pow(est_load, 2) + 2.0683E-2 * est_load);	//[-] Adjusted pump efficiency

    W_dot_pump = W_dot_htf_dP / eta_pump_adj * 1.E-6;   //[MWe] convert from W and derate for pump efficiency

    return;
}

void C_cavity_receiver::init()
{
    // ******************************************
    // Unit conversions
    // ******************************************
    m_T_htf_hot_des += 273.15;	//[K] Convert from input in [C]
    m_T_htf_cold_des += 273.15;	//[K] Convert from input in [C]
    m_q_rec_des *= 1.E6;		//[W] Convert from input in [MW]

    // ******************************************
    // Set up cavity geometry and view factors
    // ******************************************
    m_nPanels = 4;
    
    m_pipeWindings = 4;   //[-] Probably needs to be >= 2 to avoid inconsistencies in mesh calcs
    m_modelRes = 1;               //[-] Value must be = 1 until/unless modelRes code imported from Matlab
    m_is_bottomUpFlow = true;
    m_is_centerOutFlow = true;
    m_nPaths = 2;

    double surface_roughness = 4.5e-5;  //[m]

    // Create geometry(i.e.defines vertices) for a 4 - panel half - octagonal cavity receiver
    genOctCavity();

    // Meshes each surface
    meshGeometry();

    // Tube geometry calcs
    m_d_in_rec_tube = m_od_rec_tube - 2.0 * m_th_rec_tube;      //[m]
    m_A_cs_tube = 0.25*CSP::pi*pow(m_d_in_rec_tube,2);          //[m2]
    m_rel_roughness = surface_roughness / m_d_in_rec_tube;      //[-]

    // choose number of HTF tubes per route based on available space
    //    assumes each receiver panel is same area and uses same dimension tube
    m_Ntubes = std::floor(mv_rec_surfs[0].surf_elem_size*m_modelRes/m_od_rec_tube);

    // Make global elements and calculate element centroids and areas
    makeGlobalElems();

    // Assign global element solar and thermal emissivity
    surfValuesToElems();

    // Define fluid connectivity array
    zigzagRouting();

    // Calculate view factors
    VFMatrix();

    // Calculate FHat matrices
    util::matrix_t<double> rhoSol;
    FHatMatrix(m_epsilonSol, m_FHatS, rhoSol, mE_FHatS, mE_rhoSol);

    util::matrix_t<double> rhoTherm;
    FHatMatrix(m_epsilonTherm, m_FHatT, rhoTherm, mE_FHatT, mE_rhoTherm);

    matrixt_to_eigen(m_epsilonSol, mE_epsilonSol);
    matrixt_to_eigen(m_epsilonTherm, mE_epsilonTherm);
    matrixt_to_eigen(m_areas, mE_areas);

    // ********************************************
    // ********************************************

    // ********************************************
    // Complete receiver initialization

    ambient_air.SetFluid(ambient_air.Air);

    // Declare instance of fluid class for FIELD fluid
    if (m_field_fl != HTFProperties::User_defined && m_field_fl < HTFProperties::End_Library_Fluids)
    {
        if (!field_htfProps.SetFluid(m_field_fl))
        {
            throw(C_csp_exception("Receiver HTF code is not recognized", "MSPT receiver"));
        }
    }
    else if (m_field_fl == HTFProperties::User_defined)
    {
        // Check that 'm_field_fl_props' is allocated and correct dimensions
        int n_rows = (int)m_field_fl_props.nrows();
        int n_cols = (int)m_field_fl_props.ncols();
        if (n_rows > 2 && n_cols == 7)
        {
            if (!field_htfProps.SetUserDefinedFluid(m_field_fl_props))
            {
                error_msg = util::format(field_htfProps.UserFluidErrMessage(), n_rows, n_cols);
                throw(C_csp_exception(error_msg, "MSPT receiver"));
            }
        }
        else
        {
            error_msg = util::format("The user defined field HTF table must contain at least 3 rows and exactly 7 columns. The current table contains %d row(s) and %d column(s)", n_rows, n_cols);
            throw(C_csp_exception(error_msg, "MSPT receiver"));
        }
    }
    else
    {
        throw(C_csp_exception("Receiver HTF code is not recognized", "MSPT receiver"));
    }

    // Declare instance of htf class for receiver tube material
    if (m_tube_mat_code == HTFProperties::Stainless_AISI316 || m_tube_mat_code == HTFProperties::T91_Steel ||
        m_tube_mat_code == HTFProperties::N06230 || m_tube_mat_code == HTFProperties::N07740)
    {
        if (!tube_material.SetFluid(m_tube_mat_code))
        {
            throw(C_csp_exception("Tube material code not recognized", "MSPT receiver"));
        }
    }
    else if (m_tube_mat_code == HTFProperties::User_defined)
    {
        throw(C_csp_exception("Receiver material currently does not accept user defined properties", "MSPT receiver"));
    }
    else
    {
        error_msg = util::format("Receiver material code, %d, is not recognized", m_tube_mat_code);
        throw(C_csp_exception(error_msg, "MSPT receiver"));
    }

    // Calculate constant thermal piping losses to the environment
    if (m_pipe_loss_per_m > 0.0 && m_pipe_length_mult > 0.0)
        m_Q_dot_piping_loss = m_pipe_loss_per_m * (m_h_tower * m_pipe_length_mult + m_pipe_length_add);		//[Wt]
    else
        m_Q_dot_piping_loss = 0.0;

    m_mode_prev = C_csp_collector_receiver::OFF;
    m_od_control = 1.0;			                //[-] Additional defocusing for over-design conditions
    m_E_su_prev = m_q_rec_des * m_rec_qf_delay;	//[W-hr] Startup energy
    m_t_su_prev = m_rec_su_delay;				//[hr] Startup time requirement

	return;
}

void C_cavity_receiver::call(const C_csp_weatherreader::S_outputs& weather,
	const C_csp_solver_htf_1state& htf_state_in,
	const C_pt_receiver::S_inputs& inputs,
	const C_csp_solver_sim_info& sim_info)
{
    // Get inputs
    double field_eff = inputs.m_field_eff;					//[-]
    const util::matrix_t<double>* flux_map_input = inputs.m_flux_map_input;
    C_csp_collector_receiver::E_csp_cr_modes input_operation_mode = inputs.m_input_operation_mode;

    // Get sim info 
    double step = sim_info.ms_ts.m_step;			//[s]
    double time = sim_info.ms_ts.m_time;	//[s]

    // Get applicable htf state info
    double T_salt_cold_in = htf_state_in.m_temp + 273.15;		//[K] convert from C
    double T_amb = weather.m_tdry + 273.15;                     //[K] convert from C

    // Read in remaining weather inputs from weather output structure
    double zenith = weather.m_solzen;       //[deg]
    double azimuth = weather.m_solazi;      //[deg]
    double v_wind_10 = weather.m_wspd;      //[m/s]
    double I_bn = weather.m_beam;           //[W/m2]

    // Solution parameters
    
    double tol_T_HTF_node_iter = 0.1;       //[K]
    double error_T_HTF_node_iter = std::numeric_limits<double>::quiet_NaN();    //[K]
    size_t count_T_HTF_node_iter = 0;       //[-]

    double tol_rel_T_rec_node_iter = 1.E-4; //[-]
    double error_relmax_T_rec_node_iter = std::numeric_limits<double>::quiet_NaN(); //[-]
    size_t count_T_rec_node_iter = 0;       //[-]

    double tol_abs_T_htf_target = 0.1;      //[K]
    double error_T_htf_out = std::numeric_limits<double>::quiet_NaN();  //[K]
    size_t count_T_htf_out_iter = 0;        //[-]

    //[-] Guess design field efficiency
    // Used to scale hard-coded flux map
    double field_eff_des = 0.6;     

    bool debugthis = false;
    if (debugthis) {
        zenith = 0.0;
        azimuth = 0.0;
        I_bn = 950.0;
        input_operation_mode = C_csp_collector_receiver::E_csp_cr_modes::ON;
        T_salt_cold_in = 563.15;
        m_T_htf_hot_des = 848.15;
        T_amb = 20 + 273.15;        //[K] Temperature of surroundings
        field_eff = field_eff_des;  //[-] set efficiency to design to scale flux map same as Matlab code
        m_Q_dot_piping_loss = 0.0;  //[-] Matlab code doesn't consider piping losses
    }



    bool rec_is_off = false;
    bool rec_is_defocusing = false;
    double field_eff_adj = 0.0;

    // Do an initial check to make sure the solar position called is valid
    // If it's not, return the output equal to zeros. Also check to make sure
    // the solar flux is at a certain level, otherwise the correlations aren't valid
    if (input_operation_mode == C_csp_collector_receiver::OFF)
    {
        rec_is_off = true;
    }

    if (zenith > (90.0 - m_hel_stow_deploy) || I_bn <= m_f_rec_min*m_dni_des || (zenith == 0.0 && azimuth == 180.0))
    {
        m_mode = C_csp_collector_receiver::OFF;
        rec_is_off = true;
    }

    double T_coolant_prop = (m_T_htf_hot_des + T_salt_cold_in) / 2.0;		//[K] The temperature at which the coolant properties are evaluated. Validated as constant (mjw)
    double cp_htf = field_htfProps.Cp(T_coolant_prop) * 1000.0;				//[J/kg-K] Specific heat of the coolant

    double W_dot_pump, DELTAP, Pres_D, u_coolant;
    W_dot_pump = DELTAP = Pres_D = u_coolant = std::numeric_limits<double>::quiet_NaN();

    if (field_eff < m_eta_field_iter_prev && m_od_control < 1.0)
    {	// Suggests controller applied defocus, so reset *controller* defocus
        m_od_control = fmin(m_od_control + (1.0 - field_eff / m_eta_field_iter_prev), 1.0);
    }

    double m_dot_htf_tot = std::numeric_limits<double>::quiet_NaN();        //[kg/s]
    double T_htf_tower_out_calc = std::numeric_limits<double>::quiet_NaN();       //[K]
    double T_htf_rec_no_tower_losses = std::numeric_limits<double>::quiet_NaN();  //[K]
    double q_dot_inc = std::numeric_limits<double>::quiet_NaN();            //[Wt]
    double q_dot_thermal_tower_out = std::numeric_limits<double>::quiet_NaN();        //[Wt]
    double q_dot_thermal_rec_out_no_tower_losses = std::numeric_limits<double>::quiet_NaN();    //[Wt]
    double eta_thermal_calc = std::numeric_limits<double>::quiet_NaN();     //[-]
    double q_gain_net_eb = std::numeric_limits<double>::quiet_NaN();        //[Wt]
    double q_dot_conv_losses = std::numeric_limits<double>::quiet_NaN();    //[Wt]
    double q_dot_refl_losses = std::numeric_limits<double>::quiet_NaN();    //[Wt]
    double q_dot_rad_losses = std::numeric_limits<double>::quiet_NaN();     //[Wt]
    double q_dot_piping_losses = std::numeric_limits<double>::quiet_NaN();  //[Wt]
    Eigen::MatrixXd E_h;        //[W/m2-K]
    Eigen::MatrixXd E_T;        //[K]
    std::vector<double> m_dot_paths;    //[kg/s]
    std::vector<double> T_out_paths;    //[K]
    if (!rec_is_off) {

        // *************************************************
        // *************************************************
        // Get flux map and interpolate it onto EsolarFlux for each active element

        // get flux map from heliostat field class
        //util::matrix_t<double> flux_map_input = *(inputs.m_flux_map_input);

        double q_dot_inc_design_test = 120.E6;      //[Wt]
        q_dot_inc = q_dot_inc_design_test / m_dni_des * I_bn / field_eff_des * field_eff;   //[Wt]
        double flux_scale = q_dot_inc / m_area_active_total;      //[W/m2]
        //double flux_scale = flux_design / m_dni_des * I_bn / field_eff_des * field_eff;       //[W/m2]

        Eigen::MatrixXd EsolarFlux(mE_areas.rows(), 1);
        EsolarFlux.setConstant(0.0);
        for (size_t i_surf = 0; i_surf < mv_rec_surfs.size(); i_surf++) {
            if (std::isfinite(mv_rec_surfs[i_surf].vertices(0, 0)) && mv_rec_surfs[i_surf].is_active_surf) {
                for (size_t i = 0; i < m_v_elems[i_surf].nrows(); i++) {
                    EsolarFlux(m_surfIDs[i_surf](i, 0), 0) = flux_scale;        //[W/m2]
                }
            }
        }

        // *************************************************
        // *************************************************
        
        // Initialize elemental heat transfer fluid temperatures
        Eigen::MatrixXd E_T_HTF = Eigen::MatrixXd::Zero(m_nElems, 1);       //[K]

        // (using "path" here instead of "pipe" from Matlab code)
        for (size_t i_path = 0; i_path < m_nPaths; i_path++) {
            util::matrix_t<int> FCM;
            FCM = m_FCA[i_path];

            size_t nSteps = FCM.nrows();
            std::vector<double> stepTemp(nSteps);
            double deltaT = (m_T_htf_hot_des - T_salt_cold_in) / (double)(nSteps - 1);
            for (size_t i = 0; i < nSteps; i++) {
                stepTemp[i] = T_salt_cold_in + deltaT * i;      //[K]
            }

            for (size_t i = 0; i < nSteps; i++) {
                for (size_t j = 0; j < FCM.ncols(); j++) {
                    E_T_HTF(FCM(i, j), 0) = stepTemp[i];        //[K]
                }
            }
        }

        // Guess mass flow rate based on design point thermal efficiency
        double m_eta_therm_des = 0.88;
        double q_dot_htf_est = q_dot_inc * m_eta_therm_des;
        double m_dot_path_init_guess = q_dot_htf_est / (m_nPaths * cp_htf * (m_T_htf_hot_des - T_salt_cold_in));     //[kg/s]
        m_dot_paths.resize(m_nPaths, m_dot_path_init_guess);    //[kg/s]

        // Conductivity and convective heat transfer network calculates UA, but energy balance equations divide it by area
        Eigen::MatrixXd E_UA(m_nElems, 1);
        tube_UA_and_deltaP(m_dot_paths, E_T_HTF, E_UA, W_dot_pump);
        Eigen::MatrixXd E_U = E_UA.array() / mE_areas.array();

        // Assign conductance between HTF and each element
        //double UA_elemental = 4000;     //[W/K]
        //E_U.setConstant(std::numeric_limits<double>::quiet_NaN());
        //for (size_t i_surf = 0; i_surf < mv_rec_surfs.size(); i_surf++) {
        //    for (size_t i = 0; i < m_v_elems[i_surf].nrows(); i++) {
        //        E_U(m_surfIDs[i_surf](i, 0), 0) = UA_elemental * (double)mv_rec_surfs[i_surf].is_active_surf / mE_areas(m_surfIDs[i_surf](i, 0), 0);
        //    }
        //}

        Eigen::MatrixXd EqIn = mE_areas.array() * EsolarFlux.array();

        Eigen::MatrixXd eq1 = (-EsolarFlux.array() * mE_areas.array() * mE_rhoSol.array());
        Eigen::MatrixXd eq2(eq1.rows(), eq1.rows());
        Eigen::MatrixXd epsS_square(eq1.rows(), eq1.rows());
        for (size_t i = 0; i < eq2.rows(); i++) {
            for (size_t j = 0; j < eq2.cols(); j++) {
                eq2(i, j) = eq1(i, 0);
                epsS_square(i, j) = mE_epsilonSol(i, 0);
            }
        }
        Eigen::MatrixXd eq3 = (eq2.array() * mE_FHatS.array()).matrix() * mE_epsilonSol;

        Eigen::MatrixXd eq4 = (epsS_square.array() * mE_FHatS.transpose().array()).matrix() * (EsolarFlux.array() * mE_areas.array() * mE_rhoSol.array()).matrix();

        Eigen::MatrixXd EqSolOut = eq3 + eq4;

        Eigen::MatrixXd E_b(EqIn.rows() - 1, 1);
        for (size_t i = 0; i < EqIn.rows() - 1; i++) {
            E_b(i, 0) = (EqIn(i, 0) + EqSolOut(i, 0)) / (mE_epsilonTherm(i, 0) * mE_areas(i, 0) * CSP::sigma) +
                mE_epsilonTherm(m_nElems - 1, 0) * mE_FHatT(i, m_nElems - 1) * pow(T_amb, 4);
        }

        Eigen::MatrixXd E_eye = Eigen::MatrixXd::Identity(m_nElems - 1, m_nElems - 1);

        Eigen::MatrixXd E_epsT_trans = mE_epsilonTherm.transpose();
        E_epsT_trans.conservativeResize(m_nElems - 1, m_nElems - 1);
        for (size_t i = 0; i < m_nElems - 1; i++) {
            for (size_t j = 0; j < m_nElems - 1; j++) {
                E_epsT_trans(i, j) = E_epsT_trans(0, j);
            }
        }

        Eigen::MatrixXd E_A_1 = -1.0 * E_epsT_trans.array() * mE_FHatT.block(0, 0, m_nElems - 1, m_nElems - 1).array();

        Eigen::MatrixXd E_A_2 = mE_FHatT.block(0, 0, m_nElems - 1, m_nElems) * mE_epsilonTherm;
        Eigen::MatrixXd E_A_2_square(E_A_2.rows(), E_A_2.rows());
        for (size_t i = 0; i < E_A_2_square.rows(); i++) {
            for (size_t j = 0; j < E_A_2_square.rows(); j++) {
                E_A_2_square(i, j) = E_A_2(i, 0);
            }
        }
        Eigen::MatrixXd E_A_3 = E_A_2_square.array() * E_eye.array();
        Eigen::MatrixXd E_A = E_A_1 + E_A_3;

        Eigen::MatrixXd E_Tmax1 = E_A.colPivHouseholderQr().solve(E_b);
        Eigen::MatrixXd E_Tmax = Eigen::pow(E_Tmax1.array(), 0.25);

        E_Tmax.conservativeResize(m_nElems, 1);
        E_Tmax(m_nElems - 1, 0) = T_amb;


        // Start iteration between energy and fluid balances
        // Outer loop iterates on nodal/element HTF temperature
        error_T_HTF_node_iter = 10.*tol_T_HTF_node_iter;
        count_T_HTF_node_iter = 0;
        
        while (fabs(error_T_HTF_node_iter) > tol_T_HTF_node_iter && count_T_HTF_node_iter < 100) {

            count_T_HTF_node_iter++;
            
            hbarCorrelation(E_Tmax, T_amb, E_h);

            // solves the energy balance on the cavity reciver
            // system for specified incident solar flux, fluid temperatures, and
            // element view factors.The solution is achieved by first including only
            // radiation heat transfer in the energy balance to obtain an initial
            // guess for the complete iterative energy balance.This solution assumes that
            // -- fluid temperatures at every panel element are known and fixed.
            // -- mass flow rate is fixed

            E_T = E_Tmax;
            error_relmax_T_rec_node_iter = 10.0*tol_rel_T_rec_node_iter;  //[-]
            count_T_rec_node_iter = 0;

            Eigen::MatrixXd E_Tstar;
            while (error_relmax_T_rec_node_iter > tol_rel_T_rec_node_iter && count_T_rec_node_iter < 100) {

                count_T_rec_node_iter++;

                E_Tstar = E_T;

                Eigen::MatrixXd E_Tstar_2 = Eigen::pow(E_Tstar.array(), 2);
                Eigen::MatrixXd E_Tstar_trans = E_Tstar.transpose();
                Eigen::MatrixXd E_Tstar_trans2 = Eigen::pow(E_Tstar_trans.array(), 2);
                Eigen::MatrixXd A_1(m_nElems - 1, m_nElems - 1);
                Eigen::MatrixXd A_2(m_nElems - 1, m_nElems - 1);
                Eigen::MatrixXd A_3(m_nElems - 1, m_nElems);
                Eigen::MatrixXd B_1(m_nElems - 1, 1);
                Eigen::MatrixXd B_2(m_nElems - 1, 1);
                Eigen::MatrixXd B_3(m_nElems - 1, 1);
                for (size_t i = 0; i < m_nElems - 1; i++) {
                    for (size_t j = 0; j < m_nElems - 1; j++) {
                        A_1(i, j) = -mE_epsilonTherm(j, 0) * mE_FHatT(i, j) * (E_Tstar_2(i, 0) + E_Tstar_trans2(0, j)) * (E_Tstar(i, 0) + E_Tstar_trans(0, j));
                        A_2(i, j) = (E_U(i, 0) + E_h(i, 0)) / (mE_epsilonTherm(i, 0) * CSP::sigma);
                    }
                    for (size_t j = 0; j < m_nElems; j++) {
                        A_3(i, j) = (E_Tstar_2(i, 0) + E_Tstar_trans2(0, j)) * (E_Tstar(i, 0) + E_Tstar_trans(0, j)) * mE_FHatT(i, j);
                    }
                    B_1(i, 0) = (EqIn(i, 0) + EqSolOut(i, 0)) / (mE_epsilonTherm(i, 0) * mE_areas(i, 0) * CSP::sigma);
                    B_2(i, 0) = mE_epsilonTherm(m_nElems - 1, 0) * mE_FHatT(i, m_nElems - 1) * T_amb * (E_Tstar_2(i, 0) + pow(T_amb, 2)) * (E_Tstar(i, 0) + T_amb);
                    B_3(i, 0) = (E_h(i, 0) * T_amb + E_U(i, 0) * E_T_HTF(i, 0)) / (mE_epsilonTherm(i, 0) * CSP::sigma);
                }

                Eigen::MatrixXd A_4 = A_3 * mE_epsilonTherm;
                Eigen::MatrixXd A_4_square(A_4.rows(), A_4.rows());
                for (size_t i = 0; i < A_4.rows(); i++) {
                    for (size_t j = 0; j < A_4.rows(); j++) {
                        A_4_square(i, j) = A_4(i, 0);
                    }
                }

                Eigen::MatrixXd A = A_1.array() + (A_2 + A_4_square).array() * E_eye.array();

                Eigen::MatrixXd B = B_1 + B_2 + B_3;

                E_T = A.colPivHouseholderQr().solve(B);
                E_T.conservativeResize(m_nElems, 1);
                E_T(m_nElems - 1, 0) = T_amb;

                // Update convective correlation
                hbarCorrelation(E_T, T_amb, E_h);

                error_relmax_T_rec_node_iter = 0.0;
                for (size_t i = 0; i < m_nElems; i++) {
                    error_relmax_T_rec_node_iter = max(error_relmax_T_rec_node_iter, abs(E_T(i, 0) - E_Tstar(i, 0)) / E_T(i, 0));
                }
            }

            Eigen::MatrixXd E_Q_gain = E_U.array() * mE_areas.array() * (E_T.array() - E_T_HTF.array());   //[W]

            // Using net energy (q_gain) at each node,
            // solve for mass flow rate that achieves the target HTF outlet temperature

            q_gain_net_eb = E_Q_gain.sum();     //[W]
            m_dot_path_init_guess = q_gain_net_eb / (m_nPaths * cp_htf * (m_T_htf_hot_des - T_salt_cold_in));     //[kg/s]
            std::vector<util::matrix_t<double>> mt_T(m_FCA.size());

            m_dot_paths.resize(m_nPaths, m_dot_path_init_guess);    //[kg/s]
            T_out_paths.resize(m_nPaths, std::numeric_limits<double>::quiet_NaN());

            // loop through flow paths
            for (size_t k = 0; k < m_nPaths; k++) {

                error_T_htf_out = tol_abs_T_htf_target * 10.0;
                count_T_htf_out_iter = 0;

                while (fabs(error_T_htf_out) > tol_abs_T_htf_target && count_T_htf_out_iter < 100) {

                    count_T_htf_out_iter++;
                    double mc = m_dot_paths[k] * cp_htf; //[kg/s * J/kg-K] = [W/K]

                    util::matrix_t<int> FCM = m_FCA[k];

                    // padding -1 (matlab uses 0) should have been removed in zigzag method
                    size_t nSteps = FCM.nrows();
                    mt_T[k].resize(nSteps, 1);
                    Eigen::MatrixXd k_b = Eigen::MatrixXd::Zero(nSteps, 1);

                    Eigen::MatrixXd k_a1 = Eigen::MatrixXd::Zero(nSteps, nSteps);
                    k_a1.diagonal(-1).setConstant(1);
                    Eigen::MatrixXd k_a2 = Eigen::MatrixXd::Identity(nSteps, nSteps);
                    Eigen::MatrixXd k_A = mc * (k_a2.array() - k_a1.array());

                    // Extract conduction rates
                    for (size_t i = 0; i < nSteps; i++) {
                        util::matrix_t<int> step_IDs = FCM.row(i);
                        for (size_t j = 0; j < step_IDs.ncols(); j++) {
                            if (step_IDs(0, j) != -1) {
                                k_b(i, 0) += E_Q_gain(step_IDs(0, j));
                            }
                        }
                    }

                    k_b(0, 0) += mc * T_salt_cold_in;

                    // Solve energy balance for this pipe
                    Eigen::MatrixXd E_T_step = k_A.colPivHouseholderQr().solve(k_b);

                    T_out_paths[k] = E_T_step(nSteps - 1, 0);

                    for (size_t ii = 0; ii < nSteps; ii++) {
                        mt_T[k](ii, 0) = E_T_step(ii, 0);
                    }

                    error_T_htf_out = T_out_paths[k] - m_T_htf_hot_des;     //[K]

                    // adjust value of m_dot - usually not needed
                    double q_gain_total_calc = m_dot_paths[k] * cp_htf * (T_out_paths[k] - T_salt_cold_in);       //[W]
                    m_dot_paths[k] = q_gain_total_calc / (cp_htf * (m_T_htf_hot_des - T_salt_cold_in));  //[kg/s]
                }
            }
            
            Eigen::MatrixXd E_T_HTF_calc = Eigen::MatrixXd::Zero(m_nElems, 1);
            // assign fluid temperatures based on solution
            for (size_t k = 0; k < m_nPaths; k++) {
                util::matrix_t<int> FCM = m_FCA[k];

                // padding -1 (matlab uses 0) should have been removed in zigzag method
                size_t nSteps = FCM.nrows();
                for (size_t i = 0; i < nSteps; i++) {
                    util::matrix_t<int> step_IDs = FCM.row(i);
                    for (size_t j = 0; j < step_IDs.ncols(); j++) {
                        E_T_HTF_calc(step_IDs(0, j)) = mt_T[k](i, 0);
                    }
                }
            }

            // Update UA, pressure drop, and pump power calcs
            tube_UA_and_deltaP(m_dot_paths, E_T_HTF, E_UA, W_dot_pump);
            E_U = E_UA.array() / mE_areas.array();

            // Calculate difference between nodal T_HTF and T_HTF_calc
            error_T_HTF_node_iter = 0.0;
            for (size_t i = 0; i < m_nElems; i++) {
                error_T_HTF_node_iter = max(error_T_HTF_node_iter, abs(E_T_HTF(i, 0) - E_T_HTF_calc(i, 0)));
            }

            E_T_HTF = E_T_HTF_calc;
        }
        // energy and fluid balance converged

        m_dot_htf_tot = std::accumulate(m_dot_paths.begin(), m_dot_paths.end(), 0.0);   //[kg/s]
        T_htf_rec_no_tower_losses = std::accumulate(T_out_paths.begin(), T_out_paths.end(), 0.0) / (double)m_nPaths;    //[K]

        // Adjust T_hot_htf_calc for downcomer losses
        // may want to add an option to calculate piping loss with a loss coefficient like external tower model
                // proportionally apply piping loss to each path
        q_dot_piping_losses = m_Q_dot_piping_loss;      //[Wt]
        double delta_T_piping = q_dot_piping_losses / (m_dot_htf_tot * cp_htf);    //[K]

        T_htf_tower_out_calc = T_htf_rec_no_tower_losses - delta_T_piping;        //[K]

        // Calculate converged receiver performance metrics
        q_dot_thermal_tower_out = m_dot_htf_tot * cp_htf * (T_htf_tower_out_calc - T_salt_cold_in);     //[W]
        q_dot_thermal_rec_out_no_tower_losses = m_dot_htf_tot * cp_htf * (T_htf_rec_no_tower_losses - T_salt_cold_in);     //[W]

        // Check that important outputs are realistic
        if (!isfinite(q_dot_thermal_tower_out) || !isfinite(W_dot_pump) || q_dot_thermal_tower_out < 0.0 || T_htf_tower_out_calc > m_T_htf_hot_des + 50.0) {
            rec_is_off = true;
        }

        // Receiver efficiency - does not include riser/downcomer losses
        eta_thermal_calc = q_dot_thermal_rec_out_no_tower_losses / q_dot_inc;       //[-]

        E_h.conservativeResize(m_nElems, 1);
        Eigen::MatrixXd E_q_dot_conv_losses = E_h.array()*mE_areas.array()*(E_T.array() - T_amb);        //[W]
        q_dot_conv_losses = E_q_dot_conv_losses.sum();      //[W]

        q_dot_refl_losses = EqSolOut(m_nElems-1, 0);        //[W]
         
        Eigen::MatrixXd E_q_rad_loss_a1 = mE_FHatT.array()*(Eigen::pow((Eigen::MatrixXd::Zero(m_nElems, m_nElems).colwise() + E_T.col(0)).array(), 4) -
            Eigen::pow((Eigen::MatrixXd::Zero(m_nElems, m_nElems).rowwise() + E_T.col(0).transpose()).array(), 4));
        Eigen::MatrixXd E_q_rad_loss = (E_q_rad_loss_a1*mE_epsilonTherm).array()*mE_epsilonTherm.array()*mE_areas.array()*CSP::sigma;
        q_dot_rad_losses = -E_q_rad_loss(m_nElems-1, 0);    //[W]
    }

    double time_required_su = step / 3600.0;        //[hr]
    double q_startup = 0.0;     //[W-hr]

    if(!rec_is_off){
        switch (input_operation_mode)
        {
        case C_csp_collector_receiver::STARTUP:
            {
                double time_require_su_energy = m_E_su_prev / (m_dot_htf_tot * cp_htf * (T_htf_tower_out_calc - T_salt_cold_in));	//[hr]
                double time_require_su_ramping = m_t_su_prev;

                double time_required_max = fmax(time_require_su_energy, time_require_su_ramping);	//[hr]

                double time_step_hrs = step / 3600.0;		//[hr]

                if (time_required_max > time_step_hrs)		// Can't completely startup receiver in maximum allowable timestep
                {											// Need to advance timestep and try again
                    time_required_su = time_step_hrs;       //[hr]
                    m_mode = C_csp_collector_receiver::STARTUP;
                    q_startup = m_dot_htf_tot * cp_htf * (T_htf_tower_out_calc - T_salt_cold_in) * step / 3600.0; //[W]
                }
                else
                {
                    time_required_su = time_required_max;		//[hr]
                    m_mode = C_csp_collector_receiver::ON;

                    double q_startup_energy_req = m_E_su_prev;	//[W-hr]
                    double q_startup_ramping_req = m_dot_htf_tot * cp_htf * (T_htf_tower_out_calc - T_salt_cold_in) * m_t_su_prev;	//[W-hr]
                    q_startup = fmax(q_startup_energy_req, q_startup_ramping_req);      //[W-hr]
                }

                m_E_su = fmax(0.0, m_E_su_prev - m_dot_htf_tot * cp_htf * (T_htf_tower_out_calc - T_salt_cold_in) * step / 3600.0);
                m_t_su = fmax(0.0, m_t_su_prev - step / 3600.0);
            }

            rec_is_off = true;

            break;

        case C_csp_collector_receiver::ON:

            m_mode = C_csp_collector_receiver::ON;
            q_startup = 0.0;


            break;

        case C_csp_collector_receiver::STEADY_STATE:

            m_mode = C_csp_collector_receiver::STEADY_STATE;

            break;

        }	// End switch() on input_operation_mode-pr
    }
    else {
        // If receiver was off BEFORE startup deductions
        m_mode = C_csp_collector_receiver::OFF;

        // Include here outputs that are ONLY set to zero if receiver completely off, and not attempting to start-up
        W_dot_pump = 0.0;
        // Pressure drops
        DELTAP = 0.0; Pres_D = 0.0; u_coolant = 0.0;
    }

    if (rec_is_off){
        ms_outputs.m_m_dot_salt_tot = 0.0;  //[kg/hr]
        ms_outputs.m_eta_therm = 0.0;       //[-]
        ms_outputs.m_q_conv_sum = 0.0;      //[MWt]
        ms_outputs.m_q_rad_sum = 0.0;       //[MWt]
        ms_outputs.m_Q_thermal = 0.0;       //[MWt]
        ms_outputs.m_T_salt_hot = m_T_htf_cold_des - 273.15;    //[C] convert from K
        ms_outputs.m_q_dot_rec_inc = 0.0;   //[MWt]
        ms_outputs.m_m_dot_ss = 0.0;        //[kg/hr]
        ms_outputs.m_q_dot_ss = 0.0;        //[MWt]
        ms_outputs.m_Q_thermal_csky_ss = 0.0;   //[MWt]
        ms_outputs.m_Q_thermal_ss = 0.0;

        m_od_control = 1.0;
    }
    else {
        ms_outputs.m_m_dot_salt_tot = m_dot_htf_tot * 3600.0;   //[kg/hr] convert from kg/s
        ms_outputs.m_eta_therm = eta_thermal_calc;          //[-]
        ms_outputs.m_q_conv_sum = q_dot_conv_losses/1.E6;   //[MWt] convert from Wt
        ms_outputs.m_q_rad_sum = q_dot_rad_losses/1.E6;     //[MWt] convert from Wt
        ms_outputs.m_Q_thermal = q_dot_thermal_tower_out / 1.E6;      //[MWt] convert from W
        ms_outputs.m_T_salt_hot = T_htf_tower_out_calc - 273.15;               //[C] convert from K
        ms_outputs.m_q_dot_rec_inc = q_dot_inc/1.E6;        //[MWt]
        ms_outputs.m_T_salt_cold = T_salt_cold_in - 273.15; //[C] convert from K
        ms_outputs.m_m_dot_ss = std::numeric_limits<double>::quiet_NaN();
        ms_outputs.m_q_dot_ss = std::numeric_limits<double>::quiet_NaN();
        ms_outputs.m_Q_thermal_csky_ss = std::numeric_limits<double>::quiet_NaN();
        ms_outputs.m_Q_thermal_ss = std::numeric_limits<double>::quiet_NaN();
    }

    ms_outputs.m_W_dot_pump = W_dot_pump;                   //[MWe]
    ms_outputs.m_field_eff_adj = field_eff_adj;             //[-]
    ms_outputs.m_component_defocus = m_od_control;          //[-]
    ms_outputs.m_q_startup = q_startup / 1.E6;              //[MW-hr] convert from [W-hr]
    ms_outputs.m_dP_receiver = std::numeric_limits<double>::quiet_NaN();    //[bar]
    ms_outputs.m_dP_total = std::numeric_limits<double>::quiet_NaN();       //[bar]
    ms_outputs.m_vel_htf = std::numeric_limits<double>::quiet_NaN();        //[m/s]
    ms_outputs.m_f_timestep = std::numeric_limits<double>::quiet_NaN();
    ms_outputs.m_time_required_su = time_required_su * 3600.0;        //[s], convert from hr
    if (q_dot_thermal_tower_out > 0.0) {
        ms_outputs.m_q_dot_piping_loss = q_dot_piping_losses/1.E6;       //[MWt]
    }
    else {
        ms_outputs.m_q_dot_piping_loss = 0.0;       //[MWt]
    }

    ms_outputs.m_q_heattrace = 0.0;
    ms_outputs.m_clearsky = std::numeric_limits<double>::quiet_NaN();

    m_eta_field_iter_prev = field_eff;	//[-]

	return;
}

void C_cavity_receiver::off(const C_csp_weatherreader::S_outputs& weather,
	const C_csp_solver_htf_1state& htf_state_in,
	const C_csp_solver_sim_info& sim_info)
{
    // Don't currently need *any* of these inputs, but if we add recirculation or thermal capacitance it would be helpful to have in place
    m_mode = C_csp_collector_receiver::OFF;

    // Assuming no night recirculation, so... these should be zero
    ms_outputs.m_m_dot_salt_tot = 0.0;		//[kg/hr] convert from kg/s
    ms_outputs.m_eta_therm = 0.0;			//[-] RECEIVER thermal efficiency (includes radiation and convective losses. reflection losses are contained in receiver flux model)
    ms_outputs.m_W_dot_pump = 0.0;			//[MW] convert from W
    ms_outputs.m_q_conv_sum = 0.0;			//[MW] convert from W
    ms_outputs.m_q_rad_sum = 0.0;			//[MW] convert from W
    ms_outputs.m_Q_thermal = 0.0;			//[MW] convert from W
    ms_outputs.m_T_salt_hot = 0.0;			//[C] convert from K
    ms_outputs.m_field_eff_adj = 0.0;		//[-]
    ms_outputs.m_component_defocus = 1.0;	//[-]
    ms_outputs.m_q_dot_rec_inc = 0.0;		//[MW] convert from kW
    ms_outputs.m_q_startup = 0.0;			//[MW-hr] convert from W-hr
    ms_outputs.m_dP_receiver = 0.0;			//[bar] receiver pressure drop, convert from Pa
    ms_outputs.m_dP_total = 0.0;			//[bar] total pressure drop, convert from MPa
    ms_outputs.m_vel_htf = 0.0;				//[m/s]
    ms_outputs.m_T_salt_cold = 0.0;			//[C] convert from K
    ms_outputs.m_m_dot_ss = 0.0;			//[kg/hr] convert from kg/s
    ms_outputs.m_q_dot_ss = 0.0;			//[MW] convert from W
    ms_outputs.m_f_timestep = 0.0;			//[-]
    ms_outputs.m_time_required_su = sim_info.ms_ts.m_step;	//[s], convert from hr in code
    ms_outputs.m_q_dot_piping_loss = 0.0;	//[MWt]
    ms_outputs.m_q_heattrace = 0.0;
    
    ms_outputs.m_clearsky = get_clearsky(weather, sim_info.ms_ts.m_time / 3600.);  // clear-sky DNI (set to actual DNI if actual DNI is higher than computed clear-sky value)
    ms_outputs.m_Q_thermal_csky_ss = 0.0; //[MWt]
    ms_outputs.m_Q_thermal_ss = 0.0; //[MWt]

	return;
}

void C_cavity_receiver::converged()
{
	// Check HTF props?
	//!MJW 9.8.2010 :: Call the property range check subroutine with the inlet and outlet HTF temps to make sure they're in the valid range
	//call check_htf(Coolant,T_salt_hot)
	//call check_htf(Coolant,T_salt_cold)

	if (m_mode == C_csp_collector_receiver::STEADY_STATE)
	{
		throw(C_csp_exception("Receiver should only be run at STEADY STATE mode for estimating output. It must be run at a different mode before exiting a timestep",
			"MSPT receiver converged method"));
	}

    if (m_mode == C_csp_collector_receiver::OFF){
        m_E_su_prev = m_q_rec_des * m_rec_qf_delay;
        m_t_su_prev = m_rec_su_delay;
    }
    else {
        m_E_su_prev = m_E_su;
        m_t_su_prev = m_t_su;
    }

    m_mode_prev = m_mode;

    // Reset call variables
    m_od_control = 1.0;             //[-]
    m_eta_field_iter_prev = 1.0;    //[-]
}


void C_cavity_receiver::calc_pump_performance(double rho_f, double mdot, double ffact, double& PresDrop_calc, double& WdotPump_calc)
{

	

}

double C_cavity_receiver::get_pumping_parasitic_coef()
{
	return  12.3;
}

double C_cavity_receiver::area_proj()
{
	return 1.23;
}
