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

#ifndef __csp_solver_cavity_receiver_
#define __csp_solver_cavity_receiver_

#include "csp_solver_util.h"

#include "htf_props.h"
#include "csp_solver_core.h"
#include "csp_solver_pt_receiver.h"

class C_cavity_receiver : public C_pt_receiver
{
public:

    class C_rec_surface
    {
    public:
        util::matrix_t<double> vertices;    // (nr, nc) -> (vertices, xyx)
        double e_size;
        size_t type;
        bool is_active_surf;

        C_rec_surface()
        {
            e_size = 2.5;
        }
    };

    enum surf_order
    {
        PANEL1,
        PANEL2,
        PANEL3,
        PANEL4,
        FLOOR,
        COVER,
        TOPLIP,
        BOTLIP,
        APERTURE
    };

private:

    double receiverHeight; //[m] Receiver opening height in meters
    double receiverWidth; //[m] Reciever opening width in meters
    double topLipHeight;  //[m] Height of top lip in meters
    double botLipHeight;  //[m] Height of bottom lip in meters
    double e_act_sol;     //[-] Absorbtivity in short wave range for active surfaces
    double e_pass_sol;    //[-] Absorbtivity in short wave range for passive surfaces
    double e_act_therm;   //[-] Emissivity in long wave range for active surfaces
    double e_pass_therm;  //[-] Emissivity in long wave range for passive surfaces
    //double T_HTFin;       // Inlet heat transfer fluid temperature
    //double T_HTFout;      // Outlet heat transfer fluid temperature
    //double T_inf;         // Temperature of surroundings
    //double UA_elemental;  // Specified conductance from HTF to each element
    //double flux_elemental;// Specified incident solar flux on each element
    //double h;             // Convective heat transfer coefficients per element

    std::vector<C_rec_surface> mv_rec_surfs;
    std::vector<util::matrix_t<int>> m_v_elems;
    util::matrix_t<double> m_nodesGlobal;

    util::matrix_t<int> m_elements;
    std::vector<util::matrix_t<int>> m_surfIDs;
    util::matrix_t<double> m_areas;
    util::matrix_t<double> m_centroids;
    int m_nElems;

    util::matrix_t<bool> m_globalValues;

public:

	// Methods
	C_cavity_receiver();

	~C_cavity_receiver() {};

	virtual void init();

	virtual void call(const C_csp_weatherreader::S_outputs& weather,
		const C_csp_solver_htf_1state& htf_state_in,
		const C_pt_receiver::S_inputs& inputs,
		const C_csp_solver_sim_info& sim_info);

	virtual void off(const C_csp_weatherreader::S_outputs& weather,
		const C_csp_solver_htf_1state& htf_state_in,
		const C_csp_solver_sim_info& sim_info);

	virtual void converged();

	void calc_pump_performance(double rho_f, double mdot, double ffact, double& PresDrop_calc, double& WdotPump_calc);

	virtual double get_pumping_parasitic_coef();

	virtual double area_proj();

    void genOctCavity(double height /*m*/, double width /*m*/,
        double liptop /*m*/, double lipbot /*m*/);

    void meshGeometry();

    void makeGlobalElems();

    void surfValuesToElems();

    void meshMapped(const util::matrix_t<double>& poly, double elemSize,
        util::matrix_t<double>& nodes, util::matrix_t<int>& quads);

    void meshPolygon(const util::matrix_t<double>& poly, double elemSize);

    void crossproduct(const util::matrix_t<double>&, const util::matrix_t<double>&, util::matrix_t<double>& cross);

    void norm3Dvect(const util::matrix_t<double>&, util::matrix_t<double>& norm_vect);

    double mag_vect(const util::matrix_t<double>& vector_in);

    double dotprod3D(const util::matrix_t<double>&, const util::matrix_t<double>&);

    void flipup(const util::matrix_t<double>& a, util::matrix_t<double>& b);

    void sumcolumns(const util::matrix_t<double>&, util::matrix_t<double>&);

    void sum_int_columns(const util::matrix_t<int>& a, util::matrix_t<int>& summed);

    void diffrows(const util::matrix_t<double>& a, const util::matrix_t<double>& b, util::matrix_t<double>& a_less_b);

    void add_vect_rows(const util::matrix_t<double>& a, const util::matrix_t<double>& b, util::matrix_t<double>& a_plus_b);

    void scale_vect(const util::matrix_t<double>& a, double scale, util::matrix_t<double>& out_vect);

    void add_constant_to_each_element(int val, util::matrix_t<int>& a);

    void ave_columns(const util::matrix_t<double>&, util::matrix_t<double>&);

    double max_row_value(const util::matrix_t<double>& a);

    int max_row_int_value(const util::matrix_t<int>& a);

    double min_val_first_colum(const util::matrix_t<double>& a);

    int max_int_first_column(const util::matrix_t<int>& a);

    void min_max_vects_from_columns(const util::matrix_t<double>& a, util::matrix_t<double>& max_vect, util::matrix_t<double>& min_vect);

    void transpose_matrix_t(const util::matrix_t<double>& a, util::matrix_t<double>& b);

    void transpose_int_matrix_t(const util::matrix_t<int>& a, util::matrix_t<int>& b);

    void to2D(const util::matrix_t<double>& poly, const util::matrix_t<double>& center,
        const util::matrix_t<double>& normal, const util::matrix_t<double>& xaxis,
        util::matrix_t<double>& poly_xy, util::matrix_t<double>& poly_rt);

    void map(const util::matrix_t<double>& poly2D, double elemSize,
        util::matrix_t<double>& nodes, util::matrix_t<int>& quads);

    void to3D(const util::matrix_t<double>& poly_xy, const util::matrix_t<double>& origin,
        const util::matrix_t<double>& normal, const util::matrix_t<double>& xaxis,
        util::matrix_t<double>& poly3d);

    // triMesh2D(fd, fh, h0, bbox, pfix, varargin)
    void triMesh2D(double h0, const util::matrix_t<double>& bbox, const util::matrix_t<double>& pfix,
        const util::matrix_t<double>& poly_2D);

    void pointToPoly(const util::matrix_t<double>& p, const util::matrix_t<double>& POLY,
                util::matrix_t<double>& d);

    double pointToLine(const util::matrix_t<double>& p, const util::matrix_t<double>& a,
        const util::matrix_t<double>& b);

    void inpolygon(const util::matrix_t<double>& p_x, const util::matrix_t<double>& p_y,
        const util::matrix_t<double>& poly_x, const util::matrix_t<double>& poly_y,
        util::matrix_t<int>& is_in_polygon);


};

#endif // __csp_solver_cavity_receiver_
