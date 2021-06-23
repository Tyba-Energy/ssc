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

#include "../splinter/Core"
#include "../splinter/LU"
#include "../splinter/Cholesky"
#include "../splinter/QR"
#include "../splinter/SVD"
#include "../splinter/Geometry"
#include "../splinter/Eigenvalues"



class C_cavity_receiver : public C_pt_receiver
{
public:

    class C_rec_surface
    {
    public:
        util::matrix_t<double> vertices;    // (nr, nc) -> (vertex index, dimension (i.e. xyx))
        size_t type;                // mesh type: 0=triangle, 1=quad, 2=single element
        bool is_active_surf;        // True: active surface w/ HTF, False: passive surface no HTF
        bool is_flipRoute;
        double eps_sol;             //[-]
        double eps_therm;           //[-]

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

    // ************************************
    // Design parameters passed in through constructor
    // ************************************
    int m_field_fl;
    util::matrix_t<double> m_field_fl_props;

    double m_hel_stow_deploy;			//[-]

    double m_receiverHeight; //[m] Receiver opening height in meters
    double m_receiverWidth; //[m] Reciever opening width in meters
    double m_topLipHeight;  //[m] Height of top lip in meters
    double m_botLipHeight;  //[m] Height of bottom lip in meters
    double m_e_act_sol;     //[-] Absorbtivity in short wave range for active surfaces
    double m_e_pass_sol;    //[-] Absorbtivity in short wave range for passive surfaces
    double m_e_act_therm;   //[-] Emissivity in long wave range for active surfaces
    double m_e_pass_therm;  //[-] Emissivity in long wave range for passive surfaces

    double m_elemSize;      //
    // ************************************

    // ************************************
    // Calculated stored parameters
    std::vector<C_rec_surface> mv_rec_surfs;    // vector of surface classes for each surface in cavity model
    std::vector<util::matrix_t<int>> m_v_elems; // each vector index is a surface; each row lists the nodes that define a mesh element
    util::matrix_t<double> m_nodesGlobal;       // each row lists x,y,z coordinates of each node

    util::matrix_t<int> m_elements;             // global element tracker; each row, indexed by m_surfIDs, lists the nodes that define a mesh element
    std::vector<util::matrix_t<int>> m_surfIDs; // global element count for each surface
    util::matrix_t<double> m_areas;             // global element areas, each row indexed by m_surfIDs
    Eigen::MatrixXd mE_areas;                   // global element areas, each row indexed by m_surfIDs
    util::matrix_t<double> m_centroids;         // global element centroids, each row indexed by m_surfIDs
    size_t m_nElems;                            // global element centroids, each row indexed by m_surfIDs

    util::matrix_t<double> m_epsilonSol;        // global element solar emissivity 
    Eigen::MatrixXd mE_epsilonSol;              // global element solar emissivity
    util::matrix_t<double> m_epsilonTherm;      // global element thermal emissivity
    Eigen::MatrixXd mE_epsilonTherm;            // global element thermal emissivity

    std::vector<util::matrix_t<int>> m_FCA;     // fluid connectivity array

    util::matrix_t<double> m_F;                 // view factors

    util::matrix_t<double> m_FHatS;             // FHat solar
    Eigen::MatrixXd mE_FHatS;                   // FHat solar
    util::matrix_t<double> m_FHatT;             // FHat thermal
    Eigen::MatrixXd mE_FHatT;                   // FHat thermal

    Eigen::MatrixXd mE_rhoSol;                  // global element solar reflectivity
    Eigen::MatrixXd mE_rhoTherm;                // global element thermal reflectivity

    // ************************************
    // Call variables
    double m_eta_field_iter_prev;	//[-] Efficiency from heliostat on last iteration. Maybe change if CR gets defocus signal from controller
    double m_od_control;            //[-]

    // ************************************
    // State variables
    // m_mode_prev and m_mode are members of parent class
    double m_E_su_prev;         //[W-hr] Startup energy required at end of previous timestep
    double m_E_su;              //[W-hr] Startup energy required calculated at end of current timestep

    double m_t_su;          //[hr] Startup time requirement at end of previous timestep
    double m_t_su_prev;     //[hr] Startup time requirement calculated at end of current timestep

public:

	// Methods
	C_cavity_receiver(double hel_stow_deploy /*-*/, double T_htf_hot_des /*K*/, double q_dot_rec_des /*MWt*/,
        double rec_qf_delay /*-*/, double rec_su_delay /*hr*/, int field_fl /*-*/, util::matrix_t<double> field_fl_props,
        double rec_height /*m*/, double rec_width /*m*/, double toplip_height /*m*/, double botlip_height /*m*/,
        double eps_active_sol /*-*/, double eps_passive_sol /*-*/, double eps_active_therm /*-*/, double eps_passive_therm /*-*/,
        double elemSize );

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

    void genOctCavity();

    void meshGeometry();

    void makeGlobalElems();

    void surfValuesToElems();

    void zigzagRouting(size_t n_steps);

    void VFMatrix();

    void FHatMatrix(const util::matrix_t<double>& eps,
        util::matrix_t<double>& F_hat, util::matrix_t<double>& rho,
        Eigen::MatrixXd& E_F_hat, Eigen::MatrixXd& E_rho);

    void matrixt_to_eigen(const util::matrix_t<double>& matrixt,
        Eigen::MatrixXd& eigenx);

    void hbarCorrelation(const Eigen::MatrixXd& T, double T_inf, Eigen::MatrixXd& h);

    void interpSolarFlux(const util::matrix_t<double>& fluxDist);

    void edgePairParameters(const util::matrix_t<double>& Po, const util::matrix_t<double>& Pf, const util::matrix_t<double>& Qo, const util::matrix_t<double>& Qf,
        double& D, util::matrix_t<double>& sOrigin, util::matrix_t<double>& sHat, util::matrix_t<double>& lHat, util::matrix_t<double>& lOrigin, bool& skew);

    void viewFactor(const util::matrix_t<double>& a, const util::matrix_t<double>& b, double& F_AB, double& F_BA);

    double fParallel(double s, double l, double d);

    double f_skew(double s, double l, double alpha, double cosAlpha, double sinAlpha, double d);

    double imagLi_2(double mag, double angle);

    double Cl(double theta);

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

    double min_column_val(const util::matrix_t<double>& a, size_t n_c);

    double max_column_val(const util::matrix_t<double>& a, size_t n_c);

    int max_int_first_column(const util::matrix_t<int>& a);

    bool are_rows_equal(const util::matrix_t<double>& a, const util::matrix_t<double>& b, int i_row);

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

    void polygon_normal_and_area(const util::matrix_t<double>& poly_a,
        util::matrix_t<double>& norm_vect, double& area, int& n_rows);

};

#endif // __csp_solver_cavity_receiver_

