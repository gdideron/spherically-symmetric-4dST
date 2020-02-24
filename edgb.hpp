#ifndef _EDGB_HPP_
#define _EDGB_HPP_

#include <vector>
using std::vector;

#include "radial_pts.hpp"
#include "field.hpp"

/*===========================================================================*/
class EdGB
{
public:
	EdGB(	const double dt, const double dx, const double cl, const int nx,
		const double mu, const double la,
		const double gbc1, const double gbc2,
		const vector<double> &rp
	);
	~EdGB(void);

	void time_step(
		const int exc_i,
		Field &al, Field &ze, Field &phi_f, Field &phi_p, Field &phi_q
	);
	void solve_metric_fields(
		const int exc_i,
		const Field &phi_f, const Field &phi_p, const Field &phi_q,
		Field &al, Field &ze 
	);

	void compute_radial_characteristics(
		int &exc_i,
		const vector<double> &al_v, const vector<double> &ze_v,
		const vector<double>  &f_v, 
		const vector<double>  &p_v, const vector<double>  &q_v,
		vector<double> &ingoing,
		vector<double> &outgoing
	);
	void compute_eom_rr(
		const int exc_i,
		const vector<double> &al_v, const vector<double> &ze_v,
		const vector<double>  &f_v, 
		const vector<double>  &p_v, const vector<double>  &q_v,
		vector<double> &eom_rr
	);
	void compute_ncc(
		const int exc_i,
		const int nx, const double dx, const double cl, 
		const vector<double> &r_v, 
		const vector<double> &al, const vector<double> &ze,
		const vector<double> &p,  const vector<double> &q,
		vector<double> &ncc
	);
private:
	const double dt;
	const double dx;
	const int nx;

	const double mu;
	const double la;

	const double gbc1;
	const double gbc2;

	const double cl;

	vector<double> r_v;

	double ze_free_k;
	double ze_free_k1;
	double ze_free_k2;
	double ze_free_k3;
	double ze_free_k4;

	vector<double> f_k;
	vector<double> p_k;
	vector<double> q_k;

	vector<double> f_k1;
	vector<double> p_k1;
	vector<double> q_k1;

	vector<double> f_k2;
	vector<double> p_k2;
	vector<double> q_k2;

	vector<double> f_k3;
	vector<double> p_k3;
	vector<double> q_k3;

	vector<double> f_k4;
	vector<double> p_k4;
	vector<double> q_k4;

	double ingoing_c;
	double outgoing_c;

	vector<double> V_v;
	vector<double> Vp_v;

	vector<double> W_v;
	vector<double> Wp_v;
	vector<double> Wpp_v;
/*---------------------------------------------------------------------------*/
	void rescale(
		vector<double> &vec
	);
	double compute_p_k(
		const double r,
		const double Al, const double Ze,
		const double P,  const double Q,
		const double V,  const double Vp, const double Wp, const double Wpp,
		const double r_Der_Al,
		const double r_Der_Ze,
		const double r_Der_P,
		const double r_Der_Q
	);
	void compute_alma_ki(
		const int loc,
		const int val,
		const vector<double> &r3,
		const double al,
		const double ze,
		const vector<double> &p3,
		const vector<double> &q3,
		const vector<double> &r_Der_p3,
		const vector<double> &r_Der_q3
	);
	double compute_ze_free_k(
		const double r,
		const double Al, const double Ze,
		const double P,  const double Q,
		const double V,  const double Vp, const double Wp, const double Wpp,
		const double r_Der_Al,
		const double r_Der_Ze,
		const double r_Der_P,
		const double r_Der_Q
	);
	void solve_for_metric_relaxation(
		const int exc_i,
		const vector<double> &f_v,
		const vector<double> &p_v,
		const vector<double> &q_v,
		vector<double> &al_v,
		vector<double> &ze_v
	);
	void compute_scalar_potentials(const vector<double> &f_v);
	double compute_res_al(
		const double r,
		const double Al, const double Ze,
		const double P,  const double Q,
		const double V,  const double Wp, const double Wpp,
		const double r_Der_Al,
		const double r_Der_P,
		const double r_Der_Q
	);
	double compute_res_ze(
		const double r,
		const double Ze,
		const double P,  const double Q,
		const double V,  const double Wp, const double Wpp,
		const double r_Der_Ze,
		const double r_Der_P,
		const double r_Der_Q
	);
	void compute_fpq_ki(
		const int val,
		const int exc_i,
		const vector<double> &al,
		const vector<double> &ze,
		const vector<double> &f, 
		const vector<double> &p, 
		const vector<double> &q
	);

	int compute_radial_characteristic(
		const double r,
		const double Al, const double Ze,
		const double P,  const double Q,
		const double V,  const double Vp, const double Wp, const double Wpp,
		const double r_Der_Al,
		const double r_Der_Ze,
		const double r_Der_P,
		const double r_Der_Q
	);
};

#endif
