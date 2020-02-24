#ifndef _INDEP_RES_HPP_
#define _INDEP_RES_HPP_

#include <vector>
using std::vector;

/*===========================================================================*/
void compute_indep_res_q(
	const int exc_i,
	const int nx, const double dx, const double cl, 
	const vector<double> &r_v, 
	const vector<double> &phi_f, const vector<double> &phi_q,
	vector<double> &res_q
);

#endif
