#ifndef _INDEP_RES_HPP_
#define _INDEP_RES_HPP_

#include <vector>

/*===========================================================================*/
void compute_indep_res_q(
   const int exc_i,
   const int nx, const double dx, const double cl, 
   const std::vector<double> &r_v, 
   const std::vector<double> &phi_f, const std::vector<double> &phi_q,
   std::vector<double> &res_q
);

#endif
