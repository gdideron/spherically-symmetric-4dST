#include <cmath>
using std::pow ;
#include <vector>
using std::vector;

#include "indep_res.hpp"
#include "fd_stencils.hpp"

/*===========================================================================*/
void compute_indep_res_q(
   const int exc_i,
   const int nx, const double dx, const double cl, 
   const vector<double> &r_v, 
   const vector<double> &phi_f, const vector<double> &phi_q,
   vector<double> &res_q)
{
   for (int i=exc_i+2; i<nx-3; ++i) {
      double cf= 1/(1+r_v[i]/cl);
      res_q[i]= 
         pow(cf,2)*Dx_ptc_4th(phi_f[i+2],phi_f[i+1],phi_f[i-1],phi_f[i-2],dx)
      -	phi_q[i]
      ;
   }
   return;
}
