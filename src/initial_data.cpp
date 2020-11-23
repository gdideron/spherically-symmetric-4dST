#include <iomanip>
using std::setprecision;
using std::setw;
#include <iostream>
using std::cout;
using std::endl;
#include <cassert> 
#include <cmath>
using std::exp;
using std::pow;
#include <vector>
using std::vector;
#include "initial_data.hpp"
#include "sim_params.hpp"
#include "field.hpp"
#include "edgb.hpp"
/*===========================================================================*/
void set_initial_data(
   const Sim_params &sp,
   const vector<double> &rvec,
   Field &N,Field &S,
   Field &f, Field &p, Field &q)
{
   double bh_mass= sp.bh_mass;
/*-------------------------------------------------------------------------*/
   if (sp.initial_data_type=="bump_with_bh") {
      double amp= sp.amp;
      double r_l= sp.r_l;
      double r_u= sp.r_u;

      double max_vN= 0;

      assert(r_u>r_l);
      for (int i=sp.initial_exc_i; i<sp.nx-1; ++i) {
         double r= rvec[i];
         if ((r>r_l)
         &&  (r<r_u)
         ) {
            double bump= exp(-1./(r_u-r))*exp(-1./(r-r_l)); 

            f.n[i]= pow(r-r_l,2)*pow(r_u-r,2)*bump; 
            p.n[i]= 0;
            q.n[i]= (
            2*(r-r_l)*pow(r_u-r,2)
            -	2*pow(r-r_l,2)*(r_u-r)
            +	pow(r_u-r,2)
            -	pow(r-r_l,2)
            )*bump;
         } else {
            f.n[i]= 0;
            p.n[i]= 0;
            q.n[i]= 0;
         }
         max_vN= (fabs(f.n[i])>max_vN) ? fabs(f.n[i]) : max_vN;

         N.n[i]= 1;
         S.n[i]= pow(2*bh_mass/r,0.5);
      }
      /* rescale so amp is actuN maximum val */
      for (int i=0; i<sp.nx-1; ++i) {
         f.n[i]*= amp/max_vN;
         q.n[i]*= amp/max_vN;
         p.n[i]*= amp/max_vN;
      }
/*-------------------------------------------------------------------------*/
   } else 	
   if (sp.initial_data_type=="scalarized_bh") {

      double charge= sp.charge;
      double mu= sp.mu;

      for (int i=sp.initial_exc_i; i<sp.nx-1; ++i) {
         double r= rvec[i];

         f.n[i]= 
            (charge/r)*exp(-mu*(r-3*bh_mass))
         ; 
         q.n[i]= 
         -  (charge/pow(r,2))*exp(-mu*(r-3*bh_mass))
         -  mu*(charge/r)*exp(-mu*(r-3*bh_mass))
         ;
         N.n[i]= 1;
         S.n[i]= pow(2*bh_mass/r,0.5);

         p.n[i]= -S.n[i]*q.n[i];
      }
/*-------------------------------------------------------------------------*/
   } else {
      cout << "ERROR: initial_data_type " << sp.initial_data_type << " did not match" << endl;
      std::quick_exit(0);
   }
/*-------------------------------------------------------------------------*/
   for (int i=0; i<sp.nx; ++i) {
      f.inter_2[i]= f.n[i];
      f.inter_3[i]= f.n[i];
      f.inter_4[i]= f.n[i];
      f.np1[i]=     f.n[i];

      p.inter_2[i]= p.n[i];
      p.inter_3[i]= p.n[i];
      p.inter_4[i]= p.n[i];
      p.np1[i]=     p.n[i];

      q.inter_2[i]= q.n[i];
      q.inter_3[i]= q.n[i];
      q.inter_4[i]= q.n[i];
      q.np1[i]=     q.n[i];

      N.inter_2[i]= N.n[i];
      N.inter_3[i]= N.n[i];
      N.inter_4[i]= N.n[i];
      N.np1[i]=     N.n[i];

      S.inter_2[i]= S.n[i];
      S.inter_3[i]= S.n[i];
      S.inter_4[i]= S.n[i];
      S.np1[i]=     S.n[i];
   }
/*-------------------------------------------------------------------------*/
   return;
}
/*===========================================================================*/	
void time_symmetric(EdGB &edgb,
   const Sim_params &sp,
   const Field &f, const Field &q,
   Field &N, Field &S,
   Field &p)
{
   const int exc_i= sp.initial_exc_i;
   const int pt= sp.initial_exc_i+20;
   assert(pt<sp.nx-2);

   double old_p= p.n[pt];

   int iter=0;

   do {
      old_p= p.n[pt];

      edgb.solve_metric_fields(exc_i,
	 f, p, q, N, S
      );

      for (int i=0; i<sp.nx; ++i) {
         p.n[i]= -S.n[i]*q.n[i];

         p.inter_2[i]= p.n[i];
         p.inter_3[i]= p.n[i];
         p.inter_4[i]= p.n[i];
         p.np1[i]=     p.n[i];
         }
         iter++;
         if (iter>1e3) {
            cout<<"ERROR(time_symmetric): iter>1e3"<<endl;
            std::quick_exit(0);
         }
   } while ((p.n[pt]-old_p)>1e-3);

   return;
}
