#ifndef _EDGB_HPP_
#define _EDGB_HPP_

#include <vector>

#include "radial_pts.hpp"
#include "field.hpp"

/*===========================================================================*/
class EdGB
{
public:
   EdGB(const double dt, const double dx, const double cl, const int nx,
      const double V_1,  const double V_2,  const double V_3,  const double V_4,
      const double Al_0, 
      const double Al_1, const double Al_2, const double Al_3, const double Al_4,
      const double Be_1, const double Be_2, const double Be_3, const double Be_4,
      const std::vector<double> &rp
   );
   ~EdGB(void);

   void time_step(
      const int exc_i,
      Field &N, Field &S, Field &phi_f, Field &phi_p, Field &phi_q
   );
   void solve_metric_fields(
      const int exc_i,
      const Field &phi_f, const Field &phi_p, const Field &phi_q,
      Field &N, Field &S 
   );

   void compute_radial_characteristics(
      int &exc_i,
      const std::vector<double>  &N_v, const std::vector<double> &S_v,
      const std::vector<double>  &f_v, 
      const std::vector<double>  &p_v, const std::vector<double>  &q_v,
      std::vector<double> &ingoing,
      std::vector<double> &outgoing
   );
   void compute_eom_rr(
      const int exc_i,
      const std::vector<double> &N_v, const std::vector<double> &S_v,
      const std::vector<double> &f_v, 
      const std::vector<double> &p_v, const std::vector<double>  &q_v,
      std::vector<double> &eom_rr
   );
   void compute_ncc(
      const int exc_i,
      const int nx, const double dx, const double cl, 
      const std::vector<double> &r_v, 
      const std::vector<double> &N, const std::vector<double> &S,
      const std::vector<double> &f, 
      const std::vector<double> &p, const std::vector<double> &q,
      std::vector<double> &ncc
   );
private:
   const double dt;
   const double dx;
   const int nx;

   const double V_1;
   const double V_2;
   const double V_3;
   const double V_4;

   const double Al_0;
   const double Al_1;
   const double Al_2;
   const double Al_3;
   const double Al_4;

   const double Be_1;
   const double Be_2;
   const double Be_3;
   const double Be_4;

   const double cl;

   std::vector<double> r_v;

   double S_free_k;
   double S_free_k1;
   double S_free_k2;
   double S_free_k3;
   double S_free_k4;

   std::vector<double> f_k;
   std::vector<double> p_k;
   std::vector<double> q_k;

   std::vector<double> f_k1;
   std::vector<double> p_k1;
   std::vector<double> q_k1;

   std::vector<double> f_k2;
   std::vector<double> p_k2;
   std::vector<double> q_k2;

   std::vector<double> f_k3;
   std::vector<double> p_k3;
   std::vector<double> q_k3;

   std::vector<double> f_k4;
   std::vector<double> p_k4;
   std::vector<double> q_k4;

   double ingoing_c;
   double outgoing_c;

   std::vector<double> V_v;
   std::vector<double> Vp_v;

   std::vector<double> Al_v;
   std::vector<double> Alp_v;

   std::vector<double> Be_v;
   std::vector<double> Bep_v;
   std::vector<double> Bepp_v;
/*---------------------------------------------------------------------------*/
   void rescale(
      std::vector<double> &vec
   );
   double compute_p_k(
      const double r,
      const double N, const double S,
      const double P, const double Q,
      const double V, const double Vp, 
      const double Al,  const double Alp,
      const double Bep, const double Bepp,
      const double r_Der_N,
      const double r_Der_S,
      const double r_Der_P,
      const double r_Der_Q
   );
   void compute_alma_ki(
      const int loc,
      const int val,
      const std::vector<double> &r3,
      const double N,
      const double S,
      const std::vector<double> &p3,
      const std::vector<double> &q3,
      const std::vector<double> &r_Der_p3,
      const std::vector<double> &r_Der_q3
   );
   double compute_S_free_k(
      const double r,
      const double N, const double S,
      const double P, const double Q,
      const double V, const double Vp, 
      const double Al,  const double Alp,
      const double Bep, const double Bepp,
      const double r_Der_N,
      const double r_Der_S,
      const double r_Der_P,
      const double r_Der_Q
   );
   void solve_for_metric_relaxation(
      const int exc_i,
      const std::vector<double> &f_v,
      const std::vector<double> &p_v,
      const std::vector<double> &q_v,
      std::vector<double> &N_v,
      std::vector<double> &S_v
   );
   void compute_potentials(const std::vector<double> &f_v);

   double compute_res_N(
      const double r,
      const double N, const double S,
      const double P, const double Q,
      const double V,  
      const double Al,
      const double Bep, const double Bepp,
      const double r_Der_N,
      const double r_Der_P,
      const double r_Der_Q
   );
   double compute_res_S(
      const double r,
      const double S,
      const double P,  const double Q,
      const double V,
      const double Al,
      const double Bep, const double Bepp,
      const double r_Der_S,
      const double r_Der_P,
      const double r_Der_Q
   );
   void compute_fpq_ki(
      const int val,
      const int exc_i,
      const std::vector<double> &N,
      const std::vector<double> &S,
      const std::vector<double> &f, 
      const std::vector<double> &p, 
      const std::vector<double> &q
   );

   int compute_radial_characteristic(
      const double r,
      const double N, const double S,
      const double P, const double Q,
      const double V, const double Vp, 
      const double Al,  const double Alp,
      const double Bep, const double Bepp, 
      const double r_Der_N,
      const double r_Der_S,
      const double r_Der_P,
      const double r_Der_Q
   );
};

#endif
