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
      const double mu, const double la,
      const double gbc1, const double gbc2,
      const std::vector<double> &rp
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
      const std::vector<double> &al_v, const std::vector<double> &ze_v,
      const std::vector<double>  &f_v, 
      const std::vector<double>  &p_v, const std::vector<double>  &q_v,
      std::vector<double> &ingoing,
      std::vector<double> &outgoing
   );
   void compute_eom_rr(
      const int exc_i,
      const std::vector<double> &al_v, const std::vector<double> &ze_v,
      const std::vector<double>  &f_v, 
      const std::vector<double>  &p_v, const std::vector<double>  &q_v,
      std::vector<double> &eom_rr
   );
   void compute_ncc(
      const int exc_i,
      const int nx, const double dx, const double cl, 
      const std::vector<double> &r_v, 
      const std::vector<double> &al, const std::vector<double> &ze,
      const std::vector<double> &p,  const std::vector<double> &q,
      std::vector<double> &ncc
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

   std::vector<double> r_v;

   double ze_free_k;
   double ze_free_k1;
   double ze_free_k2;
   double ze_free_k3;
   double ze_free_k4;

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

   std::vector<double> A_v;

   std::vector<double> B_v;
   std::vector<double> Bp_v;
   std::vector<double> Bpp_v;
/*---------------------------------------------------------------------------*/
   void rescale(
      std::vector<double> &vec
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
      const std::vector<double> &r3,
      const double al,
      const double ze,
      const std::vector<double> &p3,
      const std::vector<double> &q3,
      const std::vector<double> &r_Der_p3,
      const std::vector<double> &r_Der_q3
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
      const std::vector<double> &f_v,
      const std::vector<double> &p_v,
      const std::vector<double> &q_v,
      std::vector<double> &al_v,
      std::vector<double> &ze_v
   );
   void compute_scalar_potentials(const std::vector<double> &f_v);
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
      const std::vector<double> &al,
      const std::vector<double> &ze,
      const std::vector<double> &f, 
      const std::vector<double> &p, 
      const std::vector<double> &q
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
