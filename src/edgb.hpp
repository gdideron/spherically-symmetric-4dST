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
      const double Be_exp2,
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
      const std::vector<double>  &p_v, const std::vector<double>  &q_v,
      std::vector<double> &ingoing,
      std::vector<double> &outgoing
   );
   void compute_eom_rr(
      const int exc_i,
      const std::vector<double> &N_v, const std::vector<double> &S_v,
      const std::vector<double> &p_v, const std::vector<double>  &q_v,
      const std::vector<double> &S_nm1_v, const std::vector<double>  &p_nm1_v,
      std::vector<double> &eom_rr
   );
   void compute_ncc(
      const int exc_i,
      const int nx, const double dx, const double cl, 
      const std::vector<double> &r_v, 
      const std::vector<double> &N, const std::vector<double> &S,
      const std::vector<double> &p, const std::vector<double> &q,
      std::vector<double> &ncc
   );
   void compute_det_effective_metric(
      const int exc_i,
      const int nx, const double dx, const double cl, 
      const vector<double> &r_v, 
      const vector<double> &n_v, const vector<double> &s_v,
      const vector<double> &p_v, const vector<double> &q_v,
      vector<double> &det
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

   const double Be_exp2;

   const double cl;

   std::vector<double> r_v;

   double S_free_k1;
   double S_free_k2;
   double S_free_k3;
   double S_free_k4;

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

   std::vector<double> res_n_v;
   std::vector<double> res_s_v;
   std::vector<double> jac_n_v;
   std::vector<double> jac_s_v;

   std::vector<double> r_avg;
   std::vector<double> dr_avg;

   std::vector<double> p_avg; 
   std::vector<double> q_avg; 
   std::vector<double> n_avg; 
   std::vector<double> s_avg; 

   std::vector<double> V_avg; 
   std::vector<double> Al_avg; 
   std::vector<double> Bep_avg;
   std::vector<double> Bepp_avg;

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
   ) const;
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
   ) const;
   void solve_for_metric_relaxation(
      const int exc_i,
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
   ) const;
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
   ) const;
   void compute_fpqS_ki(
      const int exc_i,
      const std::vector<double> &N,
      const std::vector<double> &S,
      const std::vector<double> &p, 
      const std::vector<double> &q,
      std::vector<double> &f_k,
      std::vector<double> &p_k,
      std::vector<double> &q_k,
      double &S_free_k
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

/*===========================================================================*/
/* inlined function calls: do messy computations for time derivatives
 * and residuals */
/*===========================================================================*/
inline double EdGB::compute_p_k(
   const double r,
   const double N, const double S,
   const double P, const double Q,
   const double V, const double Vp, 
   const double Al,  const double Alp,
   const double Bep, const double Bepp,
   const double r_Der_N,
   const double r_Der_S,
   const double r_Der_P,
   const double r_Der_Q) const
{
   double Qr= Q/r;
   double Sr= S/r;

   double p_k= 
      2*Qr*pow(r,4)*N - 32*Bep*pow(Qr,2)*pow(r,4)*N - 2*Al*pow(Qr,3)*pow(r,6)*N + 32*Al*Bep*pow(Qr,4)*pow(r,6)*N - 128*Al*pow(Bep,2)*pow(Qr,5)*pow(r,6)*N - 5*Al*Bep*pow(Qr,4)*pow(r,8)*pow(Sr,2)*N - 16*Bep*Bepp*pow(Qr,4)*pow(r,8)*pow(Sr,2)*N + 16*Al*pow(Bep,2)*pow(Qr,5)*pow(r,8)*pow(Sr,2)*N + 4*Bep*pow(r,4)*pow(Sr,4)*N - 32*pow(Bep,2)*pow(Qr,3)*pow(r,4)*(-4 + pow(r,2)*pow(Sr,2))*N + 24*Bep*pow(Qr,6)*pow(r,8)*(-2*Alp*Bep + Al*Bepp*pow(r,2)*pow(Sr,2))*N - pow(r,4)*Vp*N + 16*Bep*Qr*pow(r,4)*Vp*N - 64*pow(Bep,2)*pow(Qr,2)*pow(r,4)*Vp*N + 2*pow(r,4)*Sr*N*P - 60*Bep*Qr*pow(r,4)*Sr*N*P + 384*pow(Bep,2)*pow(Qr,2)*pow(r,4)*Sr*N*P - 2*Al*pow(Qr,2)*pow(r,6)*Sr*N*P + 60*Al*Bep*pow(Qr,3)*pow(r,6)*Sr*N*P - 384*Al*pow(Bep,2)*pow(Qr,4)*pow(r,6)*Sr*N*P + 32*Al*Bep*Bepp*pow(Qr,5)*pow(r,8)*Sr*N*P + 32*Bep*Bepp*Qr*pow(r,4)*pow(Sr,3)*N*P - 32*pow(Bep,2)*pow(Qr,2)*pow(r,6)*pow(Sr,3)*N*P + 16*Al*pow(Bep,2)*pow(Qr,4)*pow(r,8)*pow(Sr,3)*N*P - 32*Bep*Bepp*pow(Qr,3)*pow(r,4)*Sr*(pow(r,2) + 8*Bepp*pow(r,2)*pow(Sr,2))*N*P + 16*Bep*pow(r,4)*Sr*Vp*N*P - 128*pow(Bep,2)*Qr*pow(r,4)*Sr*Vp*N*P + 2*Al*Qr*pow(r,4)*N*pow(P,2) - 32*Al*Bep*pow(Qr,2)*pow(r,4)*N*pow(P,2) + 128*Al*pow(Bep,2)*pow(Qr,3)*pow(r,4)*N*pow(P,2) + (3*Alp*pow(Qr,2)*pow(r,6)*N*pow(P,2))/2. - 24*Alp*Bep*pow(Qr,3)*pow(r,6)*N*pow(P,2) - 34*Bep*pow(r,4)*pow(Sr,2)*N*pow(P,2) + 416*pow(Bep,2)*Qr*pow(r,4)*pow(Sr,2)*N*pow(P,2) + 38*Al*Bep*pow(Qr,2)*pow(r,6)*pow(Sr,2)*N*pow(P,2) - 416*Al*pow(Bep,2)*pow(Qr,3)*pow(r,6)*pow(Sr,2)*N*pow(P,2) - 16*Al*Bep*Bepp*pow(Qr,4)*pow(r,8)*pow(Sr,2)*N*pow(P,2) + 32*Bep*Bepp*pow(r,4)*pow(Sr,4)*N*pow(P,2) - 16*Bep*Bepp*pow(Qr,2)*pow(r,4)*pow(Sr,2)*(pow(r,2) + 16*Bepp*pow(r,2)*pow(Sr,2))*N*pow(P,2) - 64*pow(Bep,2)*pow(r,4)*pow(Sr,2)*Vp*N*pow(P,2) + 2*Al*pow(r,4)*Sr*N*pow(P,3) - 60*Al*Bep*Qr*pow(r,4)*Sr*N*pow(P,3) + 192*Alp*pow(Bep,2)*pow(Qr,3)*pow(r,6)*Sr*N*pow(P,3) - 32*Al*Bep*Bepp*pow(Qr,3)*pow(r,6)*Sr*N*pow(P,3) + 160*pow(Bep,2)*pow(r,4)*pow(Sr,3)*N*pow(P,3) - 8*Bep*pow(Qr,2)*pow(r,4)*Sr*(-48*Al*Bep + 3*Alp*pow(r,2) + 20*Al*Bep*pow(r,2)*pow(Sr,2))*N*pow(P,3) - (3*Alp*pow(r,4)*N*pow(P,4))/4. + 12*Alp*Bep*Qr*pow(r,4)*N*pow(P,4) - 33*Al*Bep*pow(r,4)*pow(Sr,2)*N*pow(P,4) + 400*Al*pow(Bep,2)*Qr*pow(r,4)*pow(Sr,2)*N*pow(P,4) + 8*Bep*pow(Qr,2)*pow(r,4)*(-6*Alp*Bep + 12*Alp*Bep*pow(r,2)*pow(Sr,2) - Al*Bepp*pow(r,2)*pow(Sr,2))*N*pow(P,4) + 12*Alp*Bep*pow(r,4)*Sr*N*pow(P,5) - 96*Alp*pow(Bep,2)*Qr*pow(r,4)*Sr*N*pow(P,5) + 144*Al*pow(Bep,2)*pow(r,4)*pow(Sr,3)*N*pow(P,5) - 48*Alp*pow(Bep,2)*pow(r,4)*pow(Sr,2)*N*pow(P,6) + (64*pow(Bep,2)*pow(r,4)*pow(r_Der_N,2)*pow(Sr,3)*(Qr*pow(r,2)*Sr + P)*(-1 + 8*Bep*Qr + 8*Bep*Sr*P))/N - 12*Alp*Bep*pow(Qr,5)*pow(r,7)*N*(-r + 8*Bep*r*Sr*P) + pow(r_Der_S,2)*(-64*pow(Bep,2)*Qr*pow(r,4)*pow(Sr,2)*N + 512*pow(Bep,3)*pow(Qr,2)*pow(r,4)*pow(Sr,2)*N + 256*pow(Bep,3)*Qr*pow(r,4)*pow(Sr,3)*N*P) - (3*Alp*pow(Qr,4)*pow(r,6)*N*(pow(r,2) - 16*Bep*pow(r,2)*Sr*P - 128*pow(Bep,2)*pow(P,2) + 64*pow(Bep,2)*pow(r,2)*pow(Sr,2)*pow(P,2)))/4. + r_Der_P*(pow(r,5)*Sr*N - 16*Bep*Qr*pow(r,5)*Sr*N - Al*pow(Qr,2)*pow(r,7)*Sr*N + 16*Al*Bep*pow(Qr,3)*pow(r,7)*Sr*N - 64*Al*pow(Bep,2)*pow(Qr,4)*pow(r,7)*Sr*N + 32*pow(Bep,2)*pow(r,3)*pow(Sr,3)*N - 32*pow(Bep,2)*pow(r,5)*pow(Sr,5)*N + 64*pow(Bep,2)*pow(Qr,2)*pow(r,3)*Sr*(pow(r,2) - 4*Bepp*pow(r,2)*pow(Sr,2) + 4*Bepp*pow(r,4)*pow(Sr,4))*N + r_Der_Q*(-256*pow(Bep,3)*pow(r,3)*pow(Sr,3)*N + 256*pow(Bep,3)*pow(r,5)*pow(Sr,5)*N) + 4*Al*Qr*pow(r,5)*N*P - 64*Al*Bep*pow(Qr,2)*pow(r,5)*N*P + 256*Al*pow(Bep,2)*pow(Qr,3)*pow(r,5)*N*P - 16*Bep*pow(r,5)*pow(Sr,2)*N*P + 128*pow(Bep,2)*Qr*pow(r,5)*pow(Sr,2)*N*P + 16*Al*Bep*pow(Qr,2)*pow(r,7)*pow(Sr,2)*N*P - 128*Al*pow(Bep,2)*pow(Qr,3)*pow(r,7)*pow(Sr,2)*N*P + 3*Al*pow(r,5)*Sr*N*pow(P,2) - 112*Al*Bep*Qr*pow(r,5)*Sr*N*pow(P,2) + 704*Al*pow(Bep,2)*pow(Qr,2)*pow(r,5)*Sr*N*pow(P,2) + 64*pow(Bep,2)*pow(r,5)*pow(Sr,3)*N*pow(P,2) - 64*Al*pow(Bep,2)*pow(Qr,2)*pow(r,7)*pow(Sr,3)*N*pow(P,2) - 48*Al*Bep*pow(r,5)*pow(Sr,2)*N*pow(P,3) + 640*Al*pow(Bep,2)*Qr*pow(r,5)*pow(Sr,2)*N*pow(P,3) + 192*Al*pow(Bep,2)*pow(r,5)*pow(Sr,3)*N*pow(P,4)) + 4*Bep*pow(r,4)*pow(Sr,2)*N*V - 64*pow(Bep,2)*Qr*pow(r,4)*pow(Sr,2)*N*V - 64*pow(Bep,2)*pow(r,4)*pow(Sr,3)*N*P*V - 2*Bep*pow(Qr,2)*pow(r,4)*pow(Sr,2)*N*(-3*pow(r,2) + 16*Bepp*pow(r,2)*pow(Sr,2) - 16*Bepp*pow(r,2)*V) + r_Der_Q*(pow(r,4)*N - 16*Bep*Qr*pow(r,4)*N - 3*Al*pow(Qr,2)*pow(r,6)*N + 48*Al*Bep*pow(Qr,3)*pow(r,6)*N - 192*Al*pow(Bep,2)*pow(Qr,4)*pow(r,6)*N + 24*Al*pow(Bep,2)*pow(Qr,4)*pow(r,8)*pow(Sr,2)*N - 32*pow(Bep,2)*pow(r,4)*pow(Sr,4)*N - 16*pow(Bep,2)*pow(Qr,2)*pow(r,4)*(-4 + pow(r,2)*pow(Sr,2))*N - 16*Bep*pow(r,4)*Sr*N*P + 96*pow(Bep,2)*Qr*pow(r,4)*Sr*N*P + 48*Al*Bep*pow(Qr,2)*pow(r,6)*Sr*N*P - 352*Al*pow(Bep,2)*pow(Qr,3)*pow(r,6)*Sr*N*P - 256*pow(Bep,2)*Bepp*Qr*pow(r,4)*pow(Sr,3)*N*P + Al*pow(r,4)*N*pow(P,2) - 16*Al*Bep*Qr*pow(r,4)*N*pow(P,2) + 64*Al*pow(Bep,2)*pow(Qr,2)*pow(r,4)*N*pow(P,2) + 48*pow(Bep,2)*pow(r,4)*pow(Sr,2)*N*pow(P,2) - 208*Al*pow(Bep,2)*pow(Qr,2)*pow(r,6)*pow(Sr,2)*N*pow(P,2) - 256*pow(Bep,2)*Bepp*pow(r,4)*pow(Sr,4)*N*pow(P,2) - 16*Al*Bep*pow(r,4)*Sr*N*pow(P,3) + 96*Al*pow(Bep,2)*Qr*pow(r,4)*Sr*N*pow(P,3) + 56*Al*pow(Bep,2)*pow(r,4)*pow(Sr,2)*N*pow(P,4) + 32*pow(Bep,2)*pow(r,4)*pow(Sr,2)*N*V) + r_Der_N*(Qr*pow(r,5) + 8*Bep*pow(r,3)*pow(Sr,2) - 64*pow(Bep,2)*Qr*pow(r,3)*pow(Sr,2) - 16*Bep*pow(r,5)*pow(Sr,4) + 128*pow(Bep,2)*Qr*pow(r,5)*pow(Sr,4) - 16*Al*pow(Bep,2)*pow(Qr,5)*pow(r,7)*(4 + 3*pow(r,2)*pow(Sr,2)) - 4*Bep*pow(Qr,2)*pow(r,3)*(4*pow(r,2) + 16*Bepp*pow(r,2)*pow(Sr,2) + pow(r,4)*pow(Sr,2)) + pow(r,5)*Sr*P - 40*Bep*Qr*pow(r,5)*Sr*P - 64*pow(Bep,2)*pow(r,3)*pow(Sr,3)*P - 192*Bep*Bepp*Qr*pow(r,5)*pow(Sr,3)*P + 96*pow(Bep,2)*pow(r,5)*pow(Sr,5)*P + pow(Qr,2)*pow(r,3)*Sr*(256*pow(Bep,2)*pow(r,2) + Al*pow(r,4) + 2048*pow(Bep,2)*Bepp*pow(r,2)*pow(Sr,2) + 16*pow(Bep,2)*pow(r,4)*pow(Sr,2))*P - 20*Bep*pow(r,5)*pow(Sr,2)*pow(P,2) - 128*Bep*Bepp*pow(r,5)*pow(Sr,4)*pow(P,2) - 4*Al*Bep*pow(Qr,2)*pow(r,5)*(12 + 5*pow(r,2)*pow(Sr,2))*pow(P,2) + Qr*r*(3*Al*pow(r,4) + 256*pow(Bep,2)*pow(r,4)*pow(Sr,2) + 2304*pow(Bep,2)*Bepp*pow(r,4)*pow(Sr,4))*pow(P,2) - 72*Al*Bep*Qr*pow(r,5)*Sr*pow(P,3) + 16*Al*pow(Bep,2)*pow(Qr,2)*pow(r,5)*Sr*(32 + 5*pow(r,2)*pow(Sr,2))*pow(P,3) + r*Sr*(Al*pow(r,4) + 80*pow(Bep,2)*pow(r,4)*pow(Sr,2) + 768*pow(Bep,2)*Bepp*pow(r,4)*pow(Sr,4))*pow(P,3) - 18*Al*Bep*pow(r,5)*pow(Sr,2)*pow(P,4) + 368*Al*pow(Bep,2)*Qr*pow(r,5)*pow(Sr,2)*pow(P,4) + 72*Al*pow(Bep,2)*pow(r,5)*pow(Sr,3)*pow(P,5) - 2*Al*Bep*pow(Qr,4)*pow(r,6)*(-8*r - 3*pow(r,3)*pow(Sr,2) + 64*Bep*r*Sr*P + 12*Bep*pow(r,3)*pow(Sr,3)*P) + r_Der_Q*(-64*pow(Bep,2)*pow(r,3)*pow(Sr,2) + 512*pow(Bep,3)*Qr*pow(r,3)*pow(Sr,2) + 512*pow(Bep,3)*pow(r,3)*pow(Sr,3)*P) + r_Der_P*(-192*pow(Bep,2)*pow(r,4)*pow(Sr,3) + 1536*pow(Bep,3)*Qr*pow(r,4)*pow(Sr,3) + 128*pow(Bep,2)*pow(r,6)*pow(Sr,5) - 1024*pow(Bep,3)*Qr*pow(r,6)*pow(Sr,5) + 1280*pow(Bep,3)*pow(r,4)*pow(Sr,4)*P - 768*pow(Bep,3)*pow(r,6)*pow(Sr,6)*P) + pow(Qr,3)*pow(r,3)*(64*pow(Bep,2)*pow(r,2) - Al*pow(r,4) + 512*pow(Bep,2)*Bepp*pow(r,2)*pow(Sr,2) + 32*pow(Bep,2)*pow(r,4)*pow(Sr,2) + 8*Al*Bep*pow(r,4)*Sr*P + 192*Al*pow(Bep,2)*pow(r,2)*pow(P,2) + 64*Al*pow(Bep,2)*pow(r,4)*pow(Sr,2)*pow(P,2)) + r_Der_S*(-16*Bep*pow(r,3)*Sr - 128*pow(Bep,2)*Qr*pow(r,3)*Sr*(-2 + pow(r,2)*pow(Sr,2)) + 1024*pow(Bep,3)*pow(Qr,2)*pow(r,3)*Sr*(-1 + pow(r,2)*pow(Sr,2)) + 192*pow(Bep,2)*pow(r,3)*pow(Sr,2)*P + 768*pow(Bep,3)*Qr*pow(r,3)*pow(Sr,2)*(-2 + pow(r,2)*pow(Sr,2))*P - 512*pow(Bep,3)*pow(r,3)*pow(Sr,3)*pow(P,2)) + 8*Bep*pow(r,5)*pow(Sr,2)*V - 64*pow(Bep,2)*Qr*pow(r,5)*pow(Sr,2)*V - 32*pow(Bep,2)*pow(r,5)*pow(Sr,3)*P*V) + r_Der_S*(-4*Bep*pow(Qr,2)*pow(r,6)*Sr*N + 6*Al*Bep*pow(Qr,4)*pow(r,8)*Sr*N - 48*Al*pow(Bep,2)*pow(Qr,5)*pow(r,8)*Sr*N - 16*Bep*pow(r,4)*pow(Sr,3)*N + 160*pow(Bep,2)*Qr*pow(r,4)*pow(Sr,3)*N - 256*pow(Bep,3)*Qr*pow(r,4)*r_Der_Q*pow(Sr,3)*N - 32*pow(Bep,2)*pow(Qr,3)*pow(r,4)*Sr*(-pow(r,2) + 8*Bepp*pow(r,2)*pow(Sr,2))*N + pow(r,4)*N*P - 24*Bep*Qr*pow(r,4)*N*P - 8*Al*Bep*pow(Qr,3)*pow(r,6)*N*P - 64*Bep*Bepp*Qr*pow(r,4)*pow(Sr,2)*N*P - 24*Al*pow(Bep,2)*pow(Qr,4)*pow(r,8)*pow(Sr,2)*N*P + 96*pow(Bep,2)*pow(r,4)*pow(Sr,4)*N*P + pow(Qr,2)*pow(r,2)*(128*pow(Bep,2)*pow(r,2) + Al*pow(r,4) + 512*pow(Bep,2)*Bepp*pow(r,2)*pow(Sr,2) + 16*pow(Bep,2)*pow(r,4)*pow(Sr,2))*N*P - 20*Bep*pow(r,4)*Sr*N*pow(P,2) - 20*Al*Bep*pow(Qr,2)*pow(r,6)*Sr*N*pow(P,2) + 128*Al*pow(Bep,2)*pow(Qr,3)*pow(r,6)*Sr*N*pow(P,2) - 128*Bep*Bepp*pow(r,4)*pow(Sr,3)*N*pow(P,2) + 64*pow(Bep,2)*Qr*pow(r,2)*Sr*(3*pow(r,2) + 20*Bepp*pow(r,2)*pow(Sr,2))*N*pow(P,2) - 24*Al*Bep*Qr*pow(r,4)*N*pow(P,3) + 16*Al*pow(Bep,2)*pow(Qr,2)*pow(r,4)*(8 + 5*pow(r,2)*pow(Sr,2))*N*pow(P,3) + (Al*pow(r,4) + 80*pow(Bep,2)*pow(r,4)*pow(Sr,2) + 768*pow(Bep,2)*Bepp*pow(r,4)*pow(Sr,4))*N*pow(P,3) - 18*Al*Bep*pow(r,4)*Sr*N*pow(P,4) + 176*Al*pow(Bep,2)*Qr*pow(r,4)*Sr*N*pow(P,4) + 72*Al*pow(Bep,2)*pow(r,4)*pow(Sr,2)*N*pow(P,5) + r_Der_P*(-64*pow(Bep,2)*pow(r,3)*pow(Sr,2)*N + 512*pow(Bep,3)*Qr*pow(r,3)*pow(Sr,2)*N + 128*pow(Bep,2)*pow(r,5)*pow(Sr,4)*N - 1024*pow(Bep,3)*Qr*pow(r,5)*pow(Sr,4)*N - 256*pow(Bep,3)*pow(r,3)*pow(Sr,3)*(-1 + 3*pow(r,2)*pow(Sr,2))*N*P) + 8*Bep*pow(r,4)*Sr*N*V - 64*pow(Bep,2)*Qr*pow(r,4)*Sr*N*V - 32*pow(Bep,2)*pow(r,4)*pow(Sr,2)*N*P*V)
   ;
   p_k/=
      pow(r,4) - 16*Bep*Qr*pow(r,4) + 16*Al*Bep*pow(Qr,3)*pow(r,6) - 64*Al*pow(Bep,2)*pow(Qr,4)*pow(r,6) - 32*pow(Bep,2)*pow(r,4)*pow(Sr,4) + 256*pow(Bep,3)*pow(r,4)*r_Der_Q*pow(Sr,4) - pow(Qr,2)*pow(r,2)*(-64*pow(Bep,2)*pow(r,2) + Al*pow(r,4) - 256*pow(Bep,2)*Bepp*pow(r,4)*pow(Sr,4)) - 16*Bep*pow(r,3)*(-r + 8*Bep*Qr*r)*(-1 + Al*pow(Qr,2)*pow(r,2))*Sr*P - pow(r,2)*(-3*Al*pow(r,2) + 48*Al*Bep*Qr*pow(r,2) - 192*Al*pow(Bep,2)*pow(Qr,2)*pow(r,2) - 64*pow(Bep,2)*pow(r,2)*pow(Sr,2) + 64*Al*pow(Bep,2)*pow(Qr,2)*pow(r,4)*pow(Sr,2))*pow(P,2) + 48*Al*Bep*pow(r,3)*(-r + 8*Bep*Qr*r)*Sr*pow(P,3) + 192*Al*pow(Bep,2)*pow(r,4)*pow(Sr,2)*pow(P,4) - (128*pow(Bep,2)*pow(r,5)*r_Der_N*pow(Sr,4)*(-1 + 8*Bep*Qr + 6*Bep*Sr*P))/N + r_Der_S*(128*pow(Bep,2)*pow(r,4)*pow(Sr,3) - 1024*pow(Bep,3)*Qr*pow(r,4)*pow(Sr,3) - 768*pow(Bep,3)*pow(r,4)*pow(Sr,4)*P)
   ;
   return p_k;
}
/*===========================================================================*/
inline double EdGB::compute_S_free_k(
   const double r,
   const double N, const double S,
   const double P, const double Q,
   const double V, const double Vp, 
   const double Al,  const double Alp,
   const double Bep, const double Bepp,
   const double r_Der_N,
   const double r_Der_S,
   const double r_Der_P,
   const double r_Der_Q) const
{
   double Sr= S/r;
   double Qr= Q/r;
   double vn=
      (r*N*pow(P,2))/4. + (7*Al*r*N*pow(P,4))/8. + (3*pow(Al,2)*r*N*pow(P,6))/8. - 6*Bep*N*pow(P,2)*Q - 23*Al*Bep*N*pow(P,4)*Q - 15*pow(Al,2)*Bep*N*pow(P,6)*Q + (r*N*pow(Q,2))/4. + (3*Al*r*N*pow(P,2)*pow(Q,2))/4. + (5*pow(Al,2)*r*N*pow(P,4)*pow(Q,2))/8. - 2*Bep*N*pow(Q,3) + 2*Al*Bep*N*pow(P,2)*pow(Q,3) + 11*pow(Al,2)*Bep*N*pow(P,4)*pow(Q,3) - (5*Al*r*N*pow(Q,4))/8. - (11*pow(Al,2)*r*N*pow(P,2)*pow(Q,4))/8. + 5*Al*Bep*N*pow(Q,5) + 7*pow(Al,2)*Bep*N*pow(P,2)*pow(Q,5) + (3*pow(Al,2)*r*N*pow(Q,6))/8. - 3*pow(Al,2)*Bep*N*pow(Q,7) + (r*N*P*Q)/(2.*S) + (2*Al*r*N*pow(P,3)*Q)/S + (3*pow(Al,2)*r*N*pow(P,5)*Q)/(2.*S) - (4*Bep*N*P*pow(Q,2))/S - (16*Al*Bep*N*pow(P,3)*pow(Q,2))/S - (12*pow(Al,2)*Bep*N*pow(P,5)*pow(Q,2))/S - (Al*r*N*P*pow(Q,3))/S - (2*pow(Al,2)*r*N*pow(P,3)*pow(Q,3))/S + (8*Al*Bep*N*P*pow(Q,4))/S + (16*pow(Al,2)*Bep*N*pow(P,3)*pow(Q,4))/S + (pow(Al,2)*r*N*P*pow(Q,5))/(2.*S) - (4*pow(Al,2)*Bep*N*P*pow(Q,6))/S - 2*Bep*N*pow(P,3)*S - 7*Al*Bep*N*pow(P,5)*S - 3*pow(Al,2)*Bep*N*pow(P,7)*S + (4*Bepp*N*P*Q*S)/r + (12*Al*Bepp*N*pow(P,3)*Q*S)/r - 2*Bep*N*P*pow(Q,2)*S - (32*Bep*Bepp*N*P*pow(Q,2)*S)/pow(r,2) - 6*Al*Bep*N*pow(P,3)*pow(Q,2)*S - (96*Al*Bep*Bepp*N*pow(P,3)*pow(Q,2)*S)/pow(r,2) - 5*pow(Al,2)*Bep*N*pow(P,5)*pow(Q,2)*S - (4*Al*Bepp*N*P*pow(Q,3)*S)/r + 5*Al*Bep*N*P*pow(Q,4)*S + (32*Al*Bep*Bepp*N*P*pow(Q,4)*S)/pow(r,2) + 11*pow(Al,2)*Bep*N*pow(P,3)*pow(Q,4)*S - 3*pow(Al,2)*Bep*N*P*pow(Q,6)*S + (N*pow(S,2))/(2.*r) - (4*Bep*Vp*N*pow(S,2))/r + (3*Al*N*pow(P,2)*pow(S,2))/(2.*r) + (4*Bepp*N*pow(P,2)*pow(S,2))/r - (3*Alp*Bep*N*pow(P,4)*pow(S,2))/r + (12*Al*Bepp*N*pow(P,4)*pow(S,2))/r + (4*Bep*N*Q*pow(S,2))/pow(r,2) + (32*pow(Bep,2)*Vp*N*Q*pow(S,2))/pow(r,2) - (4*Al*Bep*N*pow(P,2)*Q*pow(S,2))/pow(r,2) - (64*Bep*Bepp*N*pow(P,2)*Q*pow(S,2))/pow(r,2) + (24*Alp*pow(Bep,2)*N*pow(P,4)*Q*pow(S,2))/pow(r,2) - (192*Al*Bep*Bepp*N*pow(P,4)*Q*pow(S,2))/pow(r,2) - (64*pow(Bep,2)*N*pow(Q,2)*pow(S,2))/pow(r,3) - (Al*N*pow(Q,2)*pow(S,2))/(2.*r) - (64*Al*pow(Bep,2)*N*pow(P,2)*pow(Q,2)*pow(S,2))/pow(r,3) + (6*Alp*Bep*N*pow(P,2)*pow(Q,2)*pow(S,2))/r - (4*Al*Bepp*N*pow(P,2)*pow(Q,2)*pow(S,2))/r - (4*Al*Bep*N*pow(Q,3)*pow(S,2))/pow(r,2) - (48*Alp*pow(Bep,2)*N*pow(P,2)*pow(Q,3)*pow(S,2))/pow(r,2) + (64*Al*Bep*Bepp*N*pow(P,2)*pow(Q,3)*pow(S,2))/pow(r,2) + (64*Al*pow(Bep,2)*N*pow(Q,4)*pow(S,2))/pow(r,3) - (3*Alp*Bep*N*pow(Q,4)*pow(S,2))/r + (24*Alp*pow(Bep,2)*N*pow(Q,5)*pow(S,2))/pow(r,2) + (4*Bep*N*P*pow(S,3))/pow(r,2) + (32*pow(Bep,2)*Vp*N*P*pow(S,3))/pow(r,2) - (4*Al*Bep*N*pow(P,3)*pow(S,3))/pow(r,2) - (32*Bep*Bepp*N*pow(P,3)*pow(S,3))/pow(r,2) + (24*Alp*pow(Bep,2)*N*pow(P,5)*pow(S,3))/pow(r,2) - (96*Al*Bep*Bepp*N*pow(P,5)*pow(S,3))/pow(r,2) - (128*pow(Bep,2)*N*P*Q*pow(S,3))/pow(r,3) - (128*Al*pow(Bep,2)*N*pow(P,3)*Q*pow(S,3))/pow(r,3) - (4*Al*Bep*N*P*pow(Q,2)*pow(S,3))/pow(r,2) - (48*Alp*pow(Bep,2)*N*pow(P,3)*pow(Q,2)*pow(S,3))/pow(r,2) + (32*Al*Bep*Bepp*N*pow(P,3)*pow(Q,2)*pow(S,3))/pow(r,2) + (128*Al*pow(Bep,2)*N*P*pow(Q,3)*pow(S,3))/pow(r,3) + (24*Alp*pow(Bep,2)*N*P*pow(Q,4)*pow(S,3))/pow(r,2) - (80*pow(Bep,2)*N*pow(P,2)*pow(S,4))/pow(r,3) - (72*Al*pow(Bep,2)*N*pow(P,4)*pow(S,4))/pow(r,3) + (16*pow(Bep,2)*N*pow(Q,2)*pow(S,4))/pow(r,3) + (80*Al*pow(Bep,2)*N*pow(P,2)*pow(Q,2)*pow(S,4))/pow(r,3) - (8*Al*pow(Bep,2)*N*pow(Q,4)*pow(S,4))/pow(r,3) + r_Der_P*((4*Bep*N*S)/r + (12*Al*Bep*N*pow(P,2)*S)/r - (32*pow(Bep,2)*N*Q*S)/pow(r,2) - (96*Al*pow(Bep,2)*N*pow(P,2)*Q*S)/pow(r,2) - (4*Al*Bep*N*pow(Q,2)*S)/r + (32*Al*pow(Bep,2)*N*pow(Q,3)*S)/pow(r,2) - (32*pow(Bep,2)*N*P*pow(S,2))/pow(r,2) - (96*Al*pow(Bep,2)*N*pow(P,3)*pow(S,2))/pow(r,2) + (16*Al*Bep*N*P*Q*pow(S,2))/r - (96*Al*pow(Bep,2)*N*P*pow(Q,2)*pow(S,2))/pow(r,2) - (128*Al*pow(Bep,2)*N*pow(P,2)*Q*pow(S,3))/pow(r,2)) + r_Der_Q*((4*Bep*N*pow(S,2))/r + (4*Al*Bep*N*pow(P,2)*pow(S,2))/r - (32*pow(Bep,2)*N*Q*pow(S,2))/pow(r,2) - (32*Al*pow(Bep,2)*N*pow(P,2)*Q*pow(S,2))/pow(r,2) - (12*Al*Bep*N*pow(Q,2)*pow(S,2))/r + (96*Al*pow(Bep,2)*N*pow(Q,3)*pow(S,2))/pow(r,2) - (32*pow(Bep,2)*N*P*pow(S,3))/pow(r,2) - (32*Al*pow(Bep,2)*N*pow(P,3)*pow(S,3))/pow(r,2) + (96*Al*pow(Bep,2)*N*P*pow(Q,2)*pow(S,3))/pow(r,2)) + pow(r_Der_S,2)*((128*pow(Bep,2)*N*pow(S,4))/pow(r,3) - (768*pow(Bep,3)*N*Q*pow(S,4))/pow(r,4) - (768*pow(Bep,3)*N*P*pow(S,5))/pow(r,4)) + pow(r_Der_N,2)*((-256*pow(Bep,3)*P*pow(S,5))/(pow(r,4)*N) - (256*pow(Bep,3)*Q*pow(S,6))/(pow(r,4)*N)) - (r*N*V)/2. - (3*Al*r*N*pow(P,2)*V)/2. + 4*Bep*N*Q*V + 12*Al*Bep*N*pow(P,2)*Q*V + (Al*r*N*pow(Q,2)*V)/2. - 4*Al*Bep*N*pow(Q,3)*V + 4*Bep*N*P*S*V + 12*Al*Bep*N*pow(P,3)*S*V - 4*Al*Bep*N*P*pow(Q,2)*S*V + (32*pow(Bep,2)*N*pow(S,4)*V)/pow(r,3) + r_Der_S*(N*S + 3*Al*N*pow(P,2)*S - (12*Bep*N*Q*S)/r - (36*Al*Bep*N*pow(P,2)*Q*S)/r - Al*N*pow(Q,2)*S + (32*pow(Bep,2)*N*pow(Q,2)*S)/pow(r,2) + (96*Al*pow(Bep,2)*N*pow(P,2)*pow(Q,2)*S)/pow(r,2) + (12*Al*Bep*N*pow(Q,3)*S)/r - (32*Al*pow(Bep,2)*N*pow(Q,4)*S)/pow(r,2) - (12*Bep*N*P*pow(S,2))/r - (44*Al*Bep*N*pow(P,3)*pow(S,2))/r + (96*pow(Bep,2)*N*P*Q*pow(S,2))/pow(r,2) + (288*Al*pow(Bep,2)*N*pow(P,3)*Q*pow(S,2))/pow(r,2) + (20*Al*Bep*N*P*pow(Q,2)*pow(S,2))/r - (160*Al*pow(Bep,2)*N*P*pow(Q,3)*pow(S,2))/pow(r,2) + (48*pow(Bep,2)*N*pow(P,2)*pow(S,3))/pow(r,2) + (168*Al*pow(Bep,2)*N*pow(P,4)*pow(S,3))/pow(r,2) + (16*pow(Bep,2)*N*pow(Q,2)*pow(S,3))/pow(r,2) - (80*Al*pow(Bep,2)*N*pow(P,2)*pow(Q,2)*pow(S,3))/pow(r,2) - (24*Al*pow(Bep,2)*N*pow(Q,4)*pow(S,3))/pow(r,2) + (256*pow(Bep,3)*r_Der_P*N*pow(S,4))/pow(r,4) + (256*pow(Bep,2)*Bepp*N*P*Q*pow(S,4))/pow(r,4) - (32*pow(Bep,2)*N*pow(S,5))/pow(r,4) + (256*pow(Bep,3)*r_Der_Q*N*pow(S,5))/pow(r,4) + (256*pow(Bep,2)*Bepp*N*pow(Q,2)*pow(S,5))/pow(r,4) - (32*pow(Bep,2)*N*pow(S,3)*V)/pow(r,2)) + r_Der_N*((4*Bep*Q*pow(S,2))/r + (12*Al*Bep*pow(P,2)*Q*pow(S,2))/r - (32*pow(Bep,2)*pow(Q,2)*pow(S,2))/pow(r,2) - (96*Al*pow(Bep,2)*pow(P,2)*pow(Q,2)*pow(S,2))/pow(r,2) - (4*Al*Bep*pow(Q,3)*pow(S,2))/r + (32*Al*pow(Bep,2)*pow(Q,4)*pow(S,2))/pow(r,2) + (4*Bep*P*pow(S,3))/r + (4*Al*Bep*pow(P,3)*pow(S,3))/r - (32*pow(Bep,2)*P*Q*pow(S,3))/pow(r,2) - (96*Al*pow(Bep,2)*pow(P,3)*Q*pow(S,3))/pow(r,2) + (4*Al*Bep*P*pow(Q,2)*pow(S,3))/r - (32*Al*pow(Bep,2)*P*pow(Q,3)*pow(S,3))/pow(r,2) + (32*pow(Bep,2)*pow(S,4))/pow(r,4) - (256*pow(Bep,3)*r_Der_Q*pow(S,4))/pow(r,4) - (16*pow(Bep,2)*pow(P,2)*pow(S,4))/pow(r,2) - (24*Al*pow(Bep,2)*pow(P,4)*pow(S,4))/pow(r,2) - (256*pow(Bep,2)*Bepp*pow(Q,2)*pow(S,4))/pow(r,4) + (16*pow(Bep,2)*pow(Q,2)*pow(S,4))/pow(r,2) - (16*Al*pow(Bep,2)*pow(P,2)*pow(Q,2)*pow(S,4))/pow(r,2) - (24*Al*pow(Bep,2)*pow(Q,4)*pow(S,4))/pow(r,2) - (256*pow(Bep,3)*r_Der_P*pow(S,5))/pow(r,4) - (256*pow(Bep,2)*Bepp*P*Q*pow(S,5))/pow(r,4) + r_Der_S*((-64*pow(Bep,2)*pow(S,3))/pow(r,3) + (512*pow(Bep,3)*Q*pow(S,3))/pow(r,4) + (256*pow(Bep,3)*P*pow(S,4))/pow(r,4) + (128*pow(Bep,2)*pow(S,5))/pow(r,3) - (1024*pow(Bep,3)*Q*pow(S,5))/pow(r,4) - (768*pow(Bep,3)*P*pow(S,6))/pow(r,4)) - (32*pow(Bep,2)*pow(S,4)*V)/pow(r,2))
   ;
   vn/=
      1 - 16*Bep*Qr + 64*pow(Bep,2)*pow(Qr,2) - Al*pow(Qr,2)*pow(r,2) + 16*Al*Bep*pow(Qr,3)*pow(r,2) - 64*Al*pow(Bep,2)*pow(Qr,4)*pow(r,2) - 32*pow(Bep,2)*pow(Sr,4) + 256*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,2)*pow(Sr,4) + 256*pow(Bep,3)*r_Der_Q*pow(Sr,4) - 16*Bep*Sr*P + 128*pow(Bep,2)*Qr*Sr*P + 16*Al*Bep*pow(Qr,2)*pow(r,2)*Sr*P - 128*Al*pow(Bep,2)*pow(Qr,3)*pow(r,2)*Sr*P + 3*Al*pow(P,2) - 48*Al*Bep*Qr*pow(P,2) + 192*Al*pow(Bep,2)*pow(Qr,2)*pow(P,2) + 64*pow(Bep,2)*pow(Sr,2)*pow(P,2) - 64*Al*pow(Bep,2)*pow(Qr,2)*pow(r,2)*pow(Sr,2)*pow(P,2) - 48*Al*Bep*Sr*pow(P,3) + 384*Al*pow(Bep,2)*Qr*Sr*pow(P,3) + 192*Al*pow(Bep,2)*pow(Sr,2)*pow(P,4) - (128*pow(Bep,2)*r_Der_N*pow(Sr,4)*(-r + 8*Bep*Qr*r + 6*Bep*r*Sr*P))/N - 128*r_Der_S*(-(pow(Bep,2)*pow(Sr,3)) + 8*pow(Bep,3)*Qr*pow(Sr,3) + 6*pow(Bep,3)*pow(Sr,4)*P)
   ;
   return vn;
}
/*===========================================================================*/
inline double EdGB::compute_res_N(
   const double r,
   const double N, const double S,
   const double P, const double Q,
   const double V,  
   const double Al,
   const double Bep, const double Bepp,
   const double r_Der_N,
   const double r_Der_P,
   const double r_Der_Q) const
{
   return
      -5*Bep*N*pow(P,2)*Q - (9*Al*Bep*N*pow(P,4)*Q)/2. + Bep*N*pow(Q,3) + 5*Al*Bep*N*pow(P,2)*pow(Q,3) - (Al*Bep*N*pow(Q,5))/2. + (r*N*P*Q)/(2.*S) + (Al*r*N*pow(P,3)*Q)/(2.*S) - (4*Bep*N*P*pow(Q,2))/S - (4*Al*Bep*N*pow(P,3)*pow(Q,2))/S - (Al*r*N*P*pow(Q,3))/(2.*S) + (4*Al*Bep*N*P*pow(Q,4))/S + (4*Bepp*N*P*Q*S)/r - (32*Bep*Bepp*N*P*pow(Q,2)*S)/pow(r,2) - (2*Bep*N*Q*pow(S,2))/pow(r,2) + (16*pow(Bep,2)*r_Der_Q*N*Q*pow(S,2))/pow(r,2) - (48*Bep*Bepp*N*pow(P,2)*Q*pow(S,2))/pow(r,2) + (16*Bep*Bepp*N*pow(Q,3)*pow(S,2))/pow(r,2) + r_Der_P*((4*Bep*N*S)/r - (32*pow(Bep,2)*N*Q*S)/pow(r,2) - (48*pow(Bep,2)*N*P*pow(S,2))/pow(r,2)) + r_Der_N*(1 - (16*Bep*Q)/r + (64*pow(Bep,2)*pow(Q,2))/pow(r,2) - (20*Bep*P*S)/r + (160*pow(Bep,2)*P*Q*S)/pow(r,2) + (96*pow(Bep,2)*pow(P,2)*pow(S,2))/pow(r,2)) + 2*Bep*N*Q*V
   ;
}
/*===========================================================================*/
inline double EdGB::compute_res_S(
   const double r,
   const double S,
   const double P,  const double Q,
   const double V,
   const double Al,
   const double Bep, const double Bepp,
   const double r_Der_S,
   const double r_Der_P,
   const double r_Der_Q) const
{
   return
      4*Bep*pow(P,3) + 6*Al*Bep*pow(P,5) - r*P*Q - Al*r*pow(P,3)*Q + 12*Bep*P*pow(Q,2) + 4*Al*Bep*pow(P,3)*pow(Q,2) + Al*r*P*pow(Q,3) - 10*Al*Bep*P*pow(Q,4) - (r*pow(P,2))/(2.*S) - (3*Al*r*pow(P,4))/(4.*S) + (4*Bep*pow(P,2)*Q)/S + (6*Al*Bep*pow(P,4)*Q)/S - (r*pow(Q,2))/(2.*S) + (Al*r*pow(P,2)*pow(Q,2))/(2.*S) + (4*Bep*pow(Q,3))/S - (4*Al*Bep*pow(P,2)*pow(Q,3))/S + (Al*r*pow(Q,4))/(4.*S) - (2*Al*Bep*pow(Q,5))/S + S/r - (8*Bep*Q*S)/pow(r,2) + 10*Bep*pow(P,2)*Q*S + 9*Al*Bep*pow(P,4)*Q*S - (8*Bepp*pow(Q,2)*S)/r - 2*Bep*pow(Q,3)*S + (64*Bep*Bepp*pow(Q,3)*S)/pow(r,2) - 10*Al*Bep*pow(P,2)*pow(Q,3)*S + Al*Bep*pow(Q,5)*S - (8*Bep*P*pow(S,2))/pow(r,2) - (8*Bepp*P*Q*pow(S,2))/r + (128*Bep*Bepp*P*pow(Q,2)*pow(S,2))/pow(r,2) + (4*Bep*Q*pow(S,3))/pow(r,2) + (96*Bep*Bepp*pow(P,2)*Q*pow(S,3))/pow(r,2) - (32*Bep*Bepp*pow(Q,3)*pow(S,3))/pow(r,2) + r_Der_S*(2 - (32*Bep*Q)/r + (128*pow(Bep,2)*pow(Q,2))/pow(r,2) - (40*Bep*P*S)/r + (320*pow(Bep,2)*P*Q*S)/pow(r,2) + (192*pow(Bep,2)*pow(P,2)*pow(S,2))/pow(r,2)) + r_Der_P*((-8*Bep*pow(S,2))/r + (64*pow(Bep,2)*Q*pow(S,2))/pow(r,2) + (96*pow(Bep,2)*P*pow(S,3))/pow(r,2)) + r_Der_Q*((-8*Bep*S)/r + (64*pow(Bep,2)*Q*S)/pow(r,2) + (64*pow(Bep,2)*P*pow(S,2))/pow(r,2) - (32*pow(Bep,2)*Q*pow(S,3))/pow(r,2)) + 8*Bep*P*V - (r*V)/S + (8*Bep*Q*V)/S - 4*Bep*Q*S*V
   ;
}
/*===========================================================================*/
#endif
