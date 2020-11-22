#include <cmath>
using std::pow;
using std::exp;
using std::fabs;
#include <iomanip>
using std::setprecision;
using std::setw;
#include <iostream>
using std::cout; 
using std::endl;
#include <vector>
using std::vector;
#include <string>
using std::string;

#include "edgb.hpp"
#include "radial_pts.hpp"
#include "field.hpp"
#include "fd_stencils.hpp"

/*===========================================================================*/
/*===========================================================================*/
EdGB::EdGB(
	const double dt, const double dx, const double cl, const int nx,
	const double mu, const double la,
	const double gbc1, const double gbc2,
	const vector<double> &rvec)
: dt{dt},
  dx{dx},
  nx{nx},

  mu{mu},
  la{la},
  
  gbc1{gbc1},
  gbc2{gbc2},
 
  cl{cl},
 
  r_v(nx),

  S_free_k{0},
  S_free_k1{0},
  S_free_k2{0},
  S_free_k3{0},
  S_free_k4{0},

  f_k(nx,0),
  p_k(nx,0),
  q_k(nx,0),

  f_k1(nx,0),
  p_k1(nx,0),
  q_k1(nx,0),

  f_k2(nx,0),
  p_k2(nx,0),
  q_k2(nx,0),

  f_k3(nx,0),
  p_k3(nx,0),
  q_k3(nx,0),

  f_k4(nx,0),
  p_k4(nx,0),
  q_k4(nx,0),

  ingoing_c{0},
  outgoing_c{0},

  V_v(nx,0),
  Vp_v(nx,0),
  
  Al_v(nx,0),
  Alp_v(nx,0),

  Be_v(nx,0),
  Bep_v(nx,0),
  Bepp_v(nx,0)
{
   for (int i=0; i<nx; ++i) {
      r_v[i]=  rvec[i];
   }
}
/*===========================================================================*/
EdGB::~EdGB(void)
{
}
/*===========================================================================*/
void EdGB::rescale(vector<double> &vec)
{
   double vn= vec[nx-1];
   for (int i=0; i<nx; ++i) {
      vec[i]/= vn;
   }
   return;
}
/*===========================================================================*/
void EdGB::compute_potentials(const vector<double> &f_v)
{
   for (int i=0; i<nx; ++i) {
      V_v[i]= 
         (1./2.)*pow(mu,2)*pow(f_v[i],2)
      +  la*pow(f_v[i],4)
      ;
      Vp_v[i]= 
         pow(mu,2)*f_v[i]
      +  (4.)*la*pow(f_v[i],3)
      ;
      Al_v[i]=
         0.0
      ;
      Alp_v[i]=
         0.0
      ;
      Be_v[i]=
         gbc1*f_v[i]
      +  (1./8.)*gbc2*pow(f_v[i],2)
      ;
      Bep_v[i]=
         gbc1
      +  (1./4.)*gbc2*f_v[i]
      ;
      Bepp_v[i]=
         (1./4.)*gbc2
      ;
   }
   return;
}
/*===========================================================================*/
double EdGB::compute_res_N(
   const double r,
   const double N, const double S,
   const double P, const double Q,
   const double V,  
   const double Al,
   const double Bep, const double Bepp,
   const double r_Der_N,
   const double r_Der_P,
   const double r_Der_Q)
{
   return
      -5*Bep*N*pow(P,2)*Q - (9*Al*Bep*N*pow(P,4)*Q)/2. + Bep*N*pow(Q,3) + 5*Al*Bep*N*pow(P,2)*pow(Q,3) - (Al*Bep*N*pow(Q,5))/2. + (r*N*P*Q)/(2.*S) + (Al*r*N*pow(P,3)*Q)/(2.*S) - (4*Bep*N*P*pow(Q,2))/S - (4*Al*Bep*N*pow(P,3)*pow(Q,2))/S - (Al*r*N*P*pow(Q,3))/(2.*S) + (4*Al*Bep*N*P*pow(Q,4))/S + (4*Bepp*N*P*Q*S)/r - (32*Bep*Bepp*N*P*pow(Q,2)*S)/pow(r,2) - (2*Bep*N*Q*pow(S,2))/pow(r,2) + (16*pow(Bep,2)*r_Der_Q*N*Q*pow(S,2))/pow(r,2) - (48*Bep*Bepp*N*pow(P,2)*Q*pow(S,2))/pow(r,2) + (16*Bep*Bepp*N*pow(Q,3)*pow(S,2))/pow(r,2) + r_Der_P*((4*Bep*N*S)/r - (32*pow(Bep,2)*N*Q*S)/pow(r,2) - (48*pow(Bep,2)*N*P*pow(S,2))/pow(r,2)) + r_Der_N*(1 - (16*Bep*Q)/r + (64*pow(Bep,2)*pow(Q,2))/pow(r,2) - (20*Bep*P*S)/r + (160*pow(Bep,2)*P*Q*S)/pow(r,2) + (96*pow(Bep,2)*pow(P,2)*pow(S,2))/pow(r,2)) + 2*Bep*N*Q*V
   ;
}
/*===========================================================================*/
double EdGB::compute_res_S(
   const double r,
   const double S,
   const double P,  const double Q,
   const double V,
   const double Al,
   const double Bep, const double Bepp,
   const double r_Der_S,
   const double r_Der_P,
   const double r_Der_Q)
{
   return
      4*Bep*pow(P,3) + 6*Al*Bep*pow(P,5) - r*P*Q - Al*r*pow(P,3)*Q + 12*Bep*P*pow(Q,2) + 4*Al*Bep*pow(P,3)*pow(Q,2) + Al*r*P*pow(Q,3) - 10*Al*Bep*P*pow(Q,4) - (r*pow(P,2))/(2.*S) - (3*Al*r*pow(P,4))/(4.*S) + (4*Bep*pow(P,2)*Q)/S + (6*Al*Bep*pow(P,4)*Q)/S - (r*pow(Q,2))/(2.*S) + (Al*r*pow(P,2)*pow(Q,2))/(2.*S) + (4*Bep*pow(Q,3))/S - (4*Al*Bep*pow(P,2)*pow(Q,3))/S + (Al*r*pow(Q,4))/(4.*S) - (2*Al*Bep*pow(Q,5))/S + S/r - (8*Bep*Q*S)/pow(r,2) + 10*Bep*pow(P,2)*Q*S + 9*Al*Bep*pow(P,4)*Q*S - (8*Bepp*pow(Q,2)*S)/r - 2*Bep*pow(Q,3)*S + (64*Bep*Bepp*pow(Q,3)*S)/pow(r,2) - 10*Al*Bep*pow(P,2)*pow(Q,3)*S + Al*Bep*pow(Q,5)*S - (8*Bep*P*pow(S,2))/pow(r,2) - (8*Bepp*P*Q*pow(S,2))/r + (128*Bep*Bepp*P*pow(Q,2)*pow(S,2))/pow(r,2) + (4*Bep*Q*pow(S,3))/pow(r,2) + (96*Bep*Bepp*pow(P,2)*Q*pow(S,3))/pow(r,2) - (32*Bep*Bepp*pow(Q,3)*pow(S,3))/pow(r,2) + r_Der_S*(2 - (32*Bep*Q)/r + (128*pow(Bep,2)*pow(Q,2))/pow(r,2) - (40*Bep*P*S)/r + (320*pow(Bep,2)*P*Q*S)/pow(r,2) + (192*pow(Bep,2)*pow(P,2)*pow(S,2))/pow(r,2)) + r_Der_P*((-8*Bep*pow(S,2))/r + (64*pow(Bep,2)*Q*pow(S,2))/pow(r,2) + (96*pow(Bep,2)*P*pow(S,3))/pow(r,2)) + r_Der_Q*((-8*Bep*S)/r + (64*pow(Bep,2)*Q*S)/pow(r,2) + (64*pow(Bep,2)*P*pow(S,2))/pow(r,2) - (32*pow(Bep,2)*Q*pow(S,3))/pow(r,2)) + 8*Bep*P*V - (r*V)/S + (8*Bep*Q*V)/S - 4*Bep*Q*S*V
   ;
}
/*===========================================================================*/
/* solve for metric using relaxation */
/*===========================================================================*/
void EdGB::solve_for_metric_relaxation(
   const int exc_i,
   const vector<double> &f_v,
   const vector<double> &p_v,
   const vector<double> &q_v,
   vector<double> &n_v,
   vector<double> &s_v)
{
   compute_potentials(f_v);

   const double err_tolerance= 1e-12;
   double err= 0;
   do {
      for (int i=exc_i; i<nx-2; ++i) {
         double x=  dx*(2*i+1)/2;

         double cf= 1-x/cl;
         double r=  x/cf;
         double dr= dx*pow(cf,-2);

         double P= (p_v[i+1] + p_v[i])/2;
         double Q= (q_v[i+1] + q_v[i])/2;
         double N= (n_v[i+1] + n_v[i])/2;
         double S= (s_v[i+1] + s_v[i])/2;

         double V=    (   V_v[i+1]+   V_v[i])/2;
         double Al=   (  Al_v[i+1]+  Al_v[i])/2;
         double Bep=  ( Bep_v[i+1]+ Bep_v[i])/2;
         double Bepp= (Bepp_v[i+1]+Bepp_v[i])/2;

         double r_Der_P= (p_v[i+1] - p_v[i])/dr;
         double r_Der_Q= (q_v[i+1] - q_v[i])/dr;
         double r_Der_N= (n_v[i+1] - n_v[i])/dr;
         double r_Der_S= (s_v[i+1] - s_v[i])/dr;
/*---------------------------------------------------------------------------*/
         double res_N= compute_res_N(
            r,
            Al, S, P,  Q,
            V,  
            Al,
            Bep, Bepp,
            r_Der_N, r_Der_P, r_Der_Q
         );
         double jac_N= 
            (-5*Bep*pow(P,2)*Q)/2. - (9*Al*Bep*pow(P,4)*Q)/4. + (Bep*pow(Q,3))/2. + (5*Al*Bep*pow(P,2)*pow(Q,3))/2. - (Al*Bep*pow(Q,5))/4. + (r*P*Q)/(4.*S) + (Al*r*pow(P,3)*Q)/(4.*S) - (2*Bep*P*pow(Q,2))/S - (2*Al*Bep*pow(P,3)*pow(Q,2))/S - (Al*r*P*pow(Q,3))/(4.*S) + (2*Al*Bep*P*pow(Q,4))/S + (2*Bepp*P*Q*S)/r - (16*Bep*Bepp*P*pow(Q,2)*S)/pow(r,2) - (Bep*Q*pow(S,2))/pow(r,2) + (8*pow(Bep,2)*r_Der_Q*Q*pow(S,2))/pow(r,2) - (24*Bep*Bepp*pow(P,2)*Q*pow(S,2))/pow(r,2) + (8*Bep*Bepp*pow(Q,3)*pow(S,2))/pow(r,2) + r_Der_P*((2*Bep*S)/r - (16*pow(Bep,2)*Q*S)/pow(r,2) - (24*pow(Bep,2)*P*pow(S,2))/pow(r,2)) + (1 - (16*Bep*Q)/r + (64*pow(Bep,2)*pow(Q,2))/pow(r,2) - (20*Bep*P*S)/r + (160*pow(Bep,2)*P*Q*S)/pow(r,2) + (96*pow(Bep,2)*pow(P,2)*pow(S,2))/pow(r,2))/dr + Bep*Q*V
	 ;
/*---------------------------------------------------------------------------*/
	 double res_S= compute_res_S(
            r,
            S, P,  Q,
            V,
            Al,
            Bep, Bepp,
            r_Der_S, r_Der_P, r_Der_Q
         );
         double jac_S= 
            1/(2.*r) - (4*Bep*Q)/pow(r,2) + 5*Bep*pow(P,2)*Q + (9*Al*Bep*pow(P,4)*Q)/2. - (4*Bepp*pow(Q,2))/r - Bep*pow(Q,3) + (32*Bep*Bepp*pow(Q,3))/pow(r,2) - 5*Al*Bep*pow(P,2)*pow(Q,3) + (Al*Bep*pow(Q,5))/2. + (r*pow(P,2))/(4.*pow(S,2)) + (3*Al*r*pow(P,4))/(8.*pow(S,2)) - (2*Bep*pow(P,2)*Q)/pow(S,2) - (3*Al*Bep*pow(P,4)*Q)/pow(S,2) + (r*pow(Q,2))/(4.*pow(S,2)) - (Al*r*pow(P,2)*pow(Q,2))/(4.*pow(S,2)) - (2*Bep*pow(Q,3))/pow(S,2) + (2*Al*Bep*pow(P,2)*pow(Q,3))/pow(S,2) - (Al*r*pow(Q,4))/(8.*pow(S,2)) + (Al*Bep*pow(Q,5))/pow(S,2) - (8*Bep*P*S)/pow(r,2) - (8*Bepp*P*Q*S)/r + (128*Bep*Bepp*P*pow(Q,2)*S)/pow(r,2) + (6*Bep*Q*pow(S,2))/pow(r,2) + (144*Bep*Bepp*pow(P,2)*Q*pow(S,2))/pow(r,2) - (48*Bep*Bepp*pow(Q,3)*pow(S,2))/pow(r,2) + r_Der_S*((-20*Bep*P)/r + (160*pow(Bep,2)*P*Q)/pow(r,2) + (192*pow(Bep,2)*pow(P,2)*S)/pow(r,2)) + r_Der_P*((-8*Bep*S)/r + (64*pow(Bep,2)*Q*S)/pow(r,2) + (144*pow(Bep,2)*P*pow(S,2))/pow(r,2)) + (2 - (32*Bep*Q)/r + (128*pow(Bep,2)*pow(Q,2))/pow(r,2) - (40*Bep*P*S)/r + (320*pow(Bep,2)*P*Q*S)/pow(r,2) + (192*pow(Bep,2)*pow(P,2)*pow(S,2))/pow(r,2))/dr + r_Der_Q*((-4*Bep)/r + (32*pow(Bep,2)*Q)/pow(r,2) + (64*pow(Bep,2)*P*S)/pow(r,2) - (48*pow(Bep,2)*Q*pow(S,2))/pow(r,2)) - 2*Bep*Q*V + (r*V)/(2.*pow(S,2)) - (4*Bep*Q*V)/pow(S,2)
         ;
/*---------------------------------------------------------------------------*/
         n_v[i+1]-= res_N/jac_N;
         s_v[i+1]-= res_S/jac_S;
         err+= fabs(res_N);
         err+= fabs(res_S);
      }
      n_v[nx-1]= n_v[nx-2];
      s_v[nx-1]= 0;
      rescale(n_v);
      err/= nx;
   } while (err>err_tolerance);
   return;
}
/*===========================================================================*/
/* public interface to solve for metric fields (for initin data) */
/*===========================================================================*/
void EdGB::solve_metric_fields(
   const int exc_i,
   const Field &phi_f, const Field &phi_p, const Field &phi_q,
   Field &n, Field &s) 
{
   solve_for_metric_relaxation(exc_i,
      phi_f.n,phi_p.n,phi_q.n,
      n.n,s.n
   );
   for (int i=0; i<nx; ++i) {
      n.inter_2[i]= n.n[i];
      n.inter_3[i]= n.n[i];
      n.inter_4[i]= n.n[i];
      n.np1[i]=     n.n[i];

      s.inter_2[i]= s.n[i];
      s.inter_3[i]= s.n[i];
      s.inter_4[i]= s.n[i];
      s.np1[i]=     s.n[i];
   }	
   return;
}
/*===========================================================================*/
double EdGB::compute_S_free_k(
   const double r,
   const double N, const double S,
   const double P, const double Q,
   const double V, const double Vp, 
   const double Al,  const double Alp,
   const double Bep, const double Bepp,
   const double r_Der_N,
   const double r_Der_S,
   const double r_Der_P,
   const double r_Der_Q)
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
double EdGB::compute_p_k(
   const double r,
   const double N, const double S,
   const double P, const double Q,
   const double V, const double Vp, 
   const double Al,  const double Alp,
   const double Bep, const double Bepp,
   const double r_Der_N,
   const double r_Der_S,
   const double r_Der_P,
   const double r_Der_Q)
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
void EdGB::compute_fpq_ki(
   const int vn,
   const int exc_i,
   const vector<double> &n_v,
   const vector<double> &s_v,
   const vector<double> &f_v, 
   const vector<double> &p_v, 
   const vector<double> &q_v)
{
   compute_potentials(f_v);

   double 
      r, cf, 
      n, s, p, q, 
      r_Der_n, r_Der_s, r_Der_p, r_Der_q,
      V, Vp, Al, Alp, Bep, Bepp
   ;
/*---------------------------------------------------------------------------*/
   for (int i=exc_i+2; i<nx-2; ++i) {
      r= r_v[i];
      cf= 1/(1+r/cl);

      n= n_v[i];
      s= s_v[i];
      p= p_v[i];
      q= q_v[i];

      V=  V_v[i];
      Vp= Vp_v[i];
      Al=   Al_v[i];
      Alp=  Alp_v[i];
      Bep=  Bep_v[i];
      Bepp= Bepp_v[i];

      r_Der_n= pow(cf,2)*Dx_ptc_4th(n_v[i+2],n_v[i+1],n_v[i-1],n_v[i-2],dx);
      r_Der_s= pow(cf,2)*Dx_ptc_4th(s_v[i+2],s_v[i+1],s_v[i-1],s_v[i-2],dx);
      r_Der_p= pow(cf,2)*Dx_ptc_4th(p_v[i+2],p_v[i+1],p_v[i-1],p_v[i-2],dx);
      r_Der_q= pow(cf,2)*Dx_ptc_4th(q_v[i+2],q_v[i+1],q_v[i-1],q_v[i-2],dx);

      f_k[i]= 
         n*(p+s*q)
      ;
      p_k[i]= compute_p_k(
         r,
         n, s,
         p, q,
         V, Vp, 
         Al,  Alp,
         Bep, Bepp,
         r_Der_n,
         r_Der_s,
         r_Der_p,
         r_Der_q
      );
      q_k[i]= pow(cf,2)*Dx_ptc_4th(
         n_v[i+2]*(p_v[i+2]+s_v[i+2]*q_v[i+2]),
         n_v[i+1]*(p_v[i+1]+s_v[i+1]*q_v[i+1]),
         n_v[i-1]*(p_v[i-1]+s_v[i-1]*q_v[i-1]),
         n_v[i-2]*(p_v[i-2]+s_v[i-2]*q_v[i-2]),
         dx
      );
   }
/*---------------------------------------------------------------------------*/
   int i= exc_i+1;
   r= r_v[i];
   cf= 1/(1+r/cl);

   n= n_v[i];
   s= s_v[i];
   p= p_v[i];
   q= q_v[i];

   V=    V_v[i];
   Vp=  Vp_v[i];
   Al=   Al_v[i];
   Alp= Alp_v[i];
   Bep=   Bep_v[i];
   Bepp= Bepp_v[i];

   r_Der_n= pow(cf,2)*Dx_ptp1_4th(n_v[i+3],n_v[i+2],n_v[i+1],n_v[i],n_v[i-1],dx);
   r_Der_s= pow(cf,2)*Dx_ptp1_4th(s_v[i+3],s_v[i+2],s_v[i+1],s_v[i],s_v[i-1],dx);
   r_Der_p= pow(cf,2)*Dx_ptp1_4th(p_v[i+3],p_v[i+2],p_v[i+1],p_v[i],p_v[i-1],dx);
   r_Der_q= pow(cf,2)*Dx_ptp1_4th(q_v[i+3],q_v[i+2],q_v[i+1],q_v[i],q_v[i-1],dx);

   f_k[i]=
      n*(p+s*q)
   ;
   p_k[i]= compute_p_k(
      r,
      n, s,
      p,  q,
      V,   Vp, 
      Al,  Alp,
      Bep, Bepp,
      r_Der_n,  r_Der_s,
      r_Der_p,  r_Der_q
   );
   q_k[i]= pow(cf,2)*Dx_ptp1_4th(
      n_v[i+3]*(p_v[i+3]+s_v[i+3]*q_v[i+3]),
      n_v[i+2]*(p_v[i+2]+s_v[i+2]*q_v[i+2]),
      n_v[i+1]*(p_v[i+1]+s_v[i+1]*q_v[i+1]),
      n_v[i  ]*(p_v[i  ]+s_v[i  ]*q_v[i  ]),
      n_v[i-1]*(p_v[i-1]+s_v[i-1]*q_v[i-1]),
      dx
   );
/*---------------------------------------------------------------------------*/
/* if excising then freely evolve p, q, and s */
/*---------------------------------------------------------------------------*/
   if (exc_i>0) {
      i= exc_i;
      r= r_v[i];
      cf= 1/(1+r/cl);

      n= n_v[i];
      s= s_v[i];
      p= p_v[i];
      q= q_v[i];

      V=    V_v[i];
      Vp=  Vp_v[i];
      Al=   Al_v[i];
      Alp= Alp_v[i];
      Bep=   Bep_v[i];
      Bepp= Bepp_v[i];

      r_Der_n= pow(cf,2)*Dx_ptp0_4th(n_v[i+4],n_v[i+3],n_v[i+2],n_v[i+1],n_v[i],dx);
      r_Der_s= pow(cf,2)*Dx_ptp0_4th(s_v[i+4],s_v[i+3],s_v[i+2],s_v[i+1],s_v[i],dx);
      r_Der_p= pow(cf,2)*Dx_ptp0_4th(p_v[i+4],p_v[i+3],p_v[i+2],p_v[i+1],p_v[i],dx);
      r_Der_q= pow(cf,2)*Dx_ptp0_4th(q_v[i+4],q_v[i+3],q_v[i+2],q_v[i+1],q_v[i],dx);

      f_k[i]=
         n*(p+s*q)
      ;
      p_k[i]= compute_p_k(
         r,
         n, s,
         p,  q,
         V,   Vp, 
         Al,  Alp,
         Bep, Bepp,
         r_Der_n,  r_Der_s,
         r_Der_p,  r_Der_q
      );
      q_k[i]= pow(cf,2)*Dx_ptp0_4th(
         n_v[i+4]*(p_v[i+4]+s_v[i+4]*q_v[i+4]),
         n_v[i+3]*(p_v[i+3]+s_v[i+3]*q_v[i+3]),
         n_v[i+2]*(p_v[i+2]+s_v[i+2]*q_v[i+2]),
         n_v[i+1]*(p_v[i+1]+s_v[i+1]*q_v[i+1]),
         n_v[i  ]*(p_v[i]  +s_v[i  ]*q_v[i  ]),
         dx
      );
      S_free_k= compute_S_free_k(
         r,
         n, s,
         p,  q,
         V,  Vp, 
         Al,  Alp,
         Bep, Bepp,
         r_Der_n, r_Der_s,
         r_Der_p,  r_Der_q
      );
   }
/*---------------------------------------------------------------------------*/
   i= nx-2;

   r= r_v[i];
   cf= 1/(1+r/cl);

   n= n_v[i];
   s= s_v[i];
   p= p_v[i];
   q= q_v[i];

   r_Der_n= pow(cf,2)*Dx_ptm1_4th(n_v[i+1],n_v[i],n_v[i-1],n_v[i-2],n_v[i-3],dx);
   r_Der_s= pow(cf,2)*Dx_ptm1_4th(s_v[i+1],s_v[i],s_v[i-1],s_v[i-2],s_v[i-3],dx);
   r_Der_p= pow(cf,2)*Dx_ptm1_4th(p_v[i+1],p_v[i],p_v[i-1],p_v[i-2],p_v[i-3],dx);
   r_Der_q= pow(cf,2)*Dx_ptm1_4th(q_v[i+1],q_v[i],q_v[i-1],q_v[i-2],q_v[i-3],dx);

   f_k[i]=
      n*(p+s*q)
   ;

   V=    V_v[i];
   Vp=  Vp_v[i];
   Al=   Al_v[i];
   Alp= Alp_v[i];
   Bep=   Bep_v[i];
   Bepp= Bepp_v[i];

   p_k[i]= compute_p_k(
      r,
      n, s,
      p,  q,
      V,   Vp, 
      Al,  Alp,
      Bep, Bepp,
      r_Der_n,  r_Der_s,
      r_Der_p,  r_Der_q
   );
   q_k[i]= pow(cf,2)*Dx_ptm1_4th(
      n_v[i+1]*(p_v[i+1]+s_v[i+1]*q_v[i+1]),
      n_v[i  ]*(p_v[i  ]+s_v[i  ]*q_v[i  ]),
      n_v[i-1]*(p_v[i-1]+s_v[i-1]*q_v[i-1]),
      n_v[i-2]*(p_v[i-2]+s_v[i-2]*q_v[i-2]),
      n_v[i-3]*(p_v[i-3]+s_v[i-3]*q_v[i-3]),
      dx
   );
/*---------------------------------------------------------------------------*/
/* regularity condition at spatin infinity: no evolution */
/*---------------------------------------------------------------------------*/
   f_k[nx-1]= 0;
   p_k[nx-1]= 0;
   q_k[nx-1]= 0; 
/*---------------------------------------------------------------------------*/
/* set the k vectors for RK4 method */
/*---------------------------------------------------------------------------*/
   switch(vn)
   {
   case 1: 
      for (int i=0; i<nx; ++i) {
         f_k1[i]= dt*f_k[i];
         p_k1[i]= dt*p_k[i];
         q_k1[i]= dt*q_k[i];
      }
      S_free_k1= dt*S_free_k;
   break;
   case 2:
      for (int i=0; i<nx; ++i) {
         f_k2[i]= dt*f_k[i];
         p_k2[i]= dt*p_k[i];
         q_k2[i]= dt*q_k[i];
      }
      S_free_k2= dt*S_free_k;
   break;
   case 3:
      for (int i=0; i<nx; ++i) {
         f_k3[i]= dt*f_k[i];
         p_k3[i]= dt*p_k[i];
         q_k3[i]= dt*q_k[i];
      }
      S_free_k3= dt*S_free_k;
   break;
   case 4:
      for (int i=0; i<nx; ++i) {
         f_k4[i]= dt*f_k[i];
         p_k4[i]= dt*p_k[i];
         q_k4[i]= dt*q_k[i];
      }
      S_free_k4= dt*S_free_k;
   break;
   default:
      cout << "ERROR: vn = " << vn << endl;
      std::quick_exit(0);
   break;
   }
/*---------------------------------------------------------------------------*/
   return;
}
/*===========================================================================*/
/* evolve with RK4 method */
/*===========================================================================*/
void EdGB::time_step(const int exc_i, 
   Field &n, Field &s, Field &f, Field &p, Field &q)
{
   KO_filter(exc_i,nx,"even",f.n);
   KO_filter(exc_i,nx,"even",p.n);
   KO_filter(exc_i,nx,"odd" ,q.n);
/*---------------------------------------------------------------------------*/
   int start_i= (exc_i>0) ? exc_i : 1;
/*---------------------------------------------------------------------------*/
   compute_fpq_ki(1, exc_i,
      n.n, s.n, f.n, p.n, q.n
   );
   for (int i=start_i; i<nx; ++i) {
      f.inter_2[i]= f.n[i]+0.5*f_k1[i];
      p.inter_2[i]= p.n[i]+0.5*p_k1[i];
      q.inter_2[i]= q.n[i]+0.5*q_k1[i];
   }
   if (exc_i==0) {
      f.inter_2[0]= make_Dx_zero(f.inter_2[4], f.inter_2[3], f.inter_2[2], f.inter_2[1]);
      p.inter_2[0]= make_Dx_zero(p.inter_2[4], p.inter_2[3], p.inter_2[2], p.inter_2[1]);
      q.inter_2[0]= 0; 
   } else {
      s.inter_2[exc_i]= s.n[exc_i]+0.5*S_free_k1;
   }
   solve_for_metric_relaxation(exc_i,
      f.inter_2,  p.inter_2, q.inter_2,
      n.inter_2,s.inter_2
   );
/*---------------------------------------------------------------------------*/
   compute_fpq_ki(2, exc_i,
      n.inter_2, s.inter_2, f.inter_2, p.inter_2, q.inter_2
   );
   for (int i=start_i; i<nx; ++i) {
      f.inter_3[i]= f.n[i]+0.5*f_k2[i];
      p.inter_3[i]= p.n[i]+0.5*p_k2[i];
      q.inter_3[i]= q.n[i]+0.5*q_k2[i];
   }
   if (exc_i==0) {
      f.inter_3[0]= make_Dx_zero(f.inter_3[4], f.inter_3[3], f.inter_3[2], f.inter_3[1]);
      p.inter_3[0]= make_Dx_zero(p.inter_3[4], p.inter_3[3], p.inter_3[2], p.inter_3[1]);
      q.inter_3[0]= 0; 
   } else {
      s.inter_3[exc_i]= s.n[exc_i]+0.5*S_free_k2;
   }
   solve_for_metric_relaxation(exc_i,
      f.inter_3,  p.inter_3, q.inter_3,
      n.inter_3,s.inter_3
   );
/*---------------------------------------------------------------------------*/
   compute_fpq_ki(3, exc_i,
      n.inter_3, s.inter_3, f.inter_3, p.inter_3, q.inter_3
   );
   for (int i=start_i; i<nx; ++i) {
      f.inter_4[i]= f.n[i]+0.5*f_k3[i];
      p.inter_4[i]= p.n[i]+0.5*p_k3[i];
      q.inter_4[i]= q.n[i]+0.5*q_k3[i];
   }
   if (exc_i==0) {
      f.inter_4[0]= make_Dx_zero(f.inter_4[4], f.inter_4[3], f.inter_4[2], f.inter_4[1]);
      p.inter_4[0]= make_Dx_zero(p.inter_4[4], p.inter_4[3], p.inter_4[2], p.inter_4[1]);
      q.inter_4[0]= 0; 
   } else {
      s.inter_4[exc_i]= s.n[exc_i]+0.5*S_free_k3;
   }
   solve_for_metric_relaxation(exc_i,
      f.inter_4,  p.inter_4, q.inter_4,
      n.inter_4,s.inter_4
   );
/*---------------------------------------------------------------------------*/
   compute_fpq_ki(4, exc_i,
      n.inter_4, s.inter_4, f.inter_4, p.inter_4, q.inter_4
   );
   for (int i=start_i; i<nx; ++i) {
      f.np1[i]= f.n[i] + (f_k1[i]+2.*f_k2[i]+2.*f_k3[i]+f_k4[i])/6.;
      p.np1[i]= p.n[i] + (p_k1[i]+2.*p_k2[i]+2.*p_k3[i]+p_k4[i])/6.;
      q.np1[i]= q.n[i] + (q_k1[i]+2.*q_k2[i]+2.*q_k3[i]+q_k4[i])/6.;
   }
   if (exc_i==0) {
      f.np1[0]= make_Dx_zero(f.np1[4], f.np1[3], f.np1[2], f.np1[1]);
      p.np1[0]= make_Dx_zero(p.np1[4], p.np1[3], p.np1[2], p.np1[1]);
      q.np1[0]= 0; 
   } else {
      s.np1[exc_i]= s.n[exc_i]+(S_free_k1+2.*S_free_k2+2.*S_free_k3+S_free_k4)/6.;
   }
   solve_for_metric_relaxation(exc_i,
      f.np1,  p.np1, q.np1,
      n.np1,s.np1
   );
/*---------------------------------------------------------------------------*/
   return;
}
/*===========================================================================*/
int EdGB::compute_radial_characteristic(
   const double r,
   const double N, const double S,
   const double P, const double Q,
   const double V, 
   const double Al, const double Alp,
   const double Vp, const double Bep, const double Bepp,
   const double r_Der_N,
   const double r_Der_S,
   const double r_Der_P,
   const double r_Der_Q)
{
/*---------------------------------------------------------------------------*/
   double Qr= Q/r;
   double Sr= S/r;

   double t_Der_P= compute_p_k(
      r,
      N, S,
      P,  Q,
      V,  Vp, 
      Al, Alp,
      Bep, Bepp,
      r_Der_N,
      r_Der_S,
      r_Der_P,
      r_Der_Q
   );
/*---------------------------------------------------------------------------*/
   double ep_dtp= 
      (768*pow(r,4)*r_Der_Q*(-(pow(Bep,3)*pow(Sr,4)) + 8*pow(Bep,4)*Qr*pow(Sr,4) + 8*pow(Bep,4)*pow(Sr,5)*P))/(-1 + 8*Bep*Qr + 12*Bep*Sr*P) - (pow(r,4)*(1 - 24*Bep*Qr + 192*pow(Bep,2)*pow(Qr,2) - 512*pow(Bep,3)*pow(Qr,3) - Al*pow(Qr,2)*pow(r,2) + 24*Al*Bep*pow(Qr,3)*pow(r,2) - 192*Al*pow(Bep,2)*pow(Qr,4)*pow(r,2) + 512*Al*pow(Bep,3)*pow(Qr,5)*pow(r,2) + 32*pow(Bep,2)*pow(Qr,2)*pow(r,2)*pow(Sr,2) - 256*pow(Bep,3)*pow(Qr,3)*pow(r,2)*pow(Sr,2) - 16*Al*pow(Bep,2)*pow(Qr,4)*pow(r,4)*pow(Sr,2) + 128*Al*pow(Bep,3)*pow(Qr,5)*pow(r,4)*pow(Sr,2) - 96*pow(Bep,2)*pow(Sr,4) + 768*pow(Bep,3)*Qr*pow(Sr,4) + 768*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,2)*pow(Sr,4) - 6144*pow(Bep,3)*Bepp*pow(Qr,3)*pow(r,2)*pow(Sr,4) - 28*Bep*Sr*P + 448*pow(Bep,2)*Qr*Sr*P - 1792*pow(Bep,3)*pow(Qr,2)*Sr*P + 28*Al*Bep*pow(Qr,2)*pow(r,2)*Sr*P - 448*Al*pow(Bep,2)*pow(Qr,3)*pow(r,2)*Sr*P + 1792*Al*pow(Bep,3)*pow(Qr,4)*pow(r,2)*Sr*P - 192*pow(Bep,3)*pow(Qr,2)*pow(r,2)*pow(Sr,3)*P + 96*Al*pow(Bep,3)*pow(Qr,4)*pow(r,4)*pow(Sr,3)*P + 768*pow(Bep,3)*pow(Sr,5)*P - 6144*pow(Bep,3)*Bepp*pow(Qr,2)*pow(r,2)*pow(Sr,5)*P + 3*Al*pow(P,2) - 72*Al*Bep*Qr*pow(P,2) + 576*Al*pow(Bep,2)*pow(Qr,2)*pow(P,2) - 1536*Al*pow(Bep,3)*pow(Qr,3)*pow(P,2) + 288*pow(Bep,2)*pow(Sr,2)*pow(P,2) - 2304*pow(Bep,3)*Qr*pow(Sr,2)*pow(P,2) - 288*Al*pow(Bep,2)*pow(Qr,2)*pow(r,2)*pow(Sr,2)*pow(P,2) + 2304*Al*pow(Bep,3)*pow(Qr,3)*pow(r,2)*pow(Sr,2)*pow(P,2) - 84*Al*Bep*Sr*pow(P,3) + 1344*Al*pow(Bep,2)*Qr*Sr*pow(P,3) - 5376*Al*pow(Bep,3)*pow(Qr,2)*Sr*pow(P,3) - 960*pow(Bep,3)*pow(Sr,3)*pow(P,3) + 960*Al*pow(Bep,3)*pow(Qr,2)*pow(r,2)*pow(Sr,3)*pow(P,3) + 816*Al*pow(Bep,2)*pow(Sr,2)*pow(P,4) - 6528*Al*pow(Bep,3)*Qr*pow(Sr,2)*pow(P,4) - 2592*Al*pow(Bep,3)*pow(Sr,3)*pow(P,5) + 64*pow(Bep,2)*pow(Sr,2)*V - 512*pow(Bep,3)*Qr*pow(Sr,2)*V - 384*pow(Bep,3)*pow(Sr,3)*P*V))/(-1 + 8*Bep*Qr + 12*Bep*Sr*P)
   ;
   double ep_dtq= 
      0
   ;
   double ep_drp=
      -1536*pow(Bep,3)*pow(r,4)*r_Der_P*pow(Sr,4)*N - (768*pow(Bep,3)*pow(r,5)*r_Der_Q*pow(Sr,5)*N*(-1 + 4*Bep*Qr + 8*Bep*Sr*P))/(-1 + 8*Bep*Qr + 12*Bep*Sr*P) + (pow(r,5)*N*(-Sr + 36*Bep*Qr*Sr - 480*pow(Bep,2)*pow(Qr,2)*Sr + 2816*pow(Bep,3)*pow(Qr,3)*Sr - 6144*pow(Bep,4)*pow(Qr,4)*Sr + Al*pow(Qr,2)*pow(r,2)*Sr - 36*Al*Bep*pow(Qr,3)*pow(r,2)*Sr + 480*Al*pow(Bep,2)*pow(Qr,4)*pow(r,2)*Sr - 2816*Al*pow(Bep,3)*pow(Qr,5)*pow(r,2)*Sr + 6144*Al*pow(Bep,4)*pow(Qr,6)*pow(r,2)*Sr - 32*pow(Bep,2)*pow(Qr,2)*pow(r,2)*pow(Sr,3) + 256*pow(Bep,3)*pow(Qr,3)*pow(r,2)*pow(Sr,3) + 16*Al*pow(Bep,2)*pow(Qr,4)*pow(r,4)*pow(Sr,3) - 128*Al*pow(Bep,3)*pow(Qr,5)*pow(r,4)*pow(Sr,3) + 96*pow(Bep,2)*pow(Sr,5) - 1152*pow(Bep,3)*Qr*pow(Sr,5) + 3072*pow(Bep,4)*pow(Qr,2)*pow(Sr,5) - 768*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,2)*pow(Sr,5) + 9216*pow(Bep,3)*Bepp*pow(Qr,3)*pow(r,2)*pow(Sr,5) - 24576*pow(Bep,4)*Bepp*pow(Qr,4)*pow(r,2)*pow(Sr,5) - 4*Al*Qr*P + 128*Al*Bep*pow(Qr,2)*P - 1536*Al*pow(Bep,2)*pow(Qr,3)*P + 8192*Al*pow(Bep,3)*pow(Qr,4)*P - 16384*Al*pow(Bep,4)*pow(Qr,5)*P + 36*Bep*pow(Sr,2)*P - 1104*pow(Bep,2)*Qr*pow(Sr,2)*P + 10752*pow(Bep,3)*pow(Qr,2)*pow(Sr,2)*P - 33792*pow(Bep,4)*pow(Qr,3)*pow(Sr,2)*P - 36*Al*Bep*pow(Qr,2)*pow(r,2)*pow(Sr,2)*P + 1104*Al*pow(Bep,2)*pow(Qr,3)*pow(r,2)*pow(Sr,2)*P - 10752*Al*pow(Bep,3)*pow(Qr,4)*pow(r,2)*pow(Sr,2)*P + 33792*Al*pow(Bep,4)*pow(Qr,5)*pow(r,2)*pow(Sr,2)*P - 1536*pow(Bep,2)*Bepp*Qr*pow(Sr,4)*P + 24576*pow(Bep,3)*Bepp*pow(Qr,2)*pow(Sr,4)*P - 98304*pow(Bep,4)*Bepp*pow(Qr,3)*pow(Sr,4)*P + 448*pow(Bep,3)*pow(Qr,2)*pow(r,2)*pow(Sr,4)*P - 1280*pow(Bep,4)*pow(Qr,3)*pow(r,2)*pow(Sr,4)*P - 224*Al*pow(Bep,3)*pow(Qr,4)*pow(r,4)*pow(Sr,4)*P + 640*Al*pow(Bep,4)*pow(Qr,5)*pow(r,4)*pow(Sr,4)*P - 1536*pow(Bep,3)*pow(Sr,6)*P + 9216*pow(Bep,4)*Qr*pow(Sr,6)*P + 12288*pow(Bep,3)*Bepp*pow(Qr,2)*pow(r,2)*pow(Sr,6)*P - 73728*pow(Bep,4)*Bepp*pow(Qr,3)*pow(r,2)*pow(Sr,6)*P - 3*Al*Sr*pow(P,2) + 252*Al*Bep*Qr*Sr*pow(P,2) - 4896*Al*pow(Bep,2)*pow(Qr,2)*Sr*pow(P,2) + 36096*Al*pow(Bep,3)*pow(Qr,3)*Sr*pow(P,2) - 92160*Al*pow(Bep,4)*pow(Qr,4)*Sr*pow(P,2) - 512*pow(Bep,2)*pow(Sr,3)*pow(P,2) + 11520*pow(Bep,3)*Qr*pow(Sr,3)*pow(P,2) - 59392*pow(Bep,4)*pow(Qr,2)*pow(Sr,3)*pow(P,2) + 512*Al*pow(Bep,2)*pow(Qr,2)*pow(r,2)*pow(Sr,3)*pow(P,2) - 11520*Al*pow(Bep,3)*pow(Qr,3)*pow(r,2)*pow(Sr,3)*pow(P,2) + 59392*Al*pow(Bep,4)*pow(Qr,4)*pow(r,2)*pow(Sr,3)*pow(P,2) + 30720*pow(Bep,3)*Bepp*Qr*pow(Sr,5)*pow(P,2) - 245760*pow(Bep,4)*Bepp*pow(Qr,2)*pow(Sr,5)*pow(P,2) - 1536*pow(Bep,4)*pow(Qr,2)*pow(r,2)*pow(Sr,5)*pow(P,2) + 768*Al*pow(Bep,4)*pow(Qr,4)*pow(r,4)*pow(Sr,5)*pow(P,2) + 6144*pow(Bep,4)*pow(Sr,7)*pow(P,2) - 49152*pow(Bep,4)*Bepp*pow(Qr,2)*pow(r,2)*pow(Sr,7)*pow(P,2) + 108*Al*Bep*pow(Sr,2)*pow(P,3) - 4976*Al*pow(Bep,2)*Qr*pow(Sr,2)*pow(P,3) + 58880*Al*pow(Bep,3)*pow(Qr,2)*pow(Sr,2)*pow(P,3) - 207872*Al*pow(Bep,4)*pow(Qr,3)*pow(Sr,2)*pow(P,3) + 3264*pow(Bep,3)*pow(Sr,4)*pow(P,3) - 39168*pow(Bep,4)*Qr*pow(Sr,4)*pow(P,3) - 3264*Al*pow(Bep,3)*pow(Qr,2)*pow(r,2)*pow(Sr,4)*pow(P,3) + 39168*Al*pow(Bep,4)*pow(Qr,3)*pow(r,2)*pow(Sr,4)*pow(P,3) - 147456*pow(Bep,4)*Bepp*Qr*pow(Sr,6)*pow(P,3) - 1488*Al*pow(Bep,2)*pow(Sr,3)*pow(P,4) + 40320*Al*pow(Bep,3)*Qr*pow(Sr,3)*pow(P,4) - 227328*Al*pow(Bep,4)*pow(Qr,2)*pow(Sr,3)*pow(P,4) - 7680*pow(Bep,4)*pow(Sr,5)*pow(P,4) + 7680*Al*pow(Bep,4)*pow(Qr,2)*pow(r,2)*pow(Sr,5)*pow(P,4) + 9120*Al*pow(Bep,3)*pow(Sr,4)*pow(P,5) - 115584*Al*pow(Bep,4)*Qr*pow(Sr,4)*pow(P,5) - 20736*Al*pow(Bep,4)*pow(Sr,5)*pow(P,6) - 64*pow(Bep,2)*pow(Sr,3)*V + 512*pow(Bep,3)*Qr*pow(Sr,3)*V + 896*pow(Bep,3)*pow(Sr,4)*P*V - 2560*pow(Bep,4)*Qr*pow(Sr,4)*P*V - 3072*pow(Bep,4)*pow(Sr,5)*pow(P,2)*V))/((-1 + 8*Bep*Qr + 8*Bep*Sr*P)*(-1 + 8*Bep*Qr + 12*Bep*Sr*P))
   ;
   double ep_drq=
      (-768*pow(Bep,3)*pow(r,5)*r_Der_P*pow(Sr,5)*N*(-1 + 4*Bep*Qr + 8*Bep*Sr*P))/(-1 + 8*Bep*Qr + 12*Bep*Sr*P) + (768*pow(Bep,3)*pow(r,4)*pow(Sr,4)*t_Der_P*(-1 + 8*Bep*Qr + 8*Bep*Sr*P))/(-1 + 8*Bep*Qr + 12*Bep*Sr*P) + (pow(r,4)*N*(-1 + 32*Bep*Qr - 384*pow(Bep,2)*pow(Qr,2) + 2048*pow(Bep,3)*pow(Qr,3) - 4096*pow(Bep,4)*pow(Qr,4) + 3*Al*pow(Qr,2)*pow(r,2) - 96*Al*Bep*pow(Qr,3)*pow(r,2) + 1152*Al*pow(Bep,2)*pow(Qr,4)*pow(r,2) - 6144*Al*pow(Bep,3)*pow(Qr,5)*pow(r,2) + 12288*Al*pow(Bep,4)*pow(Qr,6)*pow(r,2) + 48*pow(Bep,2)*pow(Qr,2)*pow(r,2)*pow(Sr,2) - 768*pow(Bep,3)*pow(Qr,3)*pow(r,2)*pow(Sr,2) + 3072*pow(Bep,4)*pow(Qr,4)*pow(r,2)*pow(Sr,2) - 64*Al*pow(Bep,2)*pow(Qr,4)*pow(r,4)*pow(Sr,2) + 1024*Al*pow(Bep,3)*pow(Qr,5)*pow(r,4)*pow(Sr,2) - 4096*Al*pow(Bep,4)*pow(Qr,6)*pow(r,4)*pow(Sr,2) + 96*pow(Bep,2)*pow(Sr,4) - 1536*pow(Bep,3)*Qr*pow(Sr,4) + 6144*pow(Bep,4)*pow(Qr,2)*pow(Sr,4) - 256*pow(Bep,4)*pow(Qr,4)*pow(r,4)*pow(Sr,4) + 128*Al*pow(Bep,4)*pow(Qr,6)*pow(r,6)*pow(Sr,4) + 32*Bep*Sr*P - 768*pow(Bep,2)*Qr*Sr*P + 6144*pow(Bep,3)*pow(Qr,2)*Sr*P - 16384*pow(Bep,4)*pow(Qr,3)*Sr*P - 112*Al*Bep*pow(Qr,2)*pow(r,2)*Sr*P + 2688*Al*pow(Bep,2)*pow(Qr,3)*pow(r,2)*Sr*P - 21504*Al*pow(Bep,3)*pow(Qr,4)*pow(r,2)*Sr*P + 57344*Al*pow(Bep,4)*pow(Qr,5)*pow(r,2)*Sr*P - 896*pow(Bep,3)*pow(Qr,2)*pow(r,2)*pow(Sr,3)*P + 7168*pow(Bep,4)*pow(Qr,3)*pow(r,2)*pow(Sr,3)*P + 1152*Al*pow(Bep,3)*pow(Qr,4)*pow(r,4)*pow(Sr,3)*P - 9216*Al*pow(Bep,4)*pow(Qr,5)*pow(r,4)*pow(Sr,3)*P - 1536*pow(Bep,3)*pow(Sr,5)*P + 12288*pow(Bep,4)*Qr*pow(Sr,5)*P - 3072*pow(Bep,3)*Bepp*pow(Qr,2)*pow(r,2)*pow(Sr,5)*P + 24576*pow(Bep,4)*Bepp*pow(Qr,3)*pow(r,2)*pow(Sr,5)*P - Al*pow(P,2) + 32*Al*Bep*Qr*pow(P,2) - 384*Al*pow(Bep,2)*pow(Qr,2)*pow(P,2) + 2048*Al*pow(Bep,3)*pow(Qr,3)*pow(P,2) - 4096*Al*pow(Bep,4)*pow(Qr,4)*pow(P,2) - 352*pow(Bep,2)*pow(Sr,2)*pow(P,2) + 5632*pow(Bep,3)*Qr*pow(Sr,2)*pow(P,2) - 22528*pow(Bep,4)*pow(Qr,2)*pow(Sr,2)*pow(P,2) + 1616*Al*pow(Bep,2)*pow(Qr,2)*pow(r,2)*pow(Sr,2)*pow(P,2) - 25856*Al*pow(Bep,3)*pow(Qr,3)*pow(r,2)*pow(Sr,2)*pow(P,2) + 103424*Al*pow(Bep,4)*pow(Qr,4)*pow(r,2)*pow(Sr,2)*pow(P,2) + 768*pow(Bep,2)*Bepp*pow(Sr,4)*pow(P,2) - 12288*pow(Bep,3)*Bepp*Qr*pow(Sr,4)*pow(P,2) + 49152*pow(Bep,4)*Bepp*pow(Qr,2)*pow(Sr,4)*pow(P,2) + 3840*pow(Bep,4)*pow(Qr,2)*pow(r,2)*pow(Sr,4)*pow(P,2) - 4864*Al*pow(Bep,4)*pow(Qr,4)*pow(r,4)*pow(Sr,4)*pow(P,2) + 6144*pow(Bep,4)*pow(Sr,6)*pow(P,2) + 24576*pow(Bep,4)*Bepp*pow(Qr,2)*pow(r,2)*pow(Sr,6)*pow(P,2) + 32*Al*Bep*Sr*pow(P,3) - 768*Al*pow(Bep,2)*Qr*Sr*pow(P,3) + 6144*Al*pow(Bep,3)*pow(Qr,2)*Sr*pow(P,3) - 16384*Al*pow(Bep,4)*pow(Qr,3)*Sr*pow(P,3) + 1536*pow(Bep,3)*pow(Sr,3)*pow(P,3) - 12288*pow(Bep,4)*Qr*pow(Sr,3)*pow(P,3) - 10624*Al*pow(Bep,3)*pow(Qr,2)*pow(r,2)*pow(Sr,3)*pow(P,3) + 84992*Al*pow(Bep,4)*pow(Qr,3)*pow(r,2)*pow(Sr,3)*pow(P,3) - 12288*pow(Bep,3)*Bepp*pow(Sr,5)*pow(P,3) + 98304*pow(Bep,4)*Bepp*Qr*pow(Sr,5)*pow(P,3) - 368*Al*pow(Bep,2)*pow(Sr,2)*pow(P,4) + 5888*Al*pow(Bep,3)*Qr*pow(Sr,2)*pow(P,4) - 23552*Al*pow(Bep,4)*pow(Qr,2)*pow(Sr,2)*pow(P,4) - 2048*pow(Bep,4)*pow(Sr,4)*pow(P,4) + 26240*Al*pow(Bep,4)*pow(Qr,2)*pow(r,2)*pow(Sr,4)*pow(P,4) + 49152*pow(Bep,4)*Bepp*pow(Sr,6)*pow(P,4) + 1792*Al*pow(Bep,3)*pow(Sr,3)*pow(P,5) - 14336*Al*pow(Bep,4)*Qr*pow(Sr,3)*pow(P,5) - 3072*Al*pow(Bep,4)*pow(Sr,4)*pow(P,6) - 64*pow(Bep,2)*pow(Sr,2)*V + 1024*pow(Bep,3)*Qr*pow(Sr,2)*V - 4096*pow(Bep,4)*pow(Qr,2)*pow(Sr,2)*V - 512*pow(Bep,4)*pow(Qr,2)*pow(r,2)*pow(Sr,4)*V + 1024*pow(Bep,3)*pow(Sr,3)*P*V - 8192*pow(Bep,4)*Qr*pow(Sr,3)*P*V - 4096*pow(Bep,4)*pow(Sr,4)*pow(P,2)*V))/((-1 + 8*Bep*Qr + 8*Bep*Sr*P)*(-1 + 8*Bep*Qr + 12*Bep*Sr*P))
   ;
/*---------------------------------------------------------------------------*/
   double eq_dtp=
      0
   ;
   double eq_dtq=
      1
   ;
   double eq_drp=
   -   N
   ;
   double eq_drq=
   -   N*S
   ;
/*---------------------------------------------------------------------------*/
   double A= ep_dtp*eq_dtq- ep_dtq*eq_dtp;
   double B= 
   -  (ep_dtp*eq_drq-ep_dtq*eq_drp)
   -  (ep_drp*eq_dtq-ep_drq*eq_dtp)
   ;
   double C= ep_drp*eq_drq-ep_drq*eq_drp;

   double D=  pow(B,2)-4*A*C;

   if (D<0) {
      ingoing_c=  0;
      outgoing_c= 0;
      return -1;
   }
/*---------------------------------------------------------------------------*/
   ingoing_c= (-B-sqrt(D))/(2*A);
   outgoing_c=(-B+sqrt(D))/(2*A);
/*---------------------------------------------------------------------------*/
   return 0;
}
/*===========================================================================*/
void EdGB::compute_radial_characteristics(
   int &exc_i,
   const vector<double> &n_v, const vector<double> &s_v,
   const vector<double> &f_v, 
   const vector<double> &p_v, const vector<double> &q_v,
   vector<double> &ingoing,
   vector<double> &outgoing)
{
   compute_potentials(f_v);
   int new_exc_i= exc_i;
/*---------------------------------------------------------------------------*/
/* interior two grid points */
   int i=exc_i;
   double x= dx*i;
   double cf= (1-x/cl);
   double r= x/cf;

   double N= n_v[i];
   double S= s_v[i];
   double P= p_v[i];
   double Q= q_v[i];

   double V=    V_v[i];
   double Vp=   Vp_v[i];
   double Al=   Al_v[i];
   double Alp=  Alp_v[i];
   double Bep=  Bep_v[i];
   double Bepp= Bepp_v[i];

   double r_Der_N= pow(cf,2)*Dx_ptp0_4th(n_v[i+4],n_v[i+3],n_v[i+2],n_v[i+1],n_v[i],dx);
   double r_Der_S= pow(cf,2)*Dx_ptp0_4th(s_v[i+4],s_v[i+3],s_v[i+2],s_v[i+1],s_v[i],dx);
   double r_Der_P= pow(cf,2)*Dx_ptp0_4th(p_v[i+4],p_v[i+3],p_v[i+2],p_v[i+1],p_v[i],dx);
   double r_Der_Q= pow(cf,2)*Dx_ptp0_4th(q_v[i+4],q_v[i+3],q_v[i+2],q_v[i+1],q_v[i],dx);

   int status= compute_radial_characteristic(
      r,
      N,  S,
      P,   Q,
      V,   Vp,  
      Al,  Alp,
      Bep, Bepp,
      r_Der_N,
      r_Der_S,
      r_Der_P,
      r_Der_Q
   );
   if (status==-1) {
      new_exc_i= i+1;
   }
   ingoing[i]=   ingoing_c;
   outgoing[i]= outgoing_c;
/*--------------------------------------*/
   i=exc_i+1;
   x= dx*i;
   cf= (1-x/cl);
   r= x/cf;

   N= n_v[i];
   S= s_v[i];
   P= p_v[i];
   Q= q_v[i];

   V=   V_v[i];
   Vp=  Vp_v[i];
   Bep=  Bep_v[i];
   Bepp= Bepp_v[i];

   r_Der_N= pow(cf,2)*Dx_ptp1_4th(n_v[i+4],n_v[i+3],n_v[i+2],n_v[i+1],n_v[i],dx);
   r_Der_S= pow(cf,2)*Dx_ptp1_4th(s_v[i+4],s_v[i+3],s_v[i+2],s_v[i+1],s_v[i],dx);
   r_Der_P= pow(cf,2)*Dx_ptp1_4th(p_v[i+4],p_v[i+3],p_v[i+2],p_v[i+1],p_v[i],dx);
   r_Der_Q= pow(cf,2)*Dx_ptp1_4th(q_v[i+4],q_v[i+3],q_v[i+2],q_v[i+1],q_v[i],dx);

   status= compute_radial_characteristic(
      r,
      N,  S,
      P,  Q,
      V,  Vp,  
      Al, Alp,
      Bep, Bepp,
      r_Der_N,
      r_Der_S,
      r_Der_P,
      r_Der_Q
   );
   if (status==-1) {
      new_exc_i= i+1;
   }
   ingoing[i]=   ingoing_c;
   outgoing[i]= outgoing_c;
/*---------------------------------------------------------------------------*/
/* interior */
   for (i=exc_i+2; i<nx-2; ++i) {
      x= dx*i;
      cf= (1-x/cl);
      r= x/cf;

      N= n_v[i];
      S= s_v[i];
      P= p_v[i];
      Q= q_v[i];

      V=    V_v[i];
      Vp=   Vp_v[i];
      Al=   Al_v[i];
      Bep=  Bep_v[i];
      Bepp= Bepp_v[i];

      r_Der_N= pow(cf,2)*Dx_ptc_4th(n_v[i+2],n_v[i+1],n_v[i-1],n_v[i-2],dx);
      r_Der_S= pow(cf,2)*Dx_ptc_4th(s_v[i+2],s_v[i+1],s_v[i-1],s_v[i-2],dx);
      r_Der_P= pow(cf,2)*Dx_ptc_4th(p_v[i+2],p_v[i+1],p_v[i-1],p_v[i-2],dx);
      r_Der_Q= pow(cf,2)*Dx_ptc_4th(q_v[i+2],q_v[i+1],q_v[i-1],q_v[i-2],dx);

      status= compute_radial_characteristic(
         r,
         N, S,
         P, Q,
         V,  Vp,  
         Al, Alp,
         Bep,  Bepp,
         r_Der_N,
         r_Der_S,
         r_Der_P,
         r_Der_Q
      );
      if (status==-1) {
         new_exc_i= i+1;
      }
      ingoing[i]=   ingoing_c;
      outgoing[i]= outgoing_c;
   }
   exc_i= new_exc_i;

   if (outgoing[exc_i+1]>0) {
      cout<<"naked_elliptic_region"<<endl;
      std::quick_exit(0);
   }

   return;
}
/*===========================================================================*/
void EdGB::compute_eom_rr(
	const int exc_i,
	const vector<double> &n_v, const vector<double> &s_v,
	const vector<double> &f_v, 
	const vector<double> &p_v, const vector<double> &q_v,
	vector<double> &eom_rr)
{
   compute_potentials(f_v);
   for (int i=exc_i+2; i<nx-2; ++i) {
      double x= dx*i;
      double cf= (1-x/cl);
      double r= x/cf;

      double N= n_v[i];
      double S= s_v[i];
      double P= p_v[i];
      double Q= q_v[i];

      double V=    V_v[i];
      double Vp=   Vp_v[i];
      double Al=   Al_v[i];
      double Alp=  Alp_v[i];
      double Bep=  Bep_v[i];
      double Bepp= Bepp_v[i];

      double Sr= S/r;
      double Qr= Q/r;

      double r_Der_N= pow(cf,2)*Dx_ptc_4th(n_v[i+2],n_v[i+1],n_v[i-1],n_v[i-2],dx);
      double r_Der_S= pow(cf,2)*Dx_ptc_4th(s_v[i+2],s_v[i+1],s_v[i-1],s_v[i-2],dx);
      double r_Der_P= pow(cf,2)*Dx_ptc_4th(p_v[i+2],p_v[i+1],p_v[i-1],p_v[i-2],dx);
      double r_Der_Q= pow(cf,2)*Dx_ptc_4th(q_v[i+2],q_v[i+1],q_v[i-1],q_v[i-2],dx);

      double t_Der_P= compute_p_k(
         r,
         N, S,
         P, Q,
         V,   Vp, 
         Al,  Alp,
         Bep, Bepp,
         r_Der_N,
         r_Der_S,
         r_Der_P,
         r_Der_Q
      );
      double t_Der_S= compute_S_free_k(
         r,
         N, S,
         P, Q,
         V, Vp, 
         Al,  Alp,
         Bep, Bepp,
         r_Der_N, r_Der_S,
         r_Der_P, r_Der_Q
      );
      eom_rr[i]=
         -pow(Sr,2) + 8*Bep*r*r_Der_P*pow(Sr,3) - (8*Bep*pow(Sr,2)*t_Der_P)/N - 8*Bepp*pow(Sr,2)*pow(P,2) + r_Der_S*(-2*Sr + 16*Bep*Qr*Sr + 16*Bep*pow(Sr,2)*P) + t_Der_S*(2/(r*N) - (16*Bep*Qr)/(r*N) - (16*Bep*Sr*P)/(r*N)) + r_Der_N*(2/(r*N) - (16*Bep*Qr)/(r*N) + (8*Bep*Qr*r*pow(Sr,2))/N - (16*Bep*Sr*P)/(r*N)) + (-2*pow(Qr,2)*pow(r,2) + 3*Al*pow(Qr,4)*pow(r,4) - 2*pow(P,2) - 2*Al*pow(Qr,2)*pow(r,2)*pow(P,2) - Al*pow(P,4) + 4*V)/4.
      ;
   }
   return;
}
/*===========================================================================*/
/* Ricci tensor contracted on two outgoing null vectors */
/*===========================================================================*/
void EdGB::compute_ncc(
   const int exc_i,
   const int nx, const double dx, const double cl, 
   const vector<double> &r_v, 
   const vector<double> &n_v, const vector<double> &s_v,
   const vector<double> &p_v, const vector<double> &q_v,
   vector<double> &ncc)
{
   for (int i=exc_i+2; i<nx-3; ++i) {
      double r= r_v[i];
      double cf= 1/(1+r/cl);

      double N= n_v[i];
      double S= s_v[i];
      double P= p_v[i];
      double Q= q_v[i];

      double V=   V_v[i];
      double Vp=  Vp_v[i];
      double Al=  Al_v[i];
      double Alp= Alp_v[i];
      double Bep=  Bep_v[i];
      double Bepp= Bepp_v[i];

      double r_Der_N= pow(cf,2)*Dx_ptc_4th(n_v[i+2],n_v[i+1],n_v[i-1],n_v[i-2],dx);
      double r_Der_S= pow(cf,2)*Dx_ptc_4th(s_v[i+2],s_v[i+1],s_v[i-1],s_v[i-2],dx);
      double r_Der_P= pow(cf,2)*Dx_ptc_4th(p_v[i+2],p_v[i+1],p_v[i-1],p_v[i-2],dx);
      double r_Der_Q= pow(cf,2)*Dx_ptc_4th(q_v[i+2],q_v[i+1],q_v[i-1],q_v[i-2],dx);

      double t_Der_S= compute_S_free_k(
         r,
         N, S,
         P, Q,
         V, Vp, 
         Al, Alp,
         Bep, Bepp,
         r_Der_N, r_Der_S,
         r_Der_P, r_Der_Q
      );

      ncc[i]= (2*t_Der_S*N)/r + r_Der_N*((2*N)/r - (4*N*S)/r + (2*N*pow(S,2))/r);
   }
   return;
}
