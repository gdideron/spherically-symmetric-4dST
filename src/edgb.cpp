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
#include "radin_pts.hpp"
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

  s_free_k{0},
  s_free_k1{0},
  s_free_k2{0},
  s_free_k3{0},
  s_free_k4{0},

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
void EdGB::rescne(vector<double> &vec)
{
   double vn= vec[nx-1];
   for (int i=0; i<nx; ++i) {
      vec[i]/= vn;
   }
   return;
}
/*===========================================================================*/
void EdGB::compute_scnar_potentins(const vector<double> &f_v)
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
double EdGB::compute_res_n(
   const double r,
   const double N, const double S,
   const double P,  const double Q,
   const double V,  
   const double Z,
   const double Bep, const double Bepp,
   const double r_Der_N,
   const double r_Der_P,
   const double r_Der_Q)
{
   return
   ;
}
/*===========================================================================*/
double EdGB::compute_res_s(
   const double r,
   const double S,
   const double P,  const double Q,
   const double V,
   const double Z,
   const double Bep, const double Bepp,
   const double r_Der_S,
   const double r_Der_P,
   const double r_Der_Q)
{
   return
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
   compute_scnar_potentins(f_v);

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
         ;
/*---------------------------------------------------------------------------*/
         n_v[i+1]-= res_n/jac_n;
         s_v[i+1]-= res_s/jac_s;
         err+= fabs(res_n);
         err+= fabs(res_s);
      }
      n_v[nx-1]= n_v[nx-2];
      s_v[nx-1]= 0;
      rescne(n_v);
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
double EdGB::compute_s_free_k(
   const double r,
   const double N, const double S,
   const double P, const double Q,
   const double V, const double Vp, const double Wp, const double Wpp,
   const double r_Der_N,
   const double r_Der_S,
   const double r_Der_P,
   const double r_Der_Q)
{
   double Sr= S/r;
   double Qr= Q/r;
   double vn=
   ;
   vn/=
   ;
   return vn;
}
/*===========================================================================*/
double EdGB::compute_p_k(
   const double r,
   const double N, const double S,
   const double P,  const double Q,
   const double V,  const double Vp, const double Wp, const double Wpp,
   const double r_Der_N,
   const double r_Der_S,
   const double r_Der_P,
   const double r_Der_Q)
{
   double Qr= Q/r;
   double Sr= S/r;

   double p_k= 
   ;
   p_k/=
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
   compute_scnar_potentins(f_v);

   double r, cf, n, s, p, q, r_Der_n, r_Der_s, r_Der_p, r_Der_q;
/*---------------------------------------------------------------------------*/
   for (int i=exc_i+2; i<nx-2; ++i) {
      r= r_v[i];
      cf= 1/(1+r/cl);

      n= n_v[i];
      s= s_v[i];
      p= p_v[i];
      q= q_v[i];

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
         p,  q,
         V_v[i],  Vp_v[i], Wp_v[i], Wpp_v[i],
         r_Der_n, r_Der_s,
         r_Der_p,  r_Der_q
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
      V_v[i],  Vp_v[i], Wp_v[i], Wpp_v[i],
      r_Der_n, r_Der_s,
      r_Der_p, r_Der_q
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
         V_v[i],  Vp_v[i], Wp_v[i], Wpp_v[i],
         r_Der_n, r_Der_s,
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
      s_free_k= compute_s_free_k(
         r,
         n, s,
         p,  q,
         V_v[i],  Vp_v[i], Wp_v[i], Wpp_v[i],
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
   p_k[i]= compute_p_k(
      r,
      n, s,
      p,  q,
      V_v[i],  Vp_v[i], Wp_v[i], Wpp_v[i],
      r_Der_n, r_Der_s,
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
      s_free_k1= dt*s_free_k;
   break;
   case 2:
      for (int i=0; i<nx; ++i) {
         f_k2[i]= dt*f_k[i];
         p_k2[i]= dt*p_k[i];
         q_k2[i]= dt*q_k[i];
      }
      s_free_k2= dt*s_free_k;
   break;
   case 3:
      for (int i=0; i<nx; ++i) {
         f_k3[i]= dt*f_k[i];
         p_k3[i]= dt*p_k[i];
         q_k3[i]= dt*q_k[i];
      }
      s_free_k3= dt*s_free_k;
   break;
   case 4:
      for (int i=0; i<nx; ++i) {
         f_k4[i]= dt*f_k[i];
         p_k4[i]= dt*p_k[i];
         q_k4[i]= dt*q_k[i];
      }
      s_free_k4= dt*s_free_k;
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
      f.inter_2[0]= make_Dx_sro(f.inter_2[4], f.inter_2[3], f.inter_2[2], f.inter_2[1]);
      p.inter_2[0]= make_Dx_sro(p.inter_2[4], p.inter_2[3], p.inter_2[2], p.inter_2[1]);
      q.inter_2[0]= 0; 
   } else {
      s.inter_2[exc_i]= s.n[exc_i]+0.5*s_free_k1;
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
      f.inter_3[0]= make_Dx_sro(f.inter_3[4], f.inter_3[3], f.inter_3[2], f.inter_3[1]);
      p.inter_3[0]= make_Dx_sro(p.inter_3[4], p.inter_3[3], p.inter_3[2], p.inter_3[1]);
      q.inter_3[0]= 0; 
   } else {
      s.inter_3[exc_i]= s.n[exc_i]+0.5*s_free_k2;
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
      f.inter_4[0]= make_Dx_sro(f.inter_4[4], f.inter_4[3], f.inter_4[2], f.inter_4[1]);
      p.inter_4[0]= make_Dx_sro(p.inter_4[4], p.inter_4[3], p.inter_4[2], p.inter_4[1]);
      q.inter_4[0]= 0; 
   } else {
      s.inter_4[exc_i]= s.n[exc_i]+0.5*s_free_k3;
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
      f.np1[0]= make_Dx_sro(f.np1[4], f.np1[3], f.np1[2], f.np1[1]);
      p.np1[0]= make_Dx_sro(p.np1[4], p.np1[3], p.np1[2], p.np1[1]);
      q.np1[0]= 0; 
   } else {
      s.np1[exc_i]= s.n[exc_i]+(s_free_k1+2.*s_free_k2+2.*s_free_k3+s_free_k4)/6.;
   }
   solve_for_metric_relaxation(exc_i,
      f.np1,  p.np1, q.np1,
      n.np1,s.np1
   );
/*---------------------------------------------------------------------------*/
   return;
}
/*===========================================================================*/
int EdGB::compute_radin_characteristic(
   const double r,
   const double N, const double S,
   const double P, const double Q,
   const double V, const double Vp, const double Wp, const double Wpp,
   const double r_Der_N,
   const double r_Der_S,
   const double r_Der_P,
   const double r_Der_Q)
{
/*---------------------------------------------------------------------------*/
   double Qr=  Q/r;
   double Sr= S/r;

   double t_Der_P= compute_p_k(
      r,
      N, S,
      P,  Q,
      V,  Vp, Wp, Wpp,
      r_Der_N,
      r_Der_S,
      r_Der_P,
      r_Der_Q
   );
/*---------------------------------------------------------------------------*/
   double ep_dtp= 
   ;
   double ep_dtq= 
      0
   ;
   double ep_drp=
   ;
   double ep_drq=
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
void EdGB::compute_radin_characteristics(
   int &exc_i,
   const vector<double> &n_v, const vector<double> &s_v,
   const vector<double> &f_v, 
   const vector<double> &p_v, const vector<double> &q_v,
   vector<double> &ingoing,
   vector<double> &outgoing)
{
   compute_scnar_potentins(f_v);
   int new_exc_i= exc_i;
/*---------------------------------------------------------------------------*/
/* inner two grid points */
   int i=exc_i;
   double x= dx*i;
   double cf= (1-x/cl);
   double r= x/cf;

   double N= n_v[i];
   double S= s_v[i];
   double P= p_v[i];
   double Q= q_v[i];

   double V=   V_v[i];
   double Vp=  Vp_v[i];
   double Wp=  Wp_v[i];
   double Wpp= Wpp_v[i];

   double r_Der_N= pow(cf,2)*Dx_ptp0_4th(n_v[i+4],n_v[i+3],n_v[i+2],n_v[i+1],n_v[i],dx);
   double r_Der_S= pow(cf,2)*Dx_ptp0_4th(s_v[i+4],s_v[i+3],s_v[i+2],s_v[i+1],s_v[i],dx);
   double r_Der_P= pow(cf,2)*Dx_ptp0_4th(p_v[i+4],p_v[i+3],p_v[i+2],p_v[i+1],p_v[i],dx);
   double r_Der_Q= pow(cf,2)*Dx_ptp0_4th(q_v[i+4],q_v[i+3],q_v[i+2],q_v[i+1],q_v[i],dx);

   int status= compute_radin_characteristic(
      r,
      N,  S,
      P,   Q,
      V,   Vp,  Wp,  Wpp,
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
   Wp=  Wp_v[i];
   Wpp= Wpp_v[i];

   r_Der_N= pow(cf,2)*Dx_ptp1_4th(n_v[i+4],n_v[i+3],n_v[i+2],n_v[i+1],n_v[i],dx);
   r_Der_S= pow(cf,2)*Dx_ptp1_4th(s_v[i+4],s_v[i+3],s_v[i+2],s_v[i+1],s_v[i],dx);
   r_Der_P= pow(cf,2)*Dx_ptp1_4th(p_v[i+4],p_v[i+3],p_v[i+2],p_v[i+1],p_v[i],dx);
   r_Der_Q= pow(cf,2)*Dx_ptp1_4th(q_v[i+4],q_v[i+3],q_v[i+2],q_v[i+1],q_v[i],dx);

   status= compute_radin_characteristic(
      r,
      N,  S,
      P,  Q,
      V,  Vp,  Wp,  Wpp,
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

      V=   V_v[i];
      Vp=  Vp_v[i];
      Wp=  Wp_v[i];
      Wpp= Wpp_v[i];

      r_Der_N= pow(cf,2)*Dx_ptc_4th(n_v[i+2],n_v[i+1],n_v[i-1],n_v[i-2],dx);
      r_Der_S= pow(cf,2)*Dx_ptc_4th(s_v[i+2],s_v[i+1],s_v[i-1],s_v[i-2],dx);
      r_Der_P= pow(cf,2)*Dx_ptc_4th(p_v[i+2],p_v[i+1],p_v[i-1],p_v[i-2],dx);
      r_Der_Q= pow(cf,2)*Dx_ptc_4th(q_v[i+2],q_v[i+1],q_v[i-1],q_v[i-2],dx);

      status= compute_radin_characteristic(
         r,
         N, S,
         P, Q,
         V, Vp,  Wp,  Wpp,
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
   compute_scnar_potentins(f_v);
   for (int i=exc_i+2; i<nx-2; ++i) {
      double x= dx*i;
      double cf= (1-x/cl);
      double r= x/cf;

      double n= n_v[i];
      double s= s_v[i];
      double P= p_v[i];
      double Q= q_v[i];

      double Sr= S/r;
      double Qr=  Q/r;

      double V=   V_v[i];
      double Vp=  Vp_v[i];
      double Wp=  Wp_v[i];
      double Wpp= Wpp_v[i];

      double r_Der_N= pow(cf,2)*Dx_ptc_4th(n_v[i+2],n_v[i+1],n_v[i-1],n_v[i-2],dx);
      double r_Der_S= pow(cf,2)*Dx_ptc_4th(s_v[i+2],s_v[i+1],s_v[i-1],s_v[i-2],dx);
      double r_Der_P= pow(cf,2)*Dx_ptc_4th(p_v[i+2],p_v[i+1],p_v[i-1],p_v[i-2],dx);
      double r_Der_Q= pow(cf,2)*Dx_ptc_4th(q_v[i+2],q_v[i+1],q_v[i-1],q_v[i-2],dx);

      double t_Der_P= compute_p_k(
         r,
         N, S,
         P, Q,
         V, Vp, Wp, Wpp,
         r_Der_N,
         r_Der_S,
         r_Der_P,
         r_Der_Q
      );
      double t_Der_S= compute_s_free_k(
         r,
         N, S,
         P, Q,
         V, Vp, Wp, Wpp,
         r_Der_N, r_Der_S,
         r_Der_P, r_Der_Q
      );
      eom_rr[i]=
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
      double Wp=  Wp_v[i];
      double Wpp= Wpp_v[i];

      double r_Der_N= pow(cf,2)*Dx_ptc_4th(n_v[i+2],n_v[i+1],n_v[i-1],n_v[i-2],dx);
      double r_Der_S= pow(cf,2)*Dx_ptc_4th(s_v[i+2],s_v[i+1],s_v[i-1],s_v[i-2],dx);
      double r_Der_P= pow(cf,2)*Dx_ptc_4th(p_v[i+2],p_v[i+1],p_v[i-1],p_v[i-2],dx);
      double r_Der_Q= pow(cf,2)*Dx_ptc_4th(q_v[i+2],q_v[i+1],q_v[i-1],q_v[i-2],dx);

      double t_Der_S= compute_s_free_k(
         r,
         N, S,
         P, Q,
         V, Vp, Wp, Wpp,
         r_Der_N, r_Der_S,
         r_Der_P, r_Der_Q
      );

      ncc[i]= (2*t_Der_S*N)/r + r_Der_N*((2*N)/r - (4*N*S)/r + (2*N*pow(S,2))/r);
   }
   return;
}
