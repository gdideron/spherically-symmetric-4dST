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
   const double V_1,  const double V_2,  const double V_3,  const double V_4,
   const double Al_0,
   const double Al_1, const double Al_2, const double Al_3, const double Al_4,
   const double Be_1, const double Be_2, const double Be_3, const double Be_4,
   const vector<double> &rvec)
: dt{dt},
  dx{dx},
  nx{nx},

  V_1{V_1},  
  V_2{V_2},  
  V_3{V_3},  
  V_4{V_4},  
 
  Al_0{Al_0},  
  Al_1{Al_1},  
  Al_2{Al_2},  
  Al_3{Al_3},  
  Al_4{Al_4},  
 
  Be_1{Be_1},  
  Be_2{Be_2},  
  Be_3{Be_3},  
  Be_4{Be_4},  
 
  cl{cl},
 
  r_v(nx),

  S_free_k1{0},
  S_free_k2{0},
  S_free_k3{0},
  S_free_k4{0},

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
                  V_1*f_v[i]
      +  (1./2.) *V_2*pow(f_v[i],2)
      +  (1./6.) *V_3*pow(f_v[i],3)
      +  (1./24.)*V_4*pow(f_v[i],4)
      ;
      Vp_v[i]= 
                 V_1
      +          V_2*f_v[i]
      +  (1./2.)*V_3*pow(f_v[i],2)
      +  (1./6.)*V_4*pow(f_v[i],3)
      ;
      Al_v[i]=
                  Al_0
      +           Al_1*f_v[i]
      +  (1./2.) *Al_2*pow(f_v[i],2)
      +  (1./6.) *Al_3*pow(f_v[i],3)
      +  (1./24.)*Al_4*pow(f_v[i],4)
      ;
      Alp_v[i]=
                 Al_1
      +          Al_2*f_v[i]
      +  (1./2.)*Al_3*pow(f_v[i],2)
      +  (1./6.)*Al_4*pow(f_v[i],3)
      ;
      Be_v[i]=
                  Be_1*f_v[i]
      +  (1./2.) *Be_2*pow(f_v[i],2)
      +  (1./6.) *Be_3*pow(f_v[i],3)
      +  (1./24.)*Be_4*pow(f_v[i],4)
      ;
      Bep_v[i]=
                 Be_1
      +          Be_2*f_v[i]
      +  (1./2.)*Be_3*pow(f_v[i],2)
      +  (1./6.)*Be_4*pow(f_v[i],3)
      ;
      Bepp_v[i]=
                 Be_2
      +          Be_3*f_v[i]
      +  (1./2.)*Be_4*pow(f_v[i],2)
      ;
   }
   return;
}
/*===========================================================================*/
/* solve for metric using relaxation */
/*===========================================================================*/
void EdGB::solve_for_metric_relaxation(
   const int exc_i,
   const vector<double> &p_v,
   const vector<double> &q_v,
   vector<double> &n_v,
   vector<double> &s_v)
{
   const double err_tolerance= 1e-12;
   double err= 0;
   int iters= 0;
   do {
      iters++;
      for (int i=exc_i; i<nx-2; ++i) {
         double x=  dx*(2*i+1)/2;

         double cf= 1-x/cl;
         double r=  x/cf;
         double dr= dx*pow(cf,-2);

         double P= (p_v[i+1] + p_v[i])/2;
         double Q= (q_v[i+1] + q_v[i])/2;
         double N= (n_v[i+1] + n_v[i])/2;
         double S= (s_v[i+1] + s_v[i])/2;

         double V=    (   V_v[i+1] +    V_v[i])/2;
         double Al=   (  Al_v[i+1] +   Al_v[i])/2;
         double Bep=  ( Bep_v[i+1] +  Bep_v[i])/2;
         double Bepp= (Bepp_v[i+1] + Bepp_v[i])/2;

         double r_Der_P= (p_v[i+1] - p_v[i])/dr;
         double r_Der_Q= (q_v[i+1] - q_v[i])/dr;
         double r_Der_N= (n_v[i+1] - n_v[i])/dr;
         double r_Der_S= (s_v[i+1] - s_v[i])/dr;
/*---------------------------------------------------------------------------*/
         double res_N= compute_res_N(
            r,
            N, S, 
            P, Q,
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
            S, 
            P,  Q,
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
   cout<<"iters "<<iters<<endl;
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
   compute_potentials(phi_f.n);
   solve_for_metric_relaxation(exc_i,
      phi_p.n,phi_q.n,
      n.n,    s.n
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
void EdGB::compute_fpqS_ki(
   const int exc_i,
   const vector<double> &n_v,
   const vector<double> &s_v,
   const vector<double> &p_v, 
   const vector<double> &q_v,
   vector<double> &f_k,
   vector<double> &p_k,
   vector<double> &q_k,
   double &S_free_k)
{
/*---------------------------------------------------------------------------*/
   for (int i=exc_i+2; i<nx-2; ++i) {
      double cf= 1/(1+r_v[i]/cl);

      double r_Der_n= pow(cf,2)*Dx_ptc_4th(n_v[i+2],n_v[i+1],n_v[i-1],n_v[i-2],dx);
      double r_Der_s= pow(cf,2)*Dx_ptc_4th(s_v[i+2],s_v[i+1],s_v[i-1],s_v[i-2],dx);
      double r_Der_p= pow(cf,2)*Dx_ptc_4th(p_v[i+2],p_v[i+1],p_v[i-1],p_v[i-2],dx);
      double r_Der_q= pow(cf,2)*Dx_ptc_4th(q_v[i+2],q_v[i+1],q_v[i-1],q_v[i-2],dx);

      f_k[i]= 
         n_v[i]*(p_v[i]+s_v[i]*q_v[i])
      ;
      p_k[i]= compute_p_k(
         r_v[i],
         n_v[i], s_v[i],
         p_v[i], q_v[i],
         V_v[i], Vp_v[i], 
         Al_v[i],  Alp_v[i],
         Bep_v[i], Bepp_v[i],
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
   double cf= 1/(1+r_v[i]/cl);

   double r_Der_n= pow(cf,2)*Dx_ptp1_4th(n_v[i+3],n_v[i+2],n_v[i+1],n_v[i],n_v[i-1],dx);
   double r_Der_s= pow(cf,2)*Dx_ptp1_4th(s_v[i+3],s_v[i+2],s_v[i+1],s_v[i],s_v[i-1],dx);
   double r_Der_p= pow(cf,2)*Dx_ptp1_4th(p_v[i+3],p_v[i+2],p_v[i+1],p_v[i],p_v[i-1],dx);
   double r_Der_q= pow(cf,2)*Dx_ptp1_4th(q_v[i+3],q_v[i+2],q_v[i+1],q_v[i],q_v[i-1],dx);

   f_k[i]=
      n_v[i]*(p_v[i]+s_v[i]*q_v[i])
   ;
   p_k[i]= compute_p_k(
      r_v[i],
      n_v[i],   s_v[i],
      p_v[i],   q_v[i],
      V_v[i],   Vp_v[i], 
      Al_v[i],  Alp_v[i],
      Bep_v[i], Bepp_v[i],
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
      cf= 1/(1+r_v[i]/cl);

      r_Der_n= pow(cf,2)*Dx_ptp0_4th(n_v[i+4],n_v[i+3],n_v[i+2],n_v[i+1],n_v[i],dx);
      r_Der_s= pow(cf,2)*Dx_ptp0_4th(s_v[i+4],s_v[i+3],s_v[i+2],s_v[i+1],s_v[i],dx);
      r_Der_p= pow(cf,2)*Dx_ptp0_4th(p_v[i+4],p_v[i+3],p_v[i+2],p_v[i+1],p_v[i],dx);
      r_Der_q= pow(cf,2)*Dx_ptp0_4th(q_v[i+4],q_v[i+3],q_v[i+2],q_v[i+1],q_v[i],dx);

      f_k[i]=
         n_v[i]*(p_v[i]+s_v[i]*q_v[i])
      ;
      p_k[i]= compute_p_k(
         r_v[i],
         n_v[i], s_v[i],
         p_v[i], q_v[i],
         V_v[i],   Vp_v[i], 
         Al_v[i],  Alp_v[i],
         Bep_v[i], Bepp_v[i],
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
         r_v[i],
         n_v[i], s_v[i],
         p_v[i], q_v[i],
         V_v[i], Vp_v[i], 
         Al_v[i],  Alp_v[i],
         Bep_v[i], Bepp_v[i],
         r_Der_n, r_Der_s,
         r_Der_p,  r_Der_q
      );
   }
/*---------------------------------------------------------------------------*/
   i= nx-2;
   cf= 1/(1+r_v[i]/cl);

   r_Der_n= pow(cf,2)*Dx_ptm1_4th(n_v[i+1],n_v[i],n_v[i-1],n_v[i-2],n_v[i-3],dx);
   r_Der_s= pow(cf,2)*Dx_ptm1_4th(s_v[i+1],s_v[i],s_v[i-1],s_v[i-2],s_v[i-3],dx);
   r_Der_p= pow(cf,2)*Dx_ptm1_4th(p_v[i+1],p_v[i],p_v[i-1],p_v[i-2],p_v[i-3],dx);
   r_Der_q= pow(cf,2)*Dx_ptm1_4th(q_v[i+1],q_v[i],q_v[i-1],q_v[i-2],q_v[i-3],dx);

   f_k[i]=
      n_v[i]*(p_v[i]+s_v[i]*q_v[i])
   ;

   p_k[i]= compute_p_k(
      r_v[i],
      n_v[i], s_v[i],
      p_v[i], q_v[i],
      V_v[i],   Vp_v[i], 
      Al_v[i],  Alp_v[i],
      Bep_v[i], Bepp_v[i],
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
   for (int i=0; i<nx; ++i) {
      f_k[i]= dt*f_k[i];
      p_k[i]= dt*p_k[i];
      q_k[i]= dt*q_k[i];
   }
   S_free_k= dt*S_free_k;
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
   compute_potentials(f.n);

   compute_fpqS_ki(exc_i,
      n.n, s.n, p.n, q.n,
      f_k1, p_k1, q_k1, S_free_k1
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
      p.inter_2, q.inter_2,
      n.inter_2, s.inter_2
   );
/*---------------------------------------------------------------------------*/
   compute_potentials(f.inter_2);

   compute_fpqS_ki(exc_i,
      n.inter_2, s.inter_2, p.inter_2, q.inter_2,
      f_k2, p_k2, q_k2, S_free_k2 
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
      p.inter_3, q.inter_3,
      n.inter_3, s.inter_3
   );
/*---------------------------------------------------------------------------*/
   compute_potentials(f.inter_3);

   compute_fpqS_ki(exc_i,
      n.inter_3, s.inter_3, p.inter_3, q.inter_3,
      f_k3, p_k3, q_k3, S_free_k3
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
      p.inter_4, q.inter_4,
      n.inter_4, s.inter_4
   );
/*---------------------------------------------------------------------------*/
   compute_potentials(f.inter_4);

   compute_fpqS_ki(exc_i,
      n.inter_4, s.inter_4, p.inter_4, q.inter_4,
      f_k4, p_k4, q_k4, S_free_k4 
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
      p.np1, q.np1,
      n.np1, s.np1
   );
/*---------------------------------------------------------------------------*/
   compute_potentials(f.np1);

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
   const vector<double> &p_v, const vector<double> &q_v,
   vector<double> &ingoing,
   vector<double> &outgoing)
{
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
	const vector<double> &p_v, const vector<double> &q_v,
	vector<double> &eom_rr)
{
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
      double cf= 1/(1+r_v[i]/cl);

      double r_Der_N= pow(cf,2)*Dx_ptc_4th(n_v[i+2],n_v[i+1],n_v[i-1],n_v[i-2],dx);
      double r_Der_S= pow(cf,2)*Dx_ptc_4th(s_v[i+2],s_v[i+1],s_v[i-1],s_v[i-2],dx);
      double r_Der_P= pow(cf,2)*Dx_ptc_4th(p_v[i+2],p_v[i+1],p_v[i-1],p_v[i-2],dx);
      double r_Der_Q= pow(cf,2)*Dx_ptc_4th(q_v[i+2],q_v[i+1],q_v[i-1],q_v[i-2],dx);

      double t_Der_S= compute_S_free_k(
         r_v[i],
         n_v[i], s_v[i],
         p_v[i], q_v[i],
         V_v[i], Vp_v[i], 
         Al_v[i], Alp_v[i],
         Bep_v[i], Bepp_v[i],
         r_Der_N, r_Der_S,
         r_Der_P, r_Der_Q
      );

      ncc[i]= 
         (2*t_Der_S*n_v[i])/r_v[i] 
      +  r_Der_N*((2*n_v[i])/r_v[i] 
      -  (4*n_v[i]*s_v[i])/r_v[i] + (2*n_v[i]*pow(s_v[i],2))/r_v[i])
      ;
   }
   return;
}
