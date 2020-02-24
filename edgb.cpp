#include <cmath>
using std::pow;
using std::exp;
using std::fabs;
using std::isnan;
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

  ze_free_k{0},
  ze_free_k1{0},
  ze_free_k2{0},
  ze_free_k3{0},
  ze_free_k4{0},

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
  
  W_v(nx,0),
  Wp_v(nx,0),
  Wpp_v(nx,0)
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
	double val= vec[nx-1];
	for (int i=0; i<nx; ++i) {
		vec[i]/= val;
	}
	return;
}
/*===========================================================================*/
void EdGB::compute_scalar_potentials(const vector<double> &f_v)
{
	for (int i=0; i<nx; ++i) {
		V_v[i]= 
			pow(mu,2)*pow(f_v[i],2)
		+	(2.)*la*pow(f_v[i],4)
		;
		Vp_v[i]= 
			(2.)*pow(mu,2)*f_v[i]
		+	(8.)*la*pow(f_v[i],3)
		;
		W_v[i]=
			gbc1*f_v[i]
		+	(1./4.)*gbc2*pow(f_v[i],2)
		;
		Wp_v[i]=
			gbc1
		+	(1./2.)*gbc2*f_v[i]
		;
		Wpp_v[i]=
			(1./2.)*gbc2
		;
	}
	return;
}
/*===========================================================================*/
double EdGB::compute_res_al(
	const double r,
	const double Al, const double Ze,
	const double P,  const double Q,
	const double V,  const double Wp, const double Wpp,
	const double r_Der_Al,
	const double r_Der_P,
	const double r_Der_Q)
{
	return
		-5*Wp*Al*pow(P,2)*Q + Wp*Al*pow(Q,3) + 2*Wp*Al*Q*V + (r*Al*P*Q)/(2.*Ze) - (4*Wp*Al*P*pow(Q,2))/Ze + (4*Wpp*Al*P*Q*Ze)/r - (32*Wp*Wpp*Al*P*pow(Q,2)*Ze)/pow(r,2) - (2*Wp*Al*Q*pow(Ze,2))/pow(r,2) + (16*r_Der_Q*pow(Wp,2)*Al*Q*pow(Ze,2))/pow(r,2) - (48*Wp*Wpp*Al*pow(P,2)*Q*pow(Ze,2))/pow(r,2) + (16*Wp*Wpp*Al*pow(Q,3)*pow(Ze,2))/pow(r,2) + r_Der_P*((4*Wp*Al*Ze)/r - (32*pow(Wp,2)*Al*Q*Ze)/pow(r,2) - (48*pow(Wp,2)*Al*P*pow(Ze,2))/pow(r,2)) + r_Der_Al*(1 - (16*Wp*Q)/r + (64*pow(Wp,2)*pow(Q,2))/pow(r,2) - (20*Wp*P*Ze)/r + (160*pow(Wp,2)*P*Q*Ze)/pow(r,2) + (96*pow(Wp,2)*pow(P,2)*pow(Ze,2))/pow(r,2))
	;
}
/*===========================================================================*/
double EdGB::compute_res_ze(
	const double r,
	const double Ze,
	const double P,  const double Q,
	const double V,  const double Wp, const double Wpp,
	const double r_Der_Ze,
	const double r_Der_P,
	const double r_Der_Q)
{
	return
		2*Wp*pow(P,3) - (r*P*Q)/2. + 6*Wp*P*pow(Q,2) + 4*Wp*P*V - (r*pow(P,2))/(4.*Ze) + (2*Wp*pow(P,2)*Q)/Ze - (r*pow(Q,2))/(4.*Ze) + (2*Wp*pow(Q,3))/Ze - (r*V)/(2.*Ze) + (4*Wp*Q*V)/Ze + Ze/(2.*r) - (4*Wp*Q*Ze)/pow(r,2) + 5*Wp*pow(P,2)*Q*Ze - (4*Wpp*pow(Q,2)*Ze)/r - Wp*pow(Q,3)*Ze + (32*Wp*Wpp*pow(Q,3)*Ze)/pow(r,2) - 2*Wp*Q*V*Ze - (4*Wp*P*pow(Ze,2))/pow(r,2) - (4*Wpp*P*Q*pow(Ze,2))/r + (64*Wp*Wpp*P*pow(Q,2)*pow(Ze,2))/pow(r,2) + (2*Wp*Q*pow(Ze,3))/pow(r,2) + (48*Wp*Wpp*pow(P,2)*Q*pow(Ze,3))/pow(r,2) - (16*Wp*Wpp*pow(Q,3)*pow(Ze,3))/pow(r,2) + r_Der_Ze*(1 - (16*Wp*Q)/r + (64*pow(Wp,2)*pow(Q,2))/pow(r,2) - (20*Wp*P*Ze)/r + (160*pow(Wp,2)*P*Q*Ze)/pow(r,2) + (96*pow(Wp,2)*pow(P,2)*pow(Ze,2))/pow(r,2)) + r_Der_P*((-4*Wp*pow(Ze,2))/r + (32*pow(Wp,2)*Q*pow(Ze,2))/pow(r,2) + (48*pow(Wp,2)*P*pow(Ze,3))/pow(r,2)) + r_Der_Q*((-4*Wp*Ze)/r + (32*pow(Wp,2)*Q*Ze)/pow(r,2) + (32*pow(Wp,2)*P*pow(Ze,2))/pow(r,2) - (16*pow(Wp,2)*Q*pow(Ze,3))/pow(r,2))
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
	vector<double> &al_v,
	vector<double> &ze_v)
{
	compute_scalar_potentials(f_v);

	const double err_tolerance= 1e-12;
	double err= 0;
	do {
		for (int i=exc_i; i<nx-2; ++i) {
			double x=  dx*(2*i+1)/2;

			double cf= 1-x/cl;
			double r=  x/cf;
			double dr= dx*pow(cf,-2);

			double P=  ( p_v[i+1]+ p_v[i])/2;
			double Q=  ( q_v[i+1]+ q_v[i])/2;
			double Al= (al_v[i+1]+al_v[i])/2;
			double Ze= (ze_v[i+1]+ze_v[i])/2;

			double V=   (  V_v[i+1]+  V_v[i])/2;
			double Wp=  ( Wp_v[i+1]+ Wp_v[i])/2;
			double Wpp= (Wpp_v[i+1]+Wpp_v[i])/2;

			double r_Der_P=  ( p_v[i+1]- p_v[i])/dr;
			double r_Der_Q=  ( q_v[i+1]- q_v[i])/dr;
			double r_Der_Al= (al_v[i+1]-al_v[i])/dr;
			double r_Der_Ze= (ze_v[i+1]-ze_v[i])/dr;
/*---------------------------------------------------------------------------*/
			double res_al= compute_res_al(
				r,
				Al, Ze, P,  Q,
				V,  Wp, Wpp,
				r_Der_Al, r_Der_P, r_Der_Q
			);
			double jac_al= 
				(-5*Wp*pow(P,2)*Q)/2. + (Wp*pow(Q,3))/2. + Wp*Q*V + (r*P*Q)/(4.*Ze) - (2*Wp*P*pow(Q,2))/Ze + (2*Wpp*P*Q*Ze)/r - (16*Wp*Wpp*P*pow(Q,2)*Ze)/pow(r,2) - (Wp*Q*pow(Ze,2))/pow(r,2) + (8*r_Der_Q*pow(Wp,2)*Q*pow(Ze,2))/pow(r,2) - (24*Wp*Wpp*pow(P,2)*Q*pow(Ze,2))/pow(r,2) + (8*Wp*Wpp*pow(Q,3)*pow(Ze,2))/pow(r,2) + r_Der_P*((2*Wp*Ze)/r - (16*pow(Wp,2)*Q*Ze)/pow(r,2) - (24*pow(Wp,2)*P*pow(Ze,2))/pow(r,2)) + (1 - (16*Wp*Q)/r + (64*pow(Wp,2)*pow(Q,2))/pow(r,2) - (20*Wp*P*Ze)/r + (160*pow(Wp,2)*P*Q*Ze)/pow(r,2) + (96*pow(Wp,2)*pow(P,2)*pow(Ze,2))/pow(r,2))/dr
			;
/*---------------------------------------------------------------------------*/
			double res_ze= compute_res_ze(
				r,
				Ze, P,  Q,
				V,  Wp, Wpp,
				r_Der_Ze, r_Der_P, r_Der_Q
			);
			double jac_ze= 
				1/(4.*r) - (2*Wp*Q)/pow(r,2) + (5*Wp*pow(P,2)*Q)/2. - (2*Wpp*pow(Q,2))/r - (Wp*pow(Q,3))/2. + (16*Wp*Wpp*pow(Q,3))/pow(r,2) - Wp*Q*V + (r*pow(P,2))/(8.*pow(Ze,2)) - (Wp*pow(P,2)*Q)/pow(Ze,2) + (r*pow(Q,2))/(8.*pow(Ze,2)) - (Wp*pow(Q,3))/pow(Ze,2) + (r*V)/(4.*pow(Ze,2)) - (2*Wp*Q*V)/pow(Ze,2) - (4*Wp*P*Ze)/pow(r,2) - (4*Wpp*P*Q*Ze)/r + (64*Wp*Wpp*P*pow(Q,2)*Ze)/pow(r,2) + (3*Wp*Q*pow(Ze,2))/pow(r,2) + (72*Wp*Wpp*pow(P,2)*Q*pow(Ze,2))/pow(r,2) - (24*Wp*Wpp*pow(Q,3)*pow(Ze,2))/pow(r,2) + r_Der_Ze*((-10*Wp*P)/r + (80*pow(Wp,2)*P*Q)/pow(r,2) + (96*pow(Wp,2)*pow(P,2)*Ze)/pow(r,2)) + r_Der_P*((-4*Wp*Ze)/r + (32*pow(Wp,2)*Q*Ze)/pow(r,2) + (72*pow(Wp,2)*P*pow(Ze,2))/pow(r,2)) + (1 - (16*Wp*Q)/r + (64*pow(Wp,2)*pow(Q,2))/pow(r,2) - (20*Wp*P*Ze)/r + (160*pow(Wp,2)*P*Q*Ze)/pow(r,2) + (96*pow(Wp,2)*pow(P,2)*pow(Ze,2))/pow(r,2))/dr + r_Der_Q*((-2*Wp)/r + (16*pow(Wp,2)*Q)/pow(r,2) + (32*pow(Wp,2)*P*Ze)/pow(r,2) - (24*pow(Wp,2)*Q*pow(Ze,2))/pow(r,2))
			;
/*---------------------------------------------------------------------------*/
			al_v[i+1]-= res_al/jac_al;
			ze_v[i+1]-= res_ze/jac_ze;

//			cout<<setw(10)<<i<<"\t";
//			cout<<setw(10)<<Q<<"\t";
//			cout<<setw(10)<<ze_v[i+1]<<"\t";
//			cout<<setw(10)<<res_ze<<"\t";
//			cout<<setw(10)<<jac_ze<<endl;

			err+= fabs(res_al);
			err+= fabs(res_ze);
		}
		al_v[nx-1]= al_v[nx-2];
		ze_v[nx-1]= 0;
		rescale(al_v);
		err/= nx;
	} while (err>err_tolerance);
	return;
}
/*===========================================================================*/
/* public interface to solve for metric fields (for initial data) */
/*===========================================================================*/
void EdGB::solve_metric_fields(
	const int exc_i,
	const Field &phi_f, const Field &phi_p, const Field &phi_q,
	Field &al, Field &ze) 
{
	solve_for_metric_relaxation(exc_i,
		phi_f.n,phi_p.n,phi_q.n,
		al.n,ze.n
	);
	for (int i=0; i<nx; ++i) {
		al.np1[i]= al.n[i];
		ze.np1[i]= ze.n[i];
	}	
	return;
}
/*===========================================================================*/
double EdGB::compute_ze_free_k(
	const double r,
	const double Al, const double Ze,
	const double P,  const double Q,
	const double V,  const double Vp, const double Wp, const double Wpp,
	const double r_Der_Al,
	const double r_Der_Ze,
	const double r_Der_P,
	const double r_Der_Q)
{
	double Zer= Ze/r;
	double Qr= Q/r;
	double val=
		(r*Al*pow(P,2))/4. - 6*Wp*Al*pow(P,2)*Q + (r*Al*pow(Q,2))/4. - 2*Wp*Al*pow(Q,3) - (r*Al*V)/2. + 4*Wp*Al*Q*V + (r*Al*P*Q)/(2.*Ze) - (4*Wp*Al*P*pow(Q,2))/Ze - 2*Wp*Al*pow(P,3)*Ze + (4*Wpp*Al*P*Q*Ze)/r - 2*Wp*Al*P*pow(Q,2)*Ze - (32*Wp*Wpp*Al*P*pow(Q,2)*Ze)/pow(r,2) + 4*Wp*Al*P*V*Ze + (Al*pow(Ze,2))/(2.*r) - (4*Vp*Wp*Al*pow(Ze,2))/r + (4*Wpp*Al*pow(P,2)*pow(Ze,2))/r + (4*Wp*Al*Q*pow(Ze,2))/pow(r,2) + (32*Vp*pow(Wp,2)*Al*Q*pow(Ze,2))/pow(r,2) - (64*Wp*Wpp*Al*pow(P,2)*Q*pow(Ze,2))/pow(r,2) - (64*pow(Wp,2)*Al*pow(Q,2)*pow(Ze,2))/pow(r,3) + (4*Wp*Al*P*pow(Ze,3))/pow(r,2) + (32*Vp*pow(Wp,2)*Al*P*pow(Ze,3))/pow(r,2) - (32*Wp*Wpp*Al*pow(P,3)*pow(Ze,3))/pow(r,2) - (128*pow(Wp,2)*Al*P*Q*pow(Ze,3))/pow(r,3) - (80*pow(Wp,2)*Al*pow(P,2)*pow(Ze,4))/pow(r,3) + (16*pow(Wp,2)*Al*pow(Q,2)*pow(Ze,4))/pow(r,3) + (32*pow(Wp,2)*Al*V*pow(Ze,4))/pow(r,3) + r_Der_P*((4*Wp*Al*Ze)/r - (32*pow(Wp,2)*Al*Q*Ze)/pow(r,2) - (32*pow(Wp,2)*Al*P*pow(Ze,2))/pow(r,2)) + r_Der_Q*((4*Wp*Al*pow(Ze,2))/r - (32*pow(Wp,2)*Al*Q*pow(Ze,2))/pow(r,2) - (32*pow(Wp,2)*Al*P*pow(Ze,3))/pow(r,2)) + pow(r_Der_Ze,2)*((128*pow(Wp,2)*Al*pow(Ze,4))/pow(r,3) - (768*pow(Wp,3)*Al*Q*pow(Ze,4))/pow(r,4) - (768*pow(Wp,3)*Al*P*pow(Ze,5))/pow(r,4)) + r_Der_Ze*(Al*Ze - (12*Wp*Al*Q*Ze)/r + (32*pow(Wp,2)*Al*pow(Q,2)*Ze)/pow(r,2) - (12*Wp*Al*P*pow(Ze,2))/r + (96*pow(Wp,2)*Al*P*Q*pow(Ze,2))/pow(r,2) + (48*pow(Wp,2)*Al*pow(P,2)*pow(Ze,3))/pow(r,2) + (16*pow(Wp,2)*Al*pow(Q,2)*pow(Ze,3))/pow(r,2) - (32*pow(Wp,2)*Al*V*pow(Ze,3))/pow(r,2) + (256*r_Der_P*pow(Wp,3)*Al*pow(Ze,4))/pow(r,4) + (256*pow(Wp,2)*Wpp*Al*P*Q*pow(Ze,4))/pow(r,4) - (32*pow(Wp,2)*Al*pow(Ze,5))/pow(r,4) + (256*r_Der_Q*pow(Wp,3)*Al*pow(Ze,5))/pow(r,4) + (256*pow(Wp,2)*Wpp*Al*pow(Q,2)*pow(Ze,5))/pow(r,4)) + pow(r_Der_Al,2)*((-256*pow(Wp,3)*P*pow(Ze,5))/(pow(r,4)*Al) - (256*pow(Wp,3)*Q*pow(Ze,6))/(pow(r,4)*Al)) + r_Der_Al*((4*Wp*Q*pow(Ze,2))/r - (32*pow(Wp,2)*pow(Q,2)*pow(Ze,2))/pow(r,2) + (4*Wp*P*pow(Ze,3))/r - (32*pow(Wp,2)*P*Q*pow(Ze,3))/pow(r,2) + (32*pow(Wp,2)*pow(Ze,4))/pow(r,4) - (256*r_Der_Q*pow(Wp,3)*pow(Ze,4))/pow(r,4) - (16*pow(Wp,2)*pow(P,2)*pow(Ze,4))/pow(r,2) + (16*pow(Wp,2)*pow(Q,2)*pow(Ze,4))/pow(r,2) - (256*pow(Wp,2)*Wpp*pow(Q,2)*pow(Ze,4))/pow(r,4) - (32*pow(Wp,2)*V*pow(Ze,4))/pow(r,2) - (256*r_Der_P*pow(Wp,3)*pow(Ze,5))/pow(r,4) - (256*pow(Wp,2)*Wpp*P*Q*pow(Ze,5))/pow(r,4) + r_Der_Ze*((-64*pow(Wp,2)*pow(Ze,3))/pow(r,3) + (512*pow(Wp,3)*Q*pow(Ze,3))/pow(r,4) + (256*pow(Wp,3)*P*pow(Ze,4))/pow(r,4) + (128*pow(Wp,2)*pow(Ze,5))/pow(r,3) - (1024*pow(Wp,3)*Q*pow(Ze,5))/pow(r,4) - (768*pow(Wp,3)*P*pow(Ze,6))/pow(r,4)))
	;
	val/=
		1 - 16*Qr*Wp + 64*pow(Qr,2)*pow(Wp,2) - 32*pow(Wp,2)*pow(Zer,4) + 256*r_Der_Q*pow(Wp,3)*pow(Zer,4) + 256*pow(Qr,2)*pow(r,2)*pow(Wp,2)*Wpp*pow(Zer,4) - 16*Wp*Zer*P + 128*Qr*pow(Wp,2)*Zer*P + 64*pow(Wp,2)*pow(Zer,2)*pow(P,2) - 128*r_Der_Ze*(-(pow(Wp,2)*pow(Zer,3)) + 8*Qr*pow(Wp,3)*pow(Zer,3) + 6*pow(Wp,3)*pow(Zer,4)*P) - (128*r_Der_Al*(-(r*pow(Wp,2)*pow(Zer,4)) + 8*Qr*r*pow(Wp,3)*pow(Zer,4) + 6*r*pow(Wp,3)*pow(Zer,5)*P))/Al
	;
	return val;
}
/*===========================================================================*/
double EdGB::compute_p_k(
	const double r,
	const double Al, const double Ze,
	const double P,  const double Q,
	const double V,  const double Vp, const double Wp, const double Wpp,
	const double r_Der_Al,
	const double r_Der_Ze,
	const double r_Der_P,
	const double r_Der_Q)
{
	double Qr= Q/r;
	double Zer= Ze/r;

	double p_k= 
		2*Qr*Al - Vp*Al - 32*pow(Qr,2)*Wp*Al + 16*Qr*Vp*Wp*Al + 128*pow(Qr,3)*pow(Wp,2)*Al - 64*pow(Qr,2)*Vp*pow(Wp,2)*Al + 4*Wp*pow(Zer,4)*Al - 32*pow(Qr,2)*pow(r,2)*Wp*Wpp*pow(Zer,4)*Al + 2*Zer*Al*P - 64*Qr*Wp*Zer*Al*P + 16*Vp*Wp*Zer*Al*P + 384*pow(Qr,2)*pow(Wp,2)*Zer*Al*P - 128*Qr*Vp*pow(Wp,2)*Zer*Al*P + 32*Qr*Wp*Wpp*pow(Zer,3)*Al*P - 256*pow(Qr,3)*pow(r,2)*Wp*pow(Wpp,2)*pow(Zer,3)*Al*P - 32*Wp*pow(Zer,2)*Al*pow(P,2) + 384*Qr*pow(Wp,2)*pow(Zer,2)*Al*pow(P,2) - 64*Vp*pow(Wp,2)*pow(Zer,2)*Al*pow(P,2) + 32*Wp*Wpp*pow(Zer,4)*Al*pow(P,2) - 256*pow(Qr,2)*pow(r,2)*Wp*pow(Wpp,2)*pow(Zer,4)*Al*pow(P,2) + 128*pow(Wp,2)*pow(Zer,3)*Al*pow(P,3) + 64*pow(r_Der_Ze,2)*(-(Qr*pow(Wp,2)*pow(Zer,2)*Al) + 8*pow(Qr,2)*pow(Wp,3)*pow(Zer,2)*Al + 4*Qr*pow(Wp,3)*pow(Zer,3)*Al*P) + (64*pow(r_Der_Al,2)*(-(Qr*pow(r,2)*pow(Wp,2)*pow(Zer,4)) + 8*pow(Qr,2)*pow(r,2)*pow(Wp,3)*pow(Zer,4) - pow(Wp,2)*pow(Zer,3)*P + 8*Qr*pow(Wp,3)*pow(Zer,3)*P + 8*Qr*pow(r,2)*pow(Wp,3)*pow(Zer,5)*P + 8*pow(Wp,3)*pow(Zer,4)*pow(P,2)))/Al + r_Der_P*(r*Zer*Al - 16*Qr*r*Wp*Zer*Al + 64*pow(Qr,2)*r*pow(Wp,2)*Zer*Al + (32*pow(Wp,2)*pow(Zer,3)*Al)/r - 256*pow(Qr,2)*r*pow(Wp,2)*Wpp*pow(Zer,3)*Al - 32*r*pow(Wp,2)*pow(Zer,5)*Al + 256*pow(Qr,2)*pow(r,3)*pow(Wp,2)*Wpp*pow(Zer,5)*Al + (256*r_Der_Q*(-(pow(Wp,3)*pow(Zer,3)*Al) + pow(r,2)*pow(Wp,3)*pow(Zer,5)*Al))/r - 16*r*Wp*pow(Zer,2)*Al*P + 128*Qr*r*pow(Wp,2)*pow(Zer,2)*Al*P + 64*r*pow(Wp,2)*pow(Zer,3)*Al*pow(P,2)) - 4*(-(Wp*pow(Zer,2)*Al) + 8*Qr*pow(Wp,2)*pow(Zer,2)*Al + 8*pow(Wp,2)*pow(Zer,3)*Al*P)*(pow(Qr,2)*pow(r,2) - pow(P,2) + 2*V) - 2*(-(Wp*Zer) + 8*pow(Qr,2)*pow(r,2)*Wp*Wpp*Zer)*Al*(pow(Qr,2)*pow(r,2)*Zer + 2*Qr*P + Zer*pow(P,2) - 2*Zer*V) + r_Der_Q*(Al - 16*Qr*Wp*Al + 64*pow(Qr,2)*pow(Wp,2)*Al - 32*pow(Wp,2)*pow(Zer,4)*Al - 16*Wp*Zer*Al*P + 128*Qr*pow(Wp,2)*Zer*Al*P - 256*Qr*pow(Wp,2)*Wpp*pow(Zer,3)*Al*P + 64*pow(Wp,2)*pow(Zer,2)*Al*pow(P,2) - 256*pow(Wp,2)*Wpp*pow(Zer,4)*Al*pow(P,2) - 16*pow(Wp,2)*Zer*Al*(pow(Qr,2)*pow(r,2)*Zer + 2*Qr*P + Zer*pow(P,2) - 2*Zer*V)) + r_Der_Ze*(-16*Wp*pow(Zer,3)*Al + 160*Qr*pow(Wp,2)*pow(Zer,3)*Al - 256*Qr*r_Der_Q*pow(Wp,3)*pow(Zer,3)*Al - 256*pow(Qr,3)*pow(r,2)*pow(Wp,2)*Wpp*pow(Zer,3)*Al + Al*P - 16*Qr*Wp*Al*P + 64*pow(Qr,2)*pow(Wp,2)*Al*P - 64*Qr*Wp*Wpp*pow(Zer,2)*Al*P + 512*pow(Qr,2)*pow(Wp,2)*Wpp*pow(Zer,2)*Al*P + 96*pow(Wp,2)*pow(Zer,4)*Al*P - 16*Wp*Zer*Al*pow(P,2) + 128*Qr*pow(Wp,2)*Zer*Al*pow(P,2) - 128*Wp*Wpp*pow(Zer,3)*Al*pow(P,2) + 1280*Qr*pow(Wp,2)*Wpp*pow(Zer,3)*Al*pow(P,2) + 64*pow(Wp,2)*pow(Zer,2)*Al*pow(P,3) + 768*pow(Wp,2)*Wpp*pow(Zer,4)*Al*pow(P,3) - (64*r_Der_P*(pow(Wp,2)*pow(Zer,2)*Al - 8*Qr*pow(Wp,3)*pow(Zer,2)*Al - 2*pow(r,2)*pow(Wp,2)*pow(Zer,4)*Al + 16*Qr*pow(r,2)*pow(Wp,3)*pow(Zer,4)*Al - 4*pow(Wp,3)*pow(Zer,3)*Al*P + 12*pow(r,2)*pow(Wp,3)*pow(Zer,5)*Al*P))/r + 4*Al*(-Wp + 8*Qr*pow(Wp,2) + 4*pow(Wp,2)*Zer*P)*(pow(Qr,2)*pow(r,2)*Zer + 2*Qr*P + Zer*pow(P,2) - 2*Zer*V)) + r_Der_Al*(Qr*r - 16*pow(Qr,2)*r*Wp + 64*pow(Qr,3)*r*pow(Wp,2) + (8*Wp*pow(Zer,2))/r - (64*Qr*pow(Wp,2)*pow(Zer,2))/r - 64*pow(Qr,2)*r*Wp*Wpp*pow(Zer,2) + 512*pow(Qr,3)*r*pow(Wp,2)*Wpp*pow(Zer,2) - 16*r*Wp*pow(Zer,4) + 128*Qr*r*pow(Wp,2)*pow(Zer,4) + r*Zer*P - 32*Qr*r*Wp*Zer*P + 192*pow(Qr,2)*r*pow(Wp,2)*Zer*P - (64*pow(Wp,2)*pow(Zer,3)*P)/r - 192*Qr*r*Wp*Wpp*pow(Zer,3)*P + 2048*pow(Qr,2)*r*pow(Wp,2)*Wpp*pow(Zer,3)*P + 96*r*pow(Wp,2)*pow(Zer,5)*P - 16*r*Wp*pow(Zer,2)*pow(P,2) + 192*Qr*r*pow(Wp,2)*pow(Zer,2)*pow(P,2) - 128*r*Wp*Wpp*pow(Zer,4)*pow(P,2) + 2304*Qr*r*pow(Wp,2)*Wpp*pow(Zer,4)*pow(P,2) + 64*r*pow(Wp,2)*pow(Zer,3)*pow(P,3) + 768*r*pow(Wp,2)*Wpp*pow(Zer,5)*pow(P,3) + (64*r_Der_Q*(-(pow(Wp,2)*pow(Zer,2)) + 8*Qr*pow(Wp,3)*pow(Zer,2) + 8*pow(Wp,3)*pow(Zer,3)*P))/r - 64*r_Der_P*(3*pow(Wp,2)*pow(Zer,3) - 24*Qr*pow(Wp,3)*pow(Zer,3) - 2*pow(r,2)*pow(Wp,2)*pow(Zer,5) + 16*Qr*pow(r,2)*pow(Wp,3)*pow(Zer,5) - 20*pow(Wp,3)*pow(Zer,4)*P + 12*pow(r,2)*pow(Wp,3)*pow(Zer,6)*P) + (16*r_Der_Ze*(-(Wp*Zer) + 16*Qr*pow(Wp,2)*Zer - 64*pow(Qr,2)*pow(Wp,3)*Zer - 8*Qr*pow(r,2)*pow(Wp,2)*pow(Zer,3) + 64*pow(Qr,2)*pow(r,2)*pow(Wp,3)*pow(Zer,3) + 12*pow(Wp,2)*pow(Zer,2)*P - 96*Qr*pow(Wp,3)*pow(Zer,2)*P + 48*Qr*pow(r,2)*pow(Wp,3)*pow(Zer,4)*P - 32*pow(Wp,3)*pow(Zer,3)*pow(P,2)))/r + 4*r*(-(Wp*Zer) + 8*Qr*pow(Wp,2)*Zer + 4*pow(Wp,2)*pow(Zer,2)*P)*(pow(Qr,2)*pow(r,2)*Zer + 2*Qr*P + Zer*pow(P,2) - 2*Zer*V))
	;
	p_k/=
		1 - 16*Qr*Wp + 64*pow(Qr,2)*pow(Wp,2) - 32*pow(Wp,2)*pow(Zer,4) + 256*r_Der_Q*pow(Wp,3)*pow(Zer,4) + 256*pow(Qr,2)*pow(r,2)*pow(Wp,2)*Wpp*pow(Zer,4) - 16*Wp*Zer*P + 128*Qr*pow(Wp,2)*Zer*P + 64*pow(Wp,2)*pow(Zer,2)*pow(P,2) - 128*r_Der_Ze*(-(pow(Wp,2)*pow(Zer,3)) + 8*Qr*pow(Wp,3)*pow(Zer,3) + 6*pow(Wp,3)*pow(Zer,4)*P) - (128*r_Der_Al*(-(r*pow(Wp,2)*pow(Zer,4)) + 8*Qr*r*pow(Wp,3)*pow(Zer,4) + 6*r*pow(Wp,3)*pow(Zer,5)*P))/Al
	;
	return p_k;
}
/*===========================================================================*/
void EdGB::compute_fpq_ki(
	const int val,
	const int exc_i,
	const vector<double> &al_v,
	const vector<double> &ze_v,
	const vector<double> &f_v, 
	const vector<double> &p_v, 
	const vector<double> &q_v)
{
	compute_scalar_potentials(f_v);

	double r, cf, al, ze, p, q, r_Der_al, r_Der_ze, r_Der_p, r_Der_q;
/*---------------------------------------------------------------------------*/
	for (int i=exc_i+2; i<nx-2; ++i) {
		r= r_v[i];
		cf= 1/(1+r/cl);

		al= al_v[i];
		ze= ze_v[i];
		p= p_v[i];
		q= q_v[i];

		r_Der_al= pow(cf,2)*Dx_ptc_4th(al_v[i+2],al_v[i+1],al_v[i-1],al_v[i-2],dx);
		r_Der_ze= pow(cf,2)*Dx_ptc_4th(ze_v[i+2],ze_v[i+1],ze_v[i-1],ze_v[i-2],dx);
		r_Der_p=  pow(cf,2)*Dx_ptc_4th( p_v[i+2], p_v[i+1], p_v[i-1], p_v[i-2],dx);
		r_Der_q=  pow(cf,2)*Dx_ptc_4th( q_v[i+2], q_v[i+1], q_v[i-1], q_v[i-2],dx);

		f_k[i]= 
			al*(p+ze*q)
		;
		p_k[i]= compute_p_k(
			r,
			al, ze,
			p,  q,
			V_v[i],  Vp_v[i], Wp_v[i], Wpp_v[i],
			r_Der_al, r_Der_ze,
			r_Der_p,  r_Der_q
		);
		q_k[i]=
			pow(cf,2)*Dx_ptc_4th(
				al_v[i+2]*(p_v[i+2]+ze_v[i+2]*q_v[i+2]),
				al_v[i+1]*(p_v[i+1]+ze_v[i+1]*q_v[i+1]),
				al_v[i-1]*(p_v[i-1]+ze_v[i-1]*q_v[i-1]),
				al_v[i-2]*(p_v[i-2]+ze_v[i-2]*q_v[i-2]),
				dx
			)
		;
	}
/*---------------------------------------------------------------------------*/
	int i= exc_i+1;
	r= r_v[i];
	cf= 1/(1+r/cl);

	al= al_v[i];
	ze= ze_v[i];
	p= p_v[i];
	q= q_v[i];

	r_Der_al= pow(cf,2)*Dx_ptp1_4th(al_v[i+3],al_v[i+2],al_v[i+1],al_v[i],al_v[i-1],dx);
	r_Der_ze= pow(cf,2)*Dx_ptp1_4th(ze_v[i+3],ze_v[i+2],ze_v[i+1],ze_v[i],ze_v[i-1],dx);
	r_Der_p=  pow(cf,2)*Dx_ptp1_4th( p_v[i+3], p_v[i+2], p_v[i+1], p_v[i], p_v[i-1],dx);
	r_Der_q=  pow(cf,2)*Dx_ptp1_4th( q_v[i+3], q_v[i+2], q_v[i+1], q_v[i], q_v[i-1],dx);

	f_k[i]=
		al*(p+ze*q)
	;
	p_k[i]= compute_p_k(
		r,
		al, ze,
		p,  q,
		V_v[i],  Vp_v[i], Wp_v[i], Wpp_v[i],
		r_Der_al, r_Der_ze,
		r_Der_p,  r_Der_q
	);
	q_k[i]=
		pow(cf,2)*Dx_ptp1_4th(
			al_v[i+3]*(p_v[i+3]+ze_v[i+3]*q_v[i+3]),
			al_v[i+2]*(p_v[i+2]+ze_v[i+2]*q_v[i+2]),
			al_v[i+1]*(p_v[i+1]+ze_v[i+1]*q_v[i+1]),
			al_v[i  ]*(p_v[i  ]+ze_v[i  ]*q_v[i  ]),
			al_v[i-1]*(p_v[i-1]+ze_v[i-1]*q_v[i-1]),
			dx
		)
	;
/*---------------------------------------------------------------------------*/
/* if excising then freely evolve p, q, and ze */
/*---------------------------------------------------------------------------*/
	if (exc_i>0) {
		i= exc_i;
		r= r_v[i];
		cf= 1/(1+r/cl);

		al= al_v[i];
		ze= ze_v[i];
		p= p_v[i];
		q= q_v[i];

		r_Der_al= pow(cf,2)*Dx_ptp0_4th(al_v[i+4],al_v[i+3],al_v[i+2],al_v[i+1],al_v[i],dx);
		r_Der_ze= pow(cf,2)*Dx_ptp0_4th(ze_v[i+4],ze_v[i+3],ze_v[i+2],ze_v[i+1],ze_v[i],dx);
		r_Der_p=  pow(cf,2)*Dx_ptp0_4th( p_v[i+4], p_v[i+3], p_v[i+2], p_v[i+1], p_v[i],dx);
		r_Der_q=  pow(cf,2)*Dx_ptp0_4th( q_v[i+4], q_v[i+3], q_v[i+2], q_v[i+1], q_v[i],dx);

		f_k[i]=
			al*(p+ze*q)
		;
		p_k[i]= compute_p_k(
			r,
			al, ze,
			p,  q,
			V_v[i],  Vp_v[i], Wp_v[i], Wpp_v[i],
			r_Der_al, r_Der_ze,
			r_Der_p,  r_Der_q
		);
		q_k[i]=
			pow(cf,2)*Dx_ptp0_4th(
				al_v[i+4]*(p_v[i+4]+ze_v[i+4]*q_v[i+4]),
				al_v[i+3]*(p_v[i+3]+ze_v[i+3]*q_v[i+3]),
				al_v[i+2]*(p_v[i+2]+ze_v[i+2]*q_v[i+2]),
				al_v[i+1]*(p_v[i+1]+ze_v[i+1]*q_v[i+1]),
				al_v[i  ]*(p_v[i]  +ze_v[i  ]*q_v[i  ]),
				dx
			)
		;
		ze_free_k= compute_ze_free_k(
			r,
			al, ze,
			p,  q,
			V_v[i],  Vp_v[i], Wp_v[i], Wpp_v[i],
			r_Der_al, r_Der_ze,
			r_Der_p,  r_Der_q
		);
	}
/*---------------------------------------------------------------------------*/
	i= nx-2;

	r= r_v[i];
	cf= 1/(1+r/cl);

	al= al_v[i];
	ze= ze_v[i];
	p= p_v[i];
	q= q_v[i];

	r_Der_al= pow(cf,2)*Dx_ptm1_4th(al_v[i+1],al_v[i],al_v[i-1],al_v[i-2],al_v[i-3],dx);
	r_Der_ze= pow(cf,2)*Dx_ptm1_4th(ze_v[i+1],ze_v[i],ze_v[i-1],ze_v[i-2],ze_v[i-3],dx);
	r_Der_p=  pow(cf,2)*Dx_ptm1_4th( p_v[i+1], p_v[i], p_v[i-1], p_v[i-2], p_v[i-3],dx);
	r_Der_q=  pow(cf,2)*Dx_ptm1_4th( q_v[i+1], q_v[i], q_v[i-1], q_v[i-2], q_v[i-3],dx);

	f_k[i]=
		al*(p+ze*q)
	;
	p_k[i]= compute_p_k(
		r,
		al, ze,
		p,  q,
		V_v[i],  Vp_v[i], Wp_v[i], Wpp_v[i],
		r_Der_al, r_Der_ze,
		r_Der_p,  r_Der_q
	);
	q_k[i]=
		pow(cf,2)*Dx_ptm1_4th(
			al_v[i+1]*(p_v[i+1]+ze_v[i+1]*q_v[i+1]),
			al_v[i  ]*(p_v[i  ]+ze_v[i  ]*q_v[i  ]),
			al_v[i-1]*(p_v[i-1]+ze_v[i-1]*q_v[i-1]),
			al_v[i-2]*(p_v[i-2]+ze_v[i-2]*q_v[i-2]),
			al_v[i-3]*(p_v[i-3]+ze_v[i-3]*q_v[i-3]),
			dx
		)
	;
/*---------------------------------------------------------------------------*/
/* regularity condition at spatial infinity: no evolution */
/*---------------------------------------------------------------------------*/
	f_k[nx-1]= 0;
	p_k[nx-1]= 0;
	q_k[nx-1]= 0; 
/*---------------------------------------------------------------------------*/
/* set the k vectors for RK4 method */
/*---------------------------------------------------------------------------*/
	switch(val)
	{
	case 1: 
		for (int i=0; i<nx; ++i) {
			f_k1[i]= dt*f_k[i];
			p_k1[i]= dt*p_k[i];
			q_k1[i]= dt*q_k[i];
		}
		ze_free_k1= dt*ze_free_k;
	break;
	case 2:
		for (int i=0; i<nx; ++i) {
			f_k2[i]= dt*f_k[i];
			p_k2[i]= dt*p_k[i];
			q_k2[i]= dt*q_k[i];
		}
		ze_free_k2= dt*ze_free_k;
	break;
	case 3:
		for (int i=0; i<nx; ++i) {
			f_k3[i]= dt*f_k[i];
			p_k3[i]= dt*p_k[i];
			q_k3[i]= dt*q_k[i];
		}
		ze_free_k3= dt*ze_free_k;
	break;
	case 4:
		for (int i=0; i<nx; ++i) {
			f_k4[i]= dt*f_k[i];
			p_k4[i]= dt*p_k[i];
			q_k4[i]= dt*q_k[i];
		}
		ze_free_k4= dt*ze_free_k;
	break;
	default:
		cout << "ERROR: val = " << val << endl;
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
	Field &al, Field &ze, Field &f, Field &p, Field &q)
{
	KO_filter(exc_i,nx,"even",f.n);
	KO_filter(exc_i,nx,"even",p.n);
	KO_filter(exc_i,nx,"odd" ,q.n);
/*---------------------------------------------------------------------------*/
	 int start_i= (exc_i>0) ? exc_i : 1;
/*---------------------------------------------------------------------------*/
	compute_fpq_ki(1, exc_i,
		al.n, ze.n, f.n, p.n, q.n
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
		ze.inter_2[exc_i]= ze.n[exc_i]+0.5*ze_free_k1;
	}
	solve_for_metric_relaxation(exc_i,
		f.inter_2,  p.inter_2, q.inter_2,
		al.inter_2,ze.inter_2
	);
/*---------------------------------------------------------------------------*/
	compute_fpq_ki(2, exc_i,
		al.inter_2, ze.inter_2, f.inter_2, p.inter_2, q.inter_2
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
		ze.inter_3[exc_i]= ze.n[exc_i]+0.5*ze_free_k2;
	}
	solve_for_metric_relaxation(exc_i,
		f.inter_3,  p.inter_3, q.inter_3,
		al.inter_3,ze.inter_3
	);
/*---------------------------------------------------------------------------*/
	compute_fpq_ki(3, exc_i,
		al.inter_3, ze.inter_3, f.inter_3, p.inter_3, q.inter_3
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
		ze.inter_4[exc_i]= ze.n[exc_i]+0.5*ze_free_k3;
	}
	solve_for_metric_relaxation(exc_i,
		f.inter_4,  p.inter_4, q.inter_4,
		al.inter_4,ze.inter_4
	);
/*---------------------------------------------------------------------------*/
	compute_fpq_ki(4, exc_i,
		al.inter_4, ze.inter_4, f.inter_4, p.inter_4, q.inter_4
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
		ze.np1[exc_i]= ze.n[exc_i]+(ze_free_k1+2.*ze_free_k2+2.*ze_free_k3+ze_free_k4)/6.;
	}
	solve_for_metric_relaxation(exc_i,
		f.np1,  p.np1, q.np1,
		al.np1,ze.np1
	);
/*---------------------------------------------------------------------------*/
	return;
}
/*===========================================================================*/
int EdGB::compute_radial_characteristic(
	const double r,
	const double Al, const double Ze,
	const double P,  const double Q,
	const double V,  const double Vp, const double Wp, const double Wpp,
	const double r_Der_Al,
	const double r_Der_Ze,
	const double r_Der_P,
	const double r_Der_Q)
{
/*---------------------------------------------------------------------------*/
	double Qr=  Q/r;
	double Zer= Ze/r;

	double t_Der_P= compute_p_k(
		r,
		Al, Ze,
		P,  Q,
		V,  Vp, Wp, Wpp,
		r_Der_Al,
		r_Der_Ze,
		r_Der_P,
		r_Der_Q
	);
/*---------------------------------------------------------------------------*/
	double ep_dtp= 
		(768*r_Der_Q*pow(Wp,3)*pow(Zer,4)*(-1 + 8*Qr*Wp + 8*Wp*Zer*P))/(-1 + 8*Qr*Wp + 12*Wp*Zer*P) + (-1 + 24*Qr*Wp - 192*pow(Qr,2)*pow(Wp,2) + 512*pow(Qr,3)*pow(Wp,3) - 32*pow(Qr,2)*pow(r,2)*pow(Wp,2)*pow(Zer,2) + 256*pow(Qr,3)*pow(r,2)*pow(Wp,3)*pow(Zer,2) + 96*pow(Wp,2)*pow(Zer,4) - 768*Qr*pow(Wp,3)*pow(Zer,4) - 768*pow(Qr,2)*pow(r,2)*pow(Wp,2)*Wpp*pow(Zer,4) + 6144*pow(Qr,3)*pow(r,2)*pow(Wp,3)*Wpp*pow(Zer,4) + 28*Wp*Zer*P - 448*Qr*pow(Wp,2)*Zer*P + 1792*pow(Qr,2)*pow(Wp,3)*Zer*P + 192*pow(Qr,2)*pow(r,2)*pow(Wp,3)*pow(Zer,3)*P - 768*pow(Wp,3)*pow(Zer,5)*P + 6144*pow(Qr,2)*pow(r,2)*pow(Wp,3)*Wpp*pow(Zer,5)*P - 288*pow(Wp,2)*pow(Zer,2)*pow(P,2) + 2304*Qr*pow(Wp,3)*pow(Zer,2)*pow(P,2) + 960*pow(Wp,3)*pow(Zer,3)*pow(P,3) - 64*pow(Wp,2)*pow(Zer,2)*V + 512*Qr*pow(Wp,3)*pow(Zer,2)*V + 384*pow(Wp,3)*pow(Zer,3)*P*V)/(-1 + 8*Qr*Wp + 12*Wp*Zer*P)
	;
	double ep_dtq= 
		0
	;
	double ep_drp=
		-1536*r_Der_P*pow(Wp,3)*pow(Zer,4)*Al - (768*r_Der_Q*pow(Wp,3)*pow(Zer,5)*Al*(-r + 4*Qr*r*Wp + 8*r*Wp*Zer*P))/(-1 + 8*Qr*Wp + 12*Wp*Zer*P) - (r*Zer*Al*(1 - 36*Qr*Wp + 480*pow(Qr,2)*pow(Wp,2) - 2816*pow(Qr,3)*pow(Wp,3) + 6144*pow(Qr,4)*pow(Wp,4) + 32*pow(Qr,2)*pow(r,2)*pow(Wp,2)*pow(Zer,2) - 256*pow(Qr,3)*pow(r,2)*pow(Wp,3)*pow(Zer,2) - 96*pow(Wp,2)*pow(Zer,4) + 1152*Qr*pow(Wp,3)*pow(Zer,4) - 3072*pow(Qr,2)*pow(Wp,4)*pow(Zer,4) + 768*pow(Qr,2)*pow(r,2)*pow(Wp,2)*Wpp*pow(Zer,4) - 9216*pow(Qr,3)*pow(r,2)*pow(Wp,3)*Wpp*pow(Zer,4) + 24576*pow(Qr,4)*pow(r,2)*pow(Wp,4)*Wpp*pow(Zer,4) - 36*Wp*Zer*P + 1104*Qr*pow(Wp,2)*Zer*P - 10752*pow(Qr,2)*pow(Wp,3)*Zer*P + 33792*pow(Qr,3)*pow(Wp,4)*Zer*P - 448*pow(Qr,2)*pow(r,2)*pow(Wp,3)*pow(Zer,3)*P + 1280*pow(Qr,3)*pow(r,2)*pow(Wp,4)*pow(Zer,3)*P + 1536*Qr*pow(Wp,2)*Wpp*pow(Zer,3)*P - 24576*pow(Qr,2)*pow(Wp,3)*Wpp*pow(Zer,3)*P + 98304*pow(Qr,3)*pow(Wp,4)*Wpp*pow(Zer,3)*P + 1536*pow(Wp,3)*pow(Zer,5)*P - 9216*Qr*pow(Wp,4)*pow(Zer,5)*P - 12288*pow(Qr,2)*pow(r,2)*pow(Wp,3)*Wpp*pow(Zer,5)*P + 73728*pow(Qr,3)*pow(r,2)*pow(Wp,4)*Wpp*pow(Zer,5)*P + 512*pow(Wp,2)*pow(Zer,2)*pow(P,2) - 11520*Qr*pow(Wp,3)*pow(Zer,2)*pow(P,2) + 59392*pow(Qr,2)*pow(Wp,4)*pow(Zer,2)*pow(P,2) + 1536*pow(Qr,2)*pow(r,2)*pow(Wp,4)*pow(Zer,4)*pow(P,2) - 30720*Qr*pow(Wp,3)*Wpp*pow(Zer,4)*pow(P,2) + 245760*pow(Qr,2)*pow(Wp,4)*Wpp*pow(Zer,4)*pow(P,2) - 6144*pow(Wp,4)*pow(Zer,6)*pow(P,2) + 49152*pow(Qr,2)*pow(r,2)*pow(Wp,4)*Wpp*pow(Zer,6)*pow(P,2) - 3264*pow(Wp,3)*pow(Zer,3)*pow(P,3) + 39168*Qr*pow(Wp,4)*pow(Zer,3)*pow(P,3) + 147456*Qr*pow(Wp,4)*Wpp*pow(Zer,5)*pow(P,3) + 7680*pow(Wp,4)*pow(Zer,4)*pow(P,4) + 64*pow(Wp,2)*pow(Zer,2)*V - 512*Qr*pow(Wp,3)*pow(Zer,2)*V - 896*pow(Wp,3)*pow(Zer,3)*P*V + 2560*Qr*pow(Wp,4)*pow(Zer,3)*P*V + 3072*pow(Wp,4)*pow(Zer,4)*pow(P,2)*V))/(1 - 16*Qr*Wp + 64*pow(Qr,2)*pow(Wp,2) - 20*Wp*Zer*P + 160*Qr*pow(Wp,2)*Zer*P + 96*pow(Wp,2)*pow(Zer,2)*pow(P,2))
	;
	double ep_drq=
		(768*t_Der_P*pow(Wp,3)*pow(Zer,4)*(-1 + 8*Qr*Wp + 8*Wp*Zer*P))/(-1 + 8*Qr*Wp + 12*Wp*Zer*P) - (768*r_Der_P*pow(Wp,3)*pow(Zer,5)*Al*(-r + 4*Qr*r*Wp + 8*r*Wp*Zer*P))/(-1 + 8*Qr*Wp + 12*Wp*Zer*P) - (Al*(1 - 32*Qr*Wp + 384*pow(Qr,2)*pow(Wp,2) - 2048*pow(Qr,3)*pow(Wp,3) + 4096*pow(Qr,4)*pow(Wp,4) - 48*pow(Qr,2)*pow(r,2)*pow(Wp,2)*pow(Zer,2) + 768*pow(Qr,3)*pow(r,2)*pow(Wp,3)*pow(Zer,2) - 3072*pow(Qr,4)*pow(r,2)*pow(Wp,4)*pow(Zer,2) - 96*pow(Wp,2)*pow(Zer,4) + 1536*Qr*pow(Wp,3)*pow(Zer,4) - 6144*pow(Qr,2)*pow(Wp,4)*pow(Zer,4) + 256*pow(Qr,4)*pow(r,4)*pow(Wp,4)*pow(Zer,4) - 32*Wp*Zer*P + 768*Qr*pow(Wp,2)*Zer*P - 6144*pow(Qr,2)*pow(Wp,3)*Zer*P + 16384*pow(Qr,3)*pow(Wp,4)*Zer*P + 896*pow(Qr,2)*pow(r,2)*pow(Wp,3)*pow(Zer,3)*P - 7168*pow(Qr,3)*pow(r,2)*pow(Wp,4)*pow(Zer,3)*P + 1536*pow(Wp,3)*pow(Zer,5)*P - 12288*Qr*pow(Wp,4)*pow(Zer,5)*P + 3072*pow(Qr,2)*pow(r,2)*pow(Wp,3)*Wpp*pow(Zer,5)*P - 24576*pow(Qr,3)*pow(r,2)*pow(Wp,4)*Wpp*pow(Zer,5)*P + 352*pow(Wp,2)*pow(Zer,2)*pow(P,2) - 5632*Qr*pow(Wp,3)*pow(Zer,2)*pow(P,2) + 22528*pow(Qr,2)*pow(Wp,4)*pow(Zer,2)*pow(P,2) - 3840*pow(Qr,2)*pow(r,2)*pow(Wp,4)*pow(Zer,4)*pow(P,2) - 768*pow(Wp,2)*Wpp*pow(Zer,4)*pow(P,2) + 12288*Qr*pow(Wp,3)*Wpp*pow(Zer,4)*pow(P,2) - 49152*pow(Qr,2)*pow(Wp,4)*Wpp*pow(Zer,4)*pow(P,2) - 6144*pow(Wp,4)*pow(Zer,6)*pow(P,2) - 24576*pow(Qr,2)*pow(r,2)*pow(Wp,4)*Wpp*pow(Zer,6)*pow(P,2) - 1536*pow(Wp,3)*pow(Zer,3)*pow(P,3) + 12288*Qr*pow(Wp,4)*pow(Zer,3)*pow(P,3) + 12288*pow(Wp,3)*Wpp*pow(Zer,5)*pow(P,3) - 98304*Qr*pow(Wp,4)*Wpp*pow(Zer,5)*pow(P,3) + 2048*pow(Wp,4)*pow(Zer,4)*pow(P,4) - 49152*pow(Wp,4)*Wpp*pow(Zer,6)*pow(P,4) + 64*pow(Wp,2)*pow(Zer,2)*V - 1024*Qr*pow(Wp,3)*pow(Zer,2)*V + 4096*pow(Qr,2)*pow(Wp,4)*pow(Zer,2)*V + 512*pow(Qr,2)*pow(r,2)*pow(Wp,4)*pow(Zer,4)*V - 1024*pow(Wp,3)*pow(Zer,3)*P*V + 8192*Qr*pow(Wp,4)*pow(Zer,3)*P*V + 4096*pow(Wp,4)*pow(Zer,4)*pow(P,2)*V))/((-1 + 8*Qr*Wp + 8*Wp*Zer*P)*(-1 + 8*Qr*Wp + 12*Wp*Zer*P))
	;
/*---------------------------------------------------------------------------*/
	double eq_dtp=
		0
	;
	double eq_dtq=
		1
	;
	double eq_drp=
	-	Al
	;
	double eq_drq=
	-	Al*Ze
	;
/*---------------------------------------------------------------------------*/
	double A= ep_dtp*eq_dtq- ep_dtq*eq_dtp;
	double B= 
	-	(ep_dtp*eq_drq-ep_dtq*eq_drp)
	-	(ep_drp*eq_dtq-ep_drq*eq_dtp)
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
	const vector<double> &al_v, const vector<double> &ze_v,
	const vector<double>  &f_v, 
	const vector<double>  &p_v, const vector<double>  &q_v,
	vector<double> &ingoing,
	vector<double> &outgoing)
{
	compute_scalar_potentials(f_v);
	int new_exc_i= exc_i;
/*---------------------------------------------------------------------------*/
/* inner two grid points */
	int i=exc_i;
	double x= dx*i;
	double cf= (1-x/cl);
	double r= x/cf;

	double Al= al_v[i];
	double Ze= ze_v[i];
	double P=  p_v[i];
	double Q=  q_v[i];

	double V=   V_v[i];
	double Vp=  Vp_v[i];
	double Wp=  Wp_v[i];
	double Wpp= Wpp_v[i];

	double r_Der_Al= pow(cf,2)*Dx_ptp0_4th(al_v[i+4],al_v[i+3],al_v[i+2],al_v[i+1],al_v[i],dx);
	double r_Der_Ze= pow(cf,2)*Dx_ptp0_4th(ze_v[i+4],ze_v[i+3],ze_v[i+2],ze_v[i+1],ze_v[i],dx);
	double r_Der_P=  pow(cf,2)*Dx_ptp0_4th( p_v[i+4], p_v[i+3], p_v[i+2], p_v[i+1], p_v[i],dx);
	double r_Der_Q=  pow(cf,2)*Dx_ptp0_4th( q_v[i+4], q_v[i+3], q_v[i+2], q_v[i+1], q_v[i],dx);

	int status= compute_radial_characteristic(
		 r,
		 Al,  Ze,
		 P,   Q,
		 V,   Vp,  Wp,  Wpp,
		 r_Der_Al,
		 r_Der_Ze,
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

	Al= al_v[i];
	Ze= ze_v[i];
	P=  p_v[i];
	Q=  q_v[i];

	V=   V_v[i];
	Vp=  Vp_v[i];
	Wp=  Wp_v[i];
	Wpp= Wpp_v[i];

	r_Der_Al= pow(cf,2)*Dx_ptp1_4th(al_v[i+4],al_v[i+3],al_v[i+2],al_v[i+1],al_v[i],dx);
	r_Der_Ze= pow(cf,2)*Dx_ptp1_4th(ze_v[i+4],ze_v[i+3],ze_v[i+2],ze_v[i+1],ze_v[i],dx);
	r_Der_P=  pow(cf,2)*Dx_ptp1_4th( p_v[i+4], p_v[i+3], p_v[i+2], p_v[i+1], p_v[i],dx);
	r_Der_Q=  pow(cf,2)*Dx_ptp1_4th( q_v[i+4], q_v[i+3], q_v[i+2], q_v[i+1], q_v[i],dx);

	status= compute_radial_characteristic(
		 r,
		 Al,  Ze,
		 P,   Q,
		 V,   Vp,  Wp,  Wpp,
		 r_Der_Al,
		 r_Der_Ze,
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

		Al= al_v[i];
		Ze= ze_v[i];
		P=  p_v[i];
		Q=  q_v[i];

		V=   V_v[i];
		Vp=  Vp_v[i];
		Wp=  Wp_v[i];
		Wpp= Wpp_v[i];

		r_Der_Al= pow(cf,2)*Dx_ptc_4th(al_v[i+2],al_v[i+1],al_v[i-1],al_v[i-2],dx);
		r_Der_Ze= pow(cf,2)*Dx_ptc_4th(ze_v[i+2],ze_v[i+1],ze_v[i-1],ze_v[i-2],dx);
		r_Der_P=  pow(cf,2)*Dx_ptc_4th( p_v[i+2], p_v[i+1], p_v[i-1], p_v[i-2],dx);
		r_Der_Q=  pow(cf,2)*Dx_ptc_4th( q_v[i+2], q_v[i+1], q_v[i-1], q_v[i-2],dx);

		status= compute_radial_characteristic(
			 r,
			 Al,  Ze,
			 P,   Q,
			 V,   Vp,  Wp,  Wpp,
			 r_Der_Al,
			 r_Der_Ze,
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
	return;
}
/*===========================================================================*/
void EdGB::compute_eom_rr(
	const int exc_i,
	const vector<double> &al_v, const vector<double> &ze_v,
	const vector<double>  &f_v, 
	const vector<double>  &p_v, const vector<double>  &q_v,
	vector<double> &eom_rr)
{
	compute_scalar_potentials(f_v);
	for (int i=exc_i+2; i<nx-2; ++i) {
		double x= dx*i;
		double cf= (1-x/cl);
		double r= x/cf;

		double Al= al_v[i];
		double Ze= ze_v[i];
		double P=  p_v[i];
		double Q=  q_v[i];

		double Zer= Ze/r;
		double Qr=  Q/r;

		double V=   V_v[i];
		double Vp=  Vp_v[i];
		double Wp=  Wp_v[i];
		double Wpp= Wpp_v[i];

		double r_Der_Al= pow(cf,2)*Dx_ptc_4th(al_v[i+2],al_v[i+1],al_v[i-1],al_v[i-2],dx);
		double r_Der_Ze= pow(cf,2)*Dx_ptc_4th(ze_v[i+2],ze_v[i+1],ze_v[i-1],ze_v[i-2],dx);
		double r_Der_P=  pow(cf,2)*Dx_ptc_4th( p_v[i+2], p_v[i+1], p_v[i-1], p_v[i-2],dx);
		double r_Der_Q=  pow(cf,2)*Dx_ptc_4th( q_v[i+2], q_v[i+1], q_v[i-1], q_v[i-2],dx);

		double t_Der_P= compute_p_k(
			r,
			Al, Ze,
			P,  Q,
			V,  Vp, Wp, Wpp,
			r_Der_Al,
			r_Der_Ze,
			r_Der_P,
			r_Der_Q
		);
		double t_Der_Ze= compute_ze_free_k(
			r,
			Al, Ze,
			P,  Q,
			V,  Vp, Wp, Wpp,
			r_Der_Al, r_Der_Ze,
			r_Der_P,  r_Der_Q
		);
		eom_rr[i]=
			-pow(Zer,2) + 8*r*r_Der_P*Wp*pow(Zer,3) - (8*t_Der_P*Wp*pow(Zer,2))/Al - 8*Wpp*pow(Zer,2)*pow(P,2) + r_Der_Ze*(-2*Zer + 16*Qr*Wp*Zer + 16*Wp*pow(Zer,2)*P) + t_Der_Ze*(2/(r*Al) - (16*Qr*Wp)/(r*Al) - (16*Wp*Zer*P)/(r*Al)) + r_Der_Al*(2/(r*Al) - (16*Qr*Wp)/(r*Al) + (8*Qr*r*Wp*pow(Zer,2))/Al - (16*Wp*Zer*P)/(r*Al)) + (-(pow(Qr,2)*pow(r,2)) - pow(P,2) + 2*V)/2.
		;
	}
	return;
}
