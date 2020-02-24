#ifndef _FD_STENCILS_H_
#define _FD_STENCILS_H_

#include <vector>
using std::vector;
#include <string>
using std::string;

/*===========================================================================*/
void KO_filter(
	const int exc_i, const int nx,
	const string type, vector<double> &vec);
/*===========================================================================*/
inline double Dx_ptm1_4th(
	const double vp1, const double v0, const double vm1,
	const double vm2, const double vm3,
	const double dx)
{
	return ( 
		(1./4.)*  vp1
	+	(5./6.)*  v0
	+	(-3./2.)* vm1
	+	(1./2.)*  vm2
	+	(-1./12.)*vm3
	)/dx
	;
}
/*===========================================================================*/
inline double Dx_ptc_4th(
	const double vp2, const double vp1, const double vm1, const double vm2,
	const double dx)
{
	return (
		(-1./12.)*vp2
	+	(2./3.)*  vp1
	+	(-2./3.)* vm1
	+	(1./12.)* vm2
	)/dx
	;
}
/*===========================================================================*/
inline double Dx_ptp1_4th(
	const double vp3, const double vp2, const double vp1,
	const double v0, const double vm1,
	const double dx)
{
	return ( 
		(1./12.)*vp3
	+	(-1./2.)*vp2
	+	(3./2.)* vp1
	+	(-5./6.)*v0
	+	(-1./4.)*vm1
	)/dx
	;
}
/*===========================================================================*/
inline double Dx_ptp0_4th(
	const double vp4, const double vp3, const double vp2,
	const double vp1, const double v0,
	const double dx)
{
	return ( 
		(-1./4.)* vp4
	+	(4./3.)*  vp3
	+	(-3.)*    vp2
	+	(4.)*     vp1
	+	(-25./12.)*v0
	)/dx
	;
}
/*===========================================================================*/
inline double make_Dx_zero(
	const double vp4, const double vp3, const double vp2,
	const double vp1)
{
	return ( 
		(-1./4.)* vp4
	+	(4./3.)*  vp3
	+	(-3.)*    vp2
	+	(4.)*     vp1
	)*(12./25.)
	;
}

#endif
