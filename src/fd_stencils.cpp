#include <vector>
using std::vector;
#include <string>
using std::string;

#include "fd_stencils.hpp"

/*===========================================================================*/
inline double Dx_ptm1_4th(
   const double vp1, const double v0, const double vm1,
   const double vm2, const double vm3,
   const double dx);
/*===========================================================================*/
inline double Dx_ptc_4th(
   const double vp2, const double vp1, const double vm1, const double vm2,
   const double dx);
/*===========================================================================*/
inline double Dx_ptp1_4th(
   const double vp3, const double vp2, const double vp1,
   const double v0, const double vm1,
   const double dx);
/*===========================================================================*/
inline double Dx_ptp0_4th(
   const double vp4, const double vp3, const double vp2,
   const double vp1, const double v0,
   const double dx);
/*===========================================================================*/
inline double make_Dx_zero(
   const double vp4, const double vp3, const double vp2,
   const double vp1);
/*===========================================================================*/
void KO_filter(
   const int exc_i, const int nx,
   const string type, vector<double> &vec)
{
   double eps= 0.25;

   for (int i=exc_i+3; i<nx-3; ++i) {
      vec[i]+= (eps/64)*(
	 vec[i+3]
      -	 6*vec[i+2]
      +	 15*vec[i+1]
      -	 20*vec[i]
      +	 15*vec[i-1]
      -	 6*vec[i-2]
      +	 vec[i-3]
      );
   }
/*---------------------------------------------------------------------------*/
   if (exc_i>0) return;
/*---------------------------------------------------------------------------*/
   if (type=="even") {
      vec[2]+= (eps/64)*(
	 vec[5]
      -	 6*vec[4]
      +	 15*vec[3]
      -	 20*vec[2]
      +	 15*vec[1]
      -	 6*vec[0]
      +	 vec[1]
      );
      vec[1]+= (eps/64)*(
         vec[4]
      -  6*vec[3]
      +  15*vec[2]
      -  20*vec[1]
      +  15*vec[0]
      -  6*vec[1]
      +  vec[2]
      );
      vec[0]+= (eps/64)*(
         vec[3]
      -  6*vec[2]
      +  15*vec[1]
      -  20*vec[0]
      +  15*vec[1]
      -  6*vec[2]
      +  vec[3]
      );
   } else
   if (type=="odd") {
      vec[2]+= (eps/64)*(
         vec[5]
      -	 6*vec[4]
      +	 15*vec[3]
      -	 20*vec[2]
      +	 15*vec[1]
      -	 6*vec[0]
      +	 (-vec[1])
      );
      vec[1]+= (eps/64)*(
         vec[4]
      -  6*vec[3]
      +  15*vec[2]
      -  20*vec[1]
      +  15*vec[0]
      -  (-6*vec[1])
      +  (-vec[2])
      );
      vec[0]+= (eps/64)*(
         vec[3]
      -	 6*vec[2]
      +	 15*vec[1]
      -	 20*vec[0]
      +	 (-15*vec[1])
      -	 (-6*vec[2])
      +	 (-vec[3])
      );
   } else {
      /* do nothing */
   }
}   
/*===========================================================================*/
