#include <vector>
using std::vector;

#include "radial_pts.hpp"

/*===========================================================================*/
Radial_pts::Radial_pts(
	const int nx, const double dx, const double cl)
: x(nx),
  r(nx),

  cl{cl}
{
   for (int i=0; i<nx-1; ++i) {
      x[i]=  dx*i;
      r[i]=  x[i]/(1-x[i]/cl);
   }
   x[nx-1]= cl;
   r[nx-1]= 1e100;
}
/*===========================================================================*/
Radial_pts::~Radial_pts(void)
{
}
