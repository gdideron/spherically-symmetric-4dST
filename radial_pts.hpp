#ifndef _RADIAL_PTS_HPP_
#define _RADIAL_PTS_HPP_

#include <vector>
using std::vector;

/*===========================================================================*/
class Radial_pts 
{
public:
	vector<double> x;
	vector<double> r;

	const double cl;

	Radial_pts(
		const int nx, const double dx, const double cl
	);
	~Radial_pts(void);
private:
}; 

#endif
