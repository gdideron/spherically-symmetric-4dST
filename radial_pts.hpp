#ifndef _RADIAL_PTS_HPP_
#define _RADIAL_PTS_HPP_

#include <vector>

/*===========================================================================*/
class Radial_pts 
{
public:
	std::vector<double> x;
	std::vector<double> r;

	const double cl;

	Radial_pts(
		const int nx, const double dx, const double cl
	);
	~Radial_pts(void);
private:
}; 

#endif
