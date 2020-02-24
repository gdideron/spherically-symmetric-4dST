#ifndef _IO_SDF_HPP_
#define _IO_SDF_HPP_

#include <string>
using std::string;
#include <vector>
using std::vector; 

#include "field.hpp"

/*===========================================================================*/
/* see Man pages for bbhutil utilities on
   http://laplace.physics.ubc.ca/Group/Software.html 
   for more information on sdf files */
/*===========================================================================*/
class Sdf 
{
public:
	void write(const double grid_time, const class Field &field);
	
	Sdf(const string output, const int nx, const vector<double> &x_pts);
	~Sdf(void);
private:
	string output_dir;

	int nx; 

	int rank;
	int shape[1];

	char *coord_names;

	double *coords;
	double *vals;
};

#endif
