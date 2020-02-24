#include <iostream>
using std::cout; 
using std::endl;
#include <cassert>
#include <string>
using std::string;
#include <vector>
using std::vector; 
/* for sdf output */
#include <bbhutil.h>

#include "io_sdf.hpp"
#include "field.hpp"
	
/*===========================================================================*/
/* see Man pages for bbhutil utilities on
   http://laplace.physics.ubc.ca/Group/Software.html 
   for more information on sdf files */
/*===========================================================================*/
Sdf::Sdf(const string output, const int nx, const vector<double> &x_pts)
: output_dir{output},
  nx{nx},
  rank{1},
  shape{nx},
  coords((double*)calloc(nx,sizeof(double))),
  vals(  (double*)calloc(nx,sizeof(double)))
{	
	coord_names= (char *)"x";
	
	for (int i=0; i<nx; ++i) {
		coords[i]= x_pts[i];
	}
}
/*===========================================================================*/
Sdf::~Sdf(void)
{
	free(coords);
	free(vals);
}
/*===========================================================================*/
/* check==-1: did not find correct parameters
   check== 0: syntax or memory error
   check==+1: saved successfully */
/*===========================================================================*/
void Sdf::write(double grid_time, const Field &f)
{
	string file_name= output_dir+"/"+f.name;
	char *c_file_name= &file_name[0];

	int check= 0;

	assert(nx==(int)f.np1.size());
/*---------------------------------------------------------------------------*/
	for (int i=0; i<nx; ++i) {
		vals[i]= f.np1[i];
	}		
	check= gft_out_full(c_file_name, grid_time, shape, coord_names, rank, coords, vals);
	assert(check!=-1); 
	assert(check!= 0); 
	assert(check== 1); 
/*---------------------------------------------------------------------------*/
	return;
}
