#include <iostream>
using std::cout; 
using std::endl;
#include <string>
using std::string;
#include <vector>
using std::vector; 
#include <cmath>
using std::isfinite;

#include "field.hpp"

/*===========================================================================*/
Field::Field(const string file_name, const int nx, const double init_val)
: name{file_name},

  n(nx,  init_val),
  np1(nx,init_val), 
  
  inter_2(nx,init_val), 
  inter_3(nx,init_val), 
  inter_4(nx,init_val),
  
  nx{nx}
{
}
/*===========================================================================*/
Field::~Field(void)
{
}
/*===========================================================================*/
void Field::shift_time_step(void)
{
	for (int i=0; i<nx; ++i) {
		n[i]= np1[i];
	}
} 
/*===========================================================================*/
void Field::set_to_val(const int min, const int max, const double val)
{
	for (int i=min; i<=max; ++i) {
		n[i]= val;
		np1[i]= val;

		inter_2[i]= val;
		inter_3[i]= val;
		inter_4[i]= val;
	}
}
/*===========================================================================*/
void Field::check_field_isfinite(
	const double time, const string level,
	const vector<double> &vec)
{
	size_t index=0;
	for (auto x: vec) {
		if (!isfinite(x)) {
			cout<<"NaN("<<endl;
			cout<<name<<" ";
			cout<<level<<" ";
			cout<<"time "<<time;
			cout<<"index "<<index;
			cout<<endl;
			std::quick_exit(0);
		}
		index++;
	}	
}
/*===========================================================================*/
void Field::check_isfinite(const double time)
{
	check_field_isfinite(time,"n",n);

	check_field_isfinite(time,"inter_2",inter_2);
	check_field_isfinite(time,"inter_3",inter_3);
	check_field_isfinite(time,"inter_4",inter_4);
	
	check_field_isfinite(time,"np1",np1);
} 
