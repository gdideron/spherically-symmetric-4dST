#ifndef _FIELD_HPP_
#define _FIELD_HPP_

#include <string>
using std::string;
#include <vector>
using std::vector; 

/*===========================================================================*/
class Field
{ 
public:
	const string name;

	const string type;

	vector<double> n;
	vector<double> np1;
  
	vector<double> inter_2; 
	vector<double> inter_3; 
	vector<double> inter_4;

	Field(
		const string file_name,
		const int nx, const double init_val
	);
	~Field(void);

	void set_to_val(const int min, const int max, const double val);
	void shift_time_step(void);
	void check_isfinite(void);
private:
	const int nx;

	void check_field_isfinite(const string message, vector<double> &vec);
};

#endif
