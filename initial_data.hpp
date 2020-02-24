#include <vector>
using std::vector;

#include "sim_params.hpp"
#include "field.hpp"

/*===========================================================================*/
void set_initial_data(
	const Sim_params &sp,
	const vector<double> &rvec,
	Field &al,Field &ze,
	Field &f, Field &p, Field &q
);
