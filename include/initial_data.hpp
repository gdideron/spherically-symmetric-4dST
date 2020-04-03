#include <vector>

#include "sim_params.hpp"
#include "field.hpp"

/*===========================================================================*/
void set_initial_data(
	const Sim_params &sp,
	const std::vector<double> &rvec,
	Field &al,Field &ze,
	Field &f, Field &p, Field &q
);
/*===========================================================================*/
void time_symmetric_p(const Sim_params &sp,const Field &ze,const Field &q, Field &p);
