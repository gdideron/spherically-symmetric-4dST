#ifndef __INITIAL_DATA_HPP__
#define __INITIAL_DATA_HPP__

#include <vector>

#include "sim_params.hpp"
#include "field.hpp"
#include "edgb.hpp"

/*===========================================================================*/
void set_initial_data(
	const Sim_params &sp,
	const std::vector<double> &rvec,
	Field &al,Field &ze,
	Field &f, Field &p, Field &q
);
/*===========================================================================*/
void time_symmetric(EdGB &edgb,
	const Sim_params &sp,
	const Field &f, const Field &q,
	Field &al,Field &ze,
	Field &p)
;

#endif /* __INITIAL_DATA_HPP__ */
