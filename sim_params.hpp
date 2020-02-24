#ifndef _SIM_PARAMS_HPP_
#define _SIM_PARAMS_HPP_

#include <string>
using std::string;

/*===========================================================================*/
class Sim_params
{
public:	
	int initial_exc_i;

	int nt;
	int nx;
	int t_step_save;

	double dt;
	double dx;

	double cl;
	
	int phi_pt;
/*---------------------------------------------------------------------------*/
/* for the potentials */
	double mu;
	double la;

	double gbc1;
	double gbc2;
/*---------------------------------------------------------------------------*/
/* for the initial data */
	string initial_data_type;

	double bh_mass;

	double charge;

	double amp;
	double r_l;
	double r_u;

	Sim_params(const string output_dir);
	~Sim_params(void);	
private:
	string read_sim_params(const string output_dir, const string find_this_var);
};
#endif
