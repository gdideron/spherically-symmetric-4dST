#include <cassert>
#include <fstream>
using std::ifstream;
#include <iostream>
using std::cout; 
using std::endl;
#include <string>
using std::string;
#include <iomanip>
#include <fstream>

#include "sim_params.hpp"
/*===========================================================================*/
string Sim_params::read_sim_params(const string output_dir, const string find_this_var)
{
	ifstream infile(output_dir+"/sim_params.txt");
	assert(infile.good());
	
	string name, val;

	while (infile >> name >> val) {
		if (name==find_this_var) {
			return val;
		}
	}
	cout << "ERROR: did not find " << find_this_var << endl;
	std::quick_exit(0);
}
/*===========================================================================*/
Sim_params::Sim_params(const string output_dir)
{
	initial_exc_i= stoi(read_sim_params(output_dir,"initial_exc_i"));

	nt= stoi(read_sim_params(output_dir,"nt"));
	nx= stoi(read_sim_params(output_dir,"nx"));

	t_step_save= stoi(read_sim_params(output_dir,"t_step_save"));

	dt= stod(read_sim_params(output_dir,"dt"));
	dx= stod(read_sim_params(output_dir,"dx"));

	cl= stod(read_sim_params(output_dir,"compactification_length"));
	
	phi_pt= stoi(read_sim_params(output_dir,"phi_pt"));
/*---------------------------------------------------------------------------*/
/* for the potentials */
	mu= stod(read_sim_params(output_dir,"mu"));
	la= stod(read_sim_params(output_dir,"la"));

	gbc1= stod(read_sim_params(output_dir,"gbc1"));
	gbc2= stod(read_sim_params(output_dir,"gbc2"));
/*---------------------------------------------------------------------------*/
/* for the initial data */
	initial_data_type= read_sim_params(output_dir,"initial_data_type");

	bh_mass= stod(read_sim_params(output_dir,"bh_mass"));
	charge= stod(read_sim_params(output_dir,"charge"));

	amp= stod(read_sim_params(output_dir,"amp"));
	r_c= stod(read_sim_params(output_dir,"r_c"));
	r_l= stod(read_sim_params(output_dir,"r_l"));
	r_u= stod(read_sim_params(output_dir,"r_u"));
	r_w= stod(read_sim_params(output_dir,"r_w"));
}
/*===========================================================================*/
Sim_params::~Sim_params(void)
{
}
