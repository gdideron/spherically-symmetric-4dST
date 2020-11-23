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
/*---------------------------------------------------------------------------*/
/* for the potentials */
   V_1= stod(read_sim_params(output_dir,"V_1"));
   V_2= stod(read_sim_params(output_dir,"V_2"));
   V_3= stod(read_sim_params(output_dir,"V_3"));
   V_4= stod(read_sim_params(output_dir,"V_4"));

   Al_0= stod(read_sim_params(output_dir,"Al_0"));
   Al_1= stod(read_sim_params(output_dir,"Al_1"));
   Al_2= stod(read_sim_params(output_dir,"Al_2"));
   Al_3= stod(read_sim_params(output_dir,"Al_3"));
   Al_4= stod(read_sim_params(output_dir,"Al_4"));

   Be_1= stod(read_sim_params(output_dir,"Be_1"));
   Be_2= stod(read_sim_params(output_dir,"Be_2"));
   Be_3= stod(read_sim_params(output_dir,"Be_3"));
   Be_4= stod(read_sim_params(output_dir,"Be_4"));
/*---------------------------------------------------------------------------*/
/* for the initial data */
   initial_data_type= read_sim_params(output_dir,"initial_data_type");

   bh_mass= stod(read_sim_params(output_dir,"bh_mass"));

   amp= stod(read_sim_params(output_dir,"amp"));
   r_l= stod(read_sim_params(output_dir,"r_l"));
   r_u= stod(read_sim_params(output_dir,"r_u"));
}
/*===========================================================================*/
Sim_params::~Sim_params(void)
{
}
