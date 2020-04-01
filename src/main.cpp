#include <cassert>
#include <ctime>
#include <cmath>
#include <string>
using std::string;
#include <iomanip>
using std::setprecision;
using std::setw;
#include <iostream>
using std::cout; 
using std::endl;
#include <string>
using std::string;
#include <vector>
using std::vector;

#include "sim_params.hpp"
#include "field.hpp"
#include "radial_pts.hpp"
#include "edgb.hpp"
#include "initial_data.hpp"
#include "io_sdf.hpp"
#include "io_csv.hpp"
#include "indep_res.hpp"
/*===========================================================================*/
int main(int argc, char **argv)
{
	assert(argc==2);
	const string output_dir= argv[1];

	clock_t start= std::clock();
/*--------------------------------------------------------------------------*/	
	const Sim_params sp(output_dir) ;
	const Radial_pts rp(sp.nx, sp.dx, sp.cl);
/*--------------------------------------------------------------------------*/	
	EdGB edgb(
		sp.dt, sp.dx, sp.cl, sp.nx, 
		sp.mu, sp.la, 
		sp.gbc1, sp.gbc2,
		rp.r
	);
/*--------------------------------------------------------------------------*/	
/* evolution fields */
/*--------------------------------------------------------------------------*/	
	Field phi_f("phi_f",sp.nx,0);
	Field phi_q("phi_q",sp.nx,0);
	Field phi_p("phi_p",sp.nx,0);
	
	Field al("al",sp.nx,1);
	Field ze("ze",sp.nx,0);

	Field res_q( "res_q",  sp.nx,0);
	Field eom_rr("res_rr", sp.nx,0);

	Field  ingoing_c( "ingoing_c",sp.nx,0);
	Field outgoing_c("outgoing_c",sp.nx,0);

	Field ncc("ncc",sp.nx,0);
/*--------------------------------------------------------------------------*/		
/* initial data */
/*--------------------------------------------------------------------------*/		
	set_initial_data(sp, rp.r, al, ze, phi_f, phi_p, phi_q);

	edgb.solve_metric_fields(sp.initial_exc_i,
		phi_f, phi_p, phi_q, al, ze
	);
	int exc_i= sp.initial_exc_i;

	double initial_asymptotic_mass= rp.r[sp.nx-8]*pow(ze.np1[sp.nx-8],2)/2; 
/*--------------------------------------------------------------------------*/		
/* write to file */
/*--------------------------------------------------------------------------*/		
	Sdf sdf(output_dir, sp.nx, rp.x);
	
	double grid_time= 0;

	sdf.write(grid_time, phi_f);
	sdf.write(grid_time, phi_p);
	sdf.write(grid_time, phi_q);

	sdf.write(grid_time, al);
	sdf.write(grid_time, ze);

/*--------------------------------------------------------------------------*/		
	compute_indep_res_q(
		exc_i,
		sp.nx, sp.dx, sp.cl,
		rp.r,
		phi_f.np1, phi_q.np1,
		res_q.np1
	);
	edgb.compute_eom_rr(
		exc_i,
		al.np1, ze.np1,
		phi_f.np1, phi_p.np1, phi_q.np1,
		eom_rr.np1
	);
	sdf.write(grid_time, res_q);
	sdf.write(grid_time, eom_rr);
	
	edgb.compute_ncc(
		exc_i,
		sp.nx, sp.dx, sp.cl, 
		rp.r, 
		al.np1,       ze.np1,
		phi_p.np1, phi_q.np1,
		ncc.np1
	);
	sdf.write(grid_time, ncc);

	edgb.compute_radial_characteristics(
		exc_i,
		al.np1, ze.np1,
		phi_f.np1, phi_p.np1, phi_q.np1,
		ingoing_c.np1, outgoing_c.np1
	);
	sdf.write(grid_time,  ingoing_c);
	sdf.write(grid_time, outgoing_c);
/*--------------------------------------------------------------------------*/		
	cout<<setw(10)<<0<<"\t";
	cout<<setw(10)<<rp.r[sp.nx-8]*pow(ze.np1[sp.nx-8],2)/2<<"\t";
	cout<<setw(10)<<rp.r[sp.nx-8]*phi_f.np1[sp.nx-8]<<"\t";
	cout<<endl;
/*--------------------------------------------------------------------------*/		
/* 	evolve in time and write to file */
/*--------------------------------------------------------------------------*/		
	for (int tC=0; tC<sp.nt; ++tC) {
		grid_time+= sp.dt/initial_asymptotic_mass;

		edgb.compute_radial_characteristics(
			exc_i,
			al.np1, ze.np1,
			phi_f.np1, phi_p.np1, phi_q.np1,
			ingoing_c.np1, outgoing_c.np1
		);
		edgb.time_step(exc_i, al, ze, phi_f, phi_p, phi_q);

		if (tC%(sp.t_step_save)==0) {
/*--------------------------------------------------------------------------*/		
			al.check_isfinite(grid_time);
			ze.check_isfinite(grid_time);
			phi_f.check_isfinite(grid_time);
			phi_p.check_isfinite(grid_time);
			phi_q.check_isfinite(grid_time);
/*--------------------------------------------------------------------------*/		
			compute_indep_res_q(
				exc_i,
				sp.nx, sp.dx, sp.cl,
				rp.r,
				phi_f.np1, phi_q.np1,
				res_q.np1
			);
			edgb.compute_eom_rr(
				exc_i,
				al.np1, ze.np1,
				phi_f.np1, phi_p.np1, phi_q.np1,
				eom_rr.np1
			);
	
			edgb.compute_ncc(
				exc_i,
				sp.nx, sp.dx, sp.cl, 
				rp.r, 
				al.np1,       ze.np1,
				phi_p.np1, phi_q.np1,
				ncc.np1
			);
/*--------------------------------------------------------------------------*/		
			phi_f.set_to_val(0,exc_i-1,0);
			phi_p.set_to_val(0,exc_i-1,0);
			phi_q.set_to_val(0,exc_i-1,0);

			al.set_to_val(0,exc_i-1,0);
			ze.set_to_val(0,exc_i-1,0);

			res_q.set_to_val( 0,exc_i-1,0);
			eom_rr.set_to_val(0,exc_i-1,0);

			ncc.set_to_val( 0,exc_i-1,0);

			ingoing_c.set_to_val( 0,exc_i-1,0);
			outgoing_c.set_to_val(0,exc_i-1,0);
/*--------------------------------------------------------------------------*/		
			sdf.write(grid_time, phi_f);
			sdf.write(grid_time, phi_p);
			sdf.write(grid_time, phi_q);

			sdf.write(grid_time, al);
			sdf.write(grid_time, ze);

			sdf.write(grid_time, res_q);
			sdf.write(grid_time, eom_rr);
			
			sdf.write(grid_time, ncc);
			sdf.write(grid_time,  ingoing_c);
			sdf.write(grid_time, outgoing_c);
/*--------------------------------------------------------------------------*/		
			cout<<setw(10)<<tC*sp.dt/initial_asymptotic_mass<<"\t";
			cout<<setw(10)<<rp.r[sp.nx-8]*pow(ze.np1[sp.nx-8],2)/2<<"\t";
			cout<<setw(10)<<rp.r[sp.nx-8]*phi_f.np1[sp.nx-8]<<"\t";
			cout<<endl;
/*--------------------------------------------------------------------------*/		
		}
		phi_f.shift_time_step();
		phi_p.shift_time_step();
		phi_q.shift_time_step();

		al.shift_time_step();
		ze.shift_time_step();
	}
	clock_t end= std::clock();
	cout<<endl;
	cout<<"run_finished_successfully"<<endl;
	cout<<"time (to run, in sec): "<<difftime(end,start)/CLOCKS_PER_SEC<<endl;
	cout<<"initial asymptotic mass: "<<initial_asymptotic_mass<<endl;
	cout<<"time T: "<<sp.nt*sp.dt/initial_asymptotic_mass<<endl;
	return EXIT_SUCCESS;
}
