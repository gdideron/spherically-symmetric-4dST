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
      sp.V_1,  sp.V_2,  sp.V_3,  sp.V_4, 
      sp.Al_0, 
      sp.Al_1, sp.Al_2, sp.Al_3, sp.Al_4, 
      sp.Be_1, sp.Be_2, sp.Be_3, sp.Be_4, 
      rp.r
   );
/*--------------------------------------------------------------------------*/	
/* evolution fields */
/*--------------------------------------------------------------------------*/	
   Field phi_f("phi_f",sp.nx,0);
   Field phi_q("phi_q",sp.nx,0);
   Field phi_p("phi_p",sp.nx,0);

   Field N("N",sp.nx,1);
   Field S("S",sp.nx,0);

   Field res_q( "res_q",  sp.nx,0);
   Field eom_rr("res_rr", sp.nx,0);

   Field  ingoing_c( "ingoing_c",sp.nx,0);
   Field outgoing_c("outgoing_c",sp.nx,0);

   Field ncc("ncc",sp.nx,0);
/*--------------------------------------------------------------------------*/		
/* initial data */
/*--------------------------------------------------------------------------*/		
   set_initial_data(sp, rp.r, N, S, phi_f, phi_p, phi_q);

   edgb.solve_metric_fields(
      sp.initial_exc_i,
      phi_f, phi_p, phi_q, N, S
   );
   if (sp.initial_data_type=="scalarized_bh") {
      time_symmetric(
         edgb,sp,
         phi_f, phi_q,
         N, S,
         phi_p);
   }
   int exc_i= sp.initial_exc_i;
   double initial_asymptotic_mass= rp.r[sp.nx-8]*pow(S.np1[sp.nx-8],2)/2; 
/*--------------------------------------------------------------------------*/		
/* write to file */
/*--------------------------------------------------------------------------*/		
   Csv csv(output_dir);

   double grid_time= 0;

   csv.write(grid_time, phi_f);
   csv.write(grid_time, phi_p);
   csv.write(grid_time, phi_q);

   csv.write(grid_time, N);
   csv.write(grid_time, S);
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
      N.np1,     S.np1,
      phi_f.np1, phi_p.np1, phi_q.np1,
      eom_rr.np1
   );
   csv.write(grid_time, res_q);
   csv.write(grid_time, eom_rr);

   edgb.compute_ncc(
      exc_i,
      sp.nx, sp.dx, sp.cl, 
      rp.r, 
      N.np1,     S.np1,
      phi_p.np1, phi_q.np1,
      ncc.np1
   );
   csv.write(grid_time, ncc);

   edgb.compute_radial_characteristics(
      exc_i,
      N.np1, S.np1,
      phi_f.np1, phi_p.np1, phi_q.np1,
      ingoing_c.np1, outgoing_c.np1
   );
   csv.write(grid_time, ingoing_c);
   csv.write(grid_time, outgoing_c);
/*--------------------------------------------------------------------------*/		
   cout<<setw(10)<<0<<"\t";
   cout<<setw(10)<<rp.r[sp.nx-8]*pow(S.np1[sp.nx-8],2)/2<<"\t";
   cout<<setw(10)<<rp.r[sp.nx-8]*phi_f.np1[sp.nx-8]<<"\t";
   cout<<endl;
/*--------------------------------------------------------------------------*/		
/* evolve in time and write to file:
*  stop once scalar field has settled down and have evolved
*  for minimum number of time steps  */
/*--------------------------------------------------------------------------*/		
   int tC= 0;
   while (
   (tC<sp.nt) 
   || 	(((phi_f.np1[exc_i+10]-phi_f.n[exc_i+10])/sp.dt)>1e-16)
   ) {
      tC+= 1;
      grid_time+= sp.dt/initial_asymptotic_mass;

      edgb.compute_radial_characteristics(
      exc_i,
      N.np1, S.np1,
      phi_f.np1, phi_p.np1, phi_q.np1,
      ingoing_c.np1, outgoing_c.np1
      );
      edgb.time_step(exc_i, N, S, phi_f, phi_p, phi_q);

      if (tC%(sp.t_step_save)==0) {
/*--------------------------------------------------------------------------*/		
         N.check_isfinite(grid_time);
         S.check_isfinite(grid_time);
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
            N.np1, S.np1,
            phi_f.np1, phi_p.np1, phi_q.np1,
            eom_rr.np1
         );

         edgb.compute_ncc(
            exc_i,
            sp.nx, sp.dx, sp.cl, 
            rp.r, 
            N.np1,     S.np1,
            phi_p.np1, phi_q.np1,
            ncc.np1
         );
/*--------------------------------------------------------------------------*/		
         phi_f.set_to_val(0,exc_i-1,0);
         phi_p.set_to_val(0,exc_i-1,0);
         phi_q.set_to_val(0,exc_i-1,0);

         N.set_to_val(0,exc_i-1,0);
         S.set_to_val(0,exc_i-1,0);

         res_q.set_to_val( 0,exc_i-1,0);
         eom_rr.set_to_val(0,exc_i-1,0);

         ncc.set_to_val( 0,exc_i-1,0);

         ingoing_c.set_to_val( 0,exc_i-1,0);
         outgoing_c.set_to_val(0,exc_i-1,0);
/*--------------------------------------------------------------------------*/		
         csv.write(grid_time, phi_f);
         csv.write(grid_time, phi_p);
         csv.write(grid_time, phi_q);

         csv.write(grid_time, N);
         csv.write(grid_time, S);

         csv.write(grid_time, res_q);
         csv.write(grid_time, eom_rr);

         csv.write(grid_time, ncc);
         csv.write(grid_time,  ingoing_c);
         csv.write(grid_time, outgoing_c);
/*--------------------------------------------------------------------------*/		
         cout<<setw(10)<<tC*sp.dt/initial_asymptotic_mass<<"\t";
         cout<<setw(10)<<rp.r[sp.nx-8]*pow(S.np1[sp.nx-8],2)/2<<"\t";
         cout<<setw(10)<<rp.r[sp.nx-8]*phi_f.np1[sp.nx-8]<<"\t";
         cout<<setw(10)<<((phi_f.np1[exc_i+10]-phi_f.n[exc_i+10])/sp.dt)<<"\t";
         cout<<endl;
/*--------------------------------------------------------------------------*/		
      }
      phi_f.shift_time_step();
      phi_p.shift_time_step();
      phi_q.shift_time_step();

      N.shift_time_step();
      S.shift_time_step();
   }
   clock_t end= std::clock();
   cout<<endl;
   cout<<"run_finished_successfully"<<endl;
   cout<<"time (to run, in min): "<<difftime(end,start)/(CLOCKS_PER_SEC*60)<<endl;
   cout<<"initial asymptotic mass: "<<initial_asymptotic_mass<<endl;
   cout<<"time T: "<<sp.nt*sp.dt/initial_asymptotic_mass<<endl;
   return EXIT_SUCCESS;
}
