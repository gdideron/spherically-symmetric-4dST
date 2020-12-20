#!/usr/bin/env python
#=============================================================================
##
## Parameter file
## Usage: 
## ./setup [run_type] (e.g. ./setup basic_run)
##
import sys, time
from sim_class import Sim
#=============================================================================
args= sys.argv
#=============================================================================
## input paramters: set these by hand 
#=============================================================================
sim= Sim(args)
#-----------------------------------------------------------------------------
sim.binary= 'default.run'
#-----------------------------------------------------------------------------
sim.compactification_length= float(100) 
#-----------------------------------------------------------------------------
sim.evolve_time=   float(3200)  ### in units of initial black hole mass for ze field 
sim.num_saved_times= int(400)
sim.cfl= 0.25
#-----------------------------------------------------------------------------
## couplings
#-----------------------------------------------------------------------------
sim.V_1=  0.0
sim.V_2=  0.0
sim.V_3=  0.0
sim.V_4=  0.0

sim.Al_0=  0.0
sim.Al_1=  0.0
sim.Al_2=  0.0
sim.Al_3=  0.0
sim.Al_4=  0.0

sim.Be_1=  0.0
sim.Be_2=  0.0
sim.Be_3=  0.0
sim.Be_4=  0.0

sim.Be_exp2=  18.75
#sim.Be_exp2=  20.00
#-----------------------------------------------------------------------------
sim.phi_r= 15  ### where measuring phi (radial distance)
#-----------------------------------------------------------------------------
## excision point as ratio of apparent horizon
sim.excision_ratio = 0.8
#-----------------------------------------------------------------------------
## initial data
#-----------------------------------------------------------------------------
sim.bh_mass= float(5.0)
#-----------------------------------------------------------------------------
## for the noncompact scalar profile
#sim.initial_data_type= str("scalarized_bh")
#-----------------------------------------------------------------------------
## for the Gaussian-like pulse
#sim.initial_data_type= str("bump")
sim.initial_data_type= str("bump_with_bh")

sim.amp= float(1.0e-6)
sim.r_l= float(24.0)
sim.r_u= float(32.0)
#-----------------------------------------------------------------------------
sim.nx= pow(2,14)+1 
#-----------------------------------------------------------------------------
sim.set_derived_params()
#=============================================================================
## for slurm script
#=============================================================================
sim.walltime= '16:00:00' ### (hh:mm:ss)

sim.data_dir= '/home/jripley/spherically-symmetric-4dST/output'
#sim.data_dir= '/tigress/jripley/edgb'

sim.slurm= False
#=============================================================================
if (sim.run_type == 'basic_run'):
   sim.launch()
#=============================================================================
elif (sim.run_type == 'scan'):
   mu_hat=0.01
   try:
      os.makedirs(self.data_dir)
   except FileExistsError:
      pass
   for la_hat in [0,0.2,0.4,0.6,0.8]:
      for gbc2 in [300,305,310,315,317.5,320,322.5,325]:
         sim.mu_hat= mu_hat
         sim.la_hat= la_hat
         sim.gbc2=  gbc2 
         sim.launch()
         time.sleep(2*60*60)
#=============================================================================
elif (sim.run_type == 'res_study'):
   num_res= int(input('number of resolutions '))
   for i in range(num_res):	
      sim.launch()
      sim.nx= 2*(sim.nx-1)+1 
#=============================================================================
### varying eta with a fixed initial phi amplitude
#=============================================================================
elif (sim.run_type == 'search_for_elliptic'):

   sim.mu_hat= 0.1
   sim.la_hat= 3.2

   sim.data_dir= '/mnt/grtheory/jripley-data/elliptic_search_approx_scalarized'
   sim.data_dir+= '_muhat_'+str(sim.mu_hat)
   sim.data_dir+= '_lahat_'+str(sim.la_hat)

   gbc2_range=[300.0,500.0]

   sim.search_for_elliptic(gbc2_range)
#=============================================================================
else:
   raise ValueError('run_type = '+str(sim.run_type)) 
