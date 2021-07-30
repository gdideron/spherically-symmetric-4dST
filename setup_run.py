#!/usr/bin/env python
#=============================================================================
##
## Parameter file
## Usage: 
## ./setup [run_type] (e.g. ./setup basic_run)
##
import sys, time, os
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
sim.evolve_time=   float(200)  ### in units of initial black hole mass for ze field 
sim.num_saved_times= int(500)
sim.cfl= 0.2
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
sim.Be_2=  25.0
sim.Be_3=  0.0
sim.Be_4=  -7.5

sim.Be_exp2= 0.0
#-----------------------------------------------------------------------------
sim.phi_r= 15  ### where measuring phi (radial distance)
#-----------------------------------------------------------------------------
## excision point as ratio of apparent horizon
sim.excision_ratio = 0.95
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

sim.amp= float(0.01)
sim.r_l= float(12.0)
sim.r_u= float(26.0)
#-----------------------------------------------------------------------------
sim.nx= pow(2,12)+1 
#-----------------------------------------------------------------------------
sim.set_derived_params()
#=============================================================================
## for slurm script
#=============================================================================
sim.walltime= '96:00:00' ### (hh:mm:ss)

project_dir= os.path.dirname(os.path.realpath(__file__))
sim.data_dir= os.path.join(project_dir,"output")
#sim.data_dir= '/tigress/jripley/edgb'

sim.slurm= False
#=============================================================================
if (sim.run_type == 'basic_run'):
   sim.launch()
#=============================================================================
elif (sim.run_type == 'ramp'):
   sim.data_dir= '/tigress/jripley/edgb/ramp_Be_2'
   for be2 in range(18,23,1):
      sim.Be_2= be2
      sim.Be_4= -4*be2
      sim.launch()
      time.sleep(1)
#=============================================================================
elif (sim.run_type == 'scan'):
   sim.data_dir= '/tigress/jripley/edgb/scan_Be_4_large'
   sim.Be_2= 25
   for be4 in range(-60*sim.Be_2-810000,-60*sim.Be_2-810000+5000,500):
      sim.Be_4= be4
      sim.launch()
      time.sleep(5)
#=============================================================================
elif (sim.run_type == 'res_study'):
   num_res= int(input('number of resolutions '))
   for i in range(num_res):	
      sim.launch()
      sim.nx= 2*(sim.nx-1)+1 
#=============================================================================
### varying eta with a fixed initial phi amplitude
#=============================================================================
elif (sim.run_type == 'elliptic_search'):

#   sim.data_dir= '/tigress/jripley/edgb/elliptic_search_Be_exp2'
   sim.data_dir= '/tigress/jripley/edgb/elliptic_search_Be_4'

   coupling_range= [-60*sim.Be_2, 0]

   sim.search_for_elliptic(coupling_range,'Be_4')
#=============================================================================
else:
   raise ValueError('run_type = '+str(sim.run_type)) 
