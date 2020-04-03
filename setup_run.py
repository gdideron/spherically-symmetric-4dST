import sys, time
from sim_class import Sim
#############################################################################
args= sys.argv
##############################################################################
### input paramters: set these by hand 
##############################################################################
sim= Sim(args)
#-----------------------------------------------------------------------------
sim.compactification_length= float(100) 
#-----------------------------------------------------------------------------
sim.evolve_time=   float(200)  ### in units of initial black hole mass for ze field 
sim.num_saved_times= int(500)
sim.cfl= 0.25
#-----------------------------------------------------------------------------
### scalar field potentials
#-----------------------------------------------------------------------------
sim.mu_hat= 0.01
sim.la_hat= 0.2

sim.gbc1= 0.0
sim.gbc2= 300.0
#-----------------------------------------------------------------------------
sim.phi_r= 21  ### where measuring phi
#-----------------------------------------------------------------------------
### initial data
#-----------------------------------------------------------------------------
sim.bh_mass= float(10.0)
#-----------------------------------------------------------------------------
### for the noncompact scalar profile
sim.initial_data_type= str("scalarized_bh")
sim.charge_hat= float(0.0)
#-----------------------------------------------------------------------------
### for the Gaussian-like pulse
#sim.initial_data_type= str("bump")
#sim.initial_data_type= str("bump_with_bh")
sim.amp= 0#float(5.0e-3)
sim.r_l= 0#float(24.0)
sim.r_u= 0#float(32.0)
#-----------------------------------------------------------------------------
sim.nx= pow(2,11)+1 
#-----------------------------------------------------------------------------
sim.set_derived_params()
##############################################################################
### for slurm script
##############################################################################
sim.walltime= '48:00:00' ### (hh:mm:ss)
sim.memory=   '8' ### MB 
##############################################################################
if (sim.run_type == "basic_run"):
	sim.data_dir= '/mnt/grtheory/jripley-data/test'
	sim.launch()
##############################################################################
elif (sim.run_type == "scan"):
	mu_hat=0.1
#-----------------------------------------------------------------------------
	if (mu_hat==0.01):
		sim.data_dir+= '_muhat_0.01'
	elif (mu_hat==0.05):
		sim.data_dir+= '_muhat_0.05'

	elif (mu_hat==0.1):
		sim.data_dir+= '_muhat_0.1'
	else:
		pass
#-----------------------------------------------------------------------------
	for la_hat in [0,0.2,0.4,0.6,0.8]:
		for gbc2 in [300,305,310,315,317.5,320,322.5,325]:
			sim.mu_hat= mu_hat
			sim.la_hat= la_hat
			sim.gbc2=  gbc2 
			sim.set_derived_params()
			sim.launch()
			time.sleep(5)
##############################################################################
elif (sim.run_type == "convergence_test"):
	num_res= int(input("number of resolutions "))
	sim.data_dir= '/mnt/grtheory/jripley-data/convergence_test'

	sim.charge_hat=0.05

	sim.mu_hat=0.0
	sim.la_hat=0.0

	sim.gbc2= 278 

	sim.set_derived_params()
	for i in range(num_res):	
		sim.launch()
		sim.nx= 2*(sim.nx-1)+1 
		sim.set_derived_params()
		time.sleep(5)
###############################################################################
### varying eta with a fixed initial phi amplitude
###############################################################################
elif (sim.run_type == "search_for_elliptic"):
	sim.data_dir= '/mnt/grtheory/jripley-data/elliptic_search_approx_scalarized'

	sim.charge_hat= 0.05

	sim.mu_hat= 0.05
	sim.la_hat= 3.2

	gbc2_range=[100.0,330.0]

	sim.search_for_elliptic(gbc2_range)
##############################################################################
else:
	raise ValueError("run_type = "+str(sim.run_type)) 
