##############################################################################
import subprocess, os, time
from typing import List 
##############################################################################
class Sim:
##############################################################################
	def __init__(self,args:List[str])->None:
		self.home_dir= str(os.getcwd())

		assert len(args) > 1, (
			'argv[1] is empty-meed a run_type to run!'
		)	
		self.run_type = args[1]
		if (len(args)>2 and args[2]=='debug'):
			self.debug=True
		else:
			self.debug=False
##############################################################################
	def set_derived_params(self)->None:
		self.dx= float(
			self.compactification_length/(self.nx-1)
		)
		self.dt= float(
			self.cfl*self.dx
		)
		self.nt= int(
			self.bh_mass*self.evolve_time/self.dt
		)
		self.t_step_save= int(
			self.nt/float(self.num_saved_times)
		)	
		self.initial_exc_i= 0
		if (self.initial_data_type=="bump_with_bh"):
			self.initial_exc_i= int(
			(
				(1.8*self.bh_mass)/(1+(1.8*self.bh_mass/self.compactification_length))
			)/self.dx
			)
		if (self.t_step_save==0):
			self.t_step_save= 1
		self.phi_pt= int(
			(	(self.phi_r)
			/	(1+(self.phi_r/self.compactification_length))
			)/self.dx
		)
##############################################################################
	def make_output_dir(self)->None:
		self.output_dir= str(
			'_'.join('_'.join(time.asctime().split(' ')).split(':'))
		+	'_nx'+str(self.nx)
		)
		os.makedirs(self.output_dir)
##############################################################################
	def write_sim_params(self)->None:
		with open(self.output_dir+'/sim_params.txt','w') as f:
			attrs= vars(self)
			for param in attrs:
				f.write('{} {}\n'.format(param,attrs[param]))	
##############################################################################
	def launch_run(self)->None:
		output_file= self.output_dir+'/output.txt'
		run_str= (
			'./run '+self.home_dir+'/'+self.output_dir
		+	'>'+output_file+' 2>&1 &'
		)
		if (self.debug):
			run_str= 'valgrind -v '+run_str
		subprocess.call(run_str,shell='True') 
