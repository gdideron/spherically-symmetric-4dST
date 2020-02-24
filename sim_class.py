##############################################################################
import subprocess, os, time
from typing import List 

home_dir= '/home/jripley/codes_general_edgb_feynman'
data_dir= '/mnt/grtheory/jripley-data/general_edgb'

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
	def write_sim_params(self)->None:
		with open(self.output_dir+'/sim_params.txt','w') as f:
			attrs= vars(self)
			for param in attrs:
				f.write('{} {}\n'.format(param,attrs[param]))	
##############################################################################
	def make_output_dir(self)->None:
		self.output_dir= str(
			data_dir+'/'
		+	'_'.join('_'.join(time.asctime().split(' ')).split(':'))
		+	'_nx'+str(self.nx)
		+	'_mu'+str(self.mu)
		+	'_la'+str(self.la)
		+	'_gbc2'+str(self.gbc2)
		)
		os.makedirs(self.output_dir)
##############################################################################
	def write_slurm_script(self) -> None:
		with open('{}run.slurm'.format(home_dir), 'w') as f:
			f.write('#!/bin/sh\n')
			f.write('#SBATCH -N 1\t\t# nodes=1\n')
			f.write('#SBATCH --ntasks-per-node=1\t\t# ppn=1\n')
			f.write('#SBATCH -J g{:.2f}\t\t# job name\n'.format("edgb")
			f.write('#SBATCH -t {}\t\t# walltime (dd:hh:mm:ss)\n'.format(self.walltime)
			f.write('#SBATCH -p dept\t\t# partition/queue name\n')
			f.write('#SBATCH --mem={}MB\t\t# memory in MB\n'.format(self.memory)
			f.write('#SBATCH --output={}\t\t# file for STDOUT\n'.format(self.output_dir)
			f.write('#SBATCH --mail-user=jripley@princeton.edu\t\t# Mail  id of the user\n')
	#		f.write('#SBATCH --mail-type=begin\t\t# Slurm will send mail at the beginning of the job\n')
	#		f.write('#SBATCH --mail-type=end\t\t# Slurm will send at the completion of your job\n')
			f.write('\n./collapse_excise_elliptic\n\n')

		shutil.copyfile(
			'{}run_TEdGB_collapse.slurm'.format(JRIPLEY_DIR),
			'{}run_TEdGB_collapse.slurm'.format(file_data['output_dir'])
		)
			
		return 	
