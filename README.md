# Evolution codes for EdGB gravity in spherical symmetry.

For help please contact jripley [at] princeton [dot] edu. 

This setup files (setup.py and sim_class.py)
are configured to run on the Feynman cluster at Princeton University.
They write and run a slurm script
(see https://slurm.schedmd.com/documentation.html),
which launches the 'run' binary (see the Makefile).
The run binary also can be run on a personal computer.

The code ouputs to csv files.
There is also a .cpp and .hpp file for saving to .sdf files
using the  bbhutil library; see 
http://laplace.physics.ubc.ca/Group/Software.html
for more information. Remove the bbhutil headers if you do not want to save
to .sdf.

# Citation

If you use this code in any way, please cite the
original paper that introduces it:
 
https://inspirehep.net/literature/1795898 
