This is a protocol to set up eq and non-eq simulations. 
Set up the simulations in the same regular way using HTMD as using the jupyter notebooks in the BUILD folder. 

Make sure that you have all the input files in the cwd that can run short bursts of MD simulation in bins of 5ns each. 

First run the 1_eq_run.sh file.

At the end of the 50ns simulation, it will make a directory and the copy all the files in the directory.
The second simulation will begin from 50ns onwards for the next 10ns and so on till the end of the run at 500ns.

Now each directory will contain the conventional files
structure.* output.* restart.* 

The input into the non-eq simulations will be:
the perturned x_out.coor, x_out.vel 
x_out.psf and x_out.pdb
The size of the cell will come from output.xsc
The parameters will come from the same parameters file used to run the eq simulations

essential files to be here in this folder:

input.xsc/input.coor. These files are from output.xsc/output.coor from equilibration

x_out.pdb  without ligand...edit it using psfgen in vmd
x_out.psf without ligand...edit it using using psfgen in vmd

vmd -dispdev text -e modify_psf_pdb.tcl 

Then run the 2_non-eq_run.sh

Dont forget to change the number of device in 1*.sh and 2*.sh

