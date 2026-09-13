This is a protocol to set up eq and non-eq simulations. 
Set up the simulations in the same regular way using Xleap. 

Make sure that you have all the input files in the cwd that can run short bursts of MD simulation in bins of 50ns each. 

usage: Run_Simulation.py [options]

options:
  --device <id>          GPU device ID for ACEMD (default: 0)
  --first-len <ns>       Length of the first run in ns (default: 50)
  --interval <ns>        Interval added for subsequent runs (default: 5)
  --n-intervals <n>      Total number of runs including the first (default: 40)
  --ligand-resid <id>    Residue index of the ligand to remove (required)
  --skip-equil           Skip Stage 1 (equilibrium runs)
  --skip-noneq           Skip Stage 2 (non-equilibrium runs)
  --dry-run              Print commands without executing

example usage:
python Run_Simulation.py --first-len 5 --interval 5 --n-intervals 3 --ligand-resid 264

Now each directory will contain the conventional files
x.* output.* restart.* 

The input into the non-eq simulations will be the perturned x_out.coor file
x_out.prmtop and x_out.pdb
The size of the cell will come from output.xsc

essential files to be here in this folder:

input.xsc/input.coor. These files are from output.xsc/output.coor from equilibration
x_out.pdb  without ligand...edit it directly
x_out.prmtop without ligand...edit it using parmed

python parmed x.prmtop
strip :262
writeout x_out.prmtop
quit

