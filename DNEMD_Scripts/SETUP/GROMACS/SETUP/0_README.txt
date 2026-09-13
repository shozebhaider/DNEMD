This is a protocol to run eq and non-eq simulations. 


Make sure that you have all the input files in the cwd that can run short bursts of MD simulation in bins of 50ns each. 

usage: Run_Simulation_GROMACS.py [options]

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
python Run_Simulation_GROMACS.py --first-len 5 --interval 5 --n-intervals 3 --ligand-resid 264

