#!/bin/bash

device=0
psf_file="x_out.psf"
pdb_file="x_out.pdb"
velocity_coor_manipulation_script="velocity_coor_manipulation.py"
input_file="input-noneqb_MD1"
parameters_file="parameters"

for step in {1..46}; do
  echo "Starting Non-EQ Simulation for Step $step"
  cd $step || { echo "Failed to enter directory $step. Exiting."; exit 1; }
  
  cp ../$psf_file ../$pdb_file ../$velocity_coor_manipulation_script ../$input_file ../$parameters_file ./ || { echo "Failed to copy files. Exiting."; exit 1; }
  cp output.coor start_frame.coor || { echo "Failed to copy output.coor to start_frame.coor. Exiting."; exit 1; }
  
  python $velocity_coor_manipulation_script || { echo "Python script failed. Exiting."; exit 1; }
  
  acemd3 --device $device $input_file > MD.log || { echo "ACEMD3 failed. Exiting."; exit 1; }
  
  echo "END OF MD SIMULATION $step"
  cd ..
done
