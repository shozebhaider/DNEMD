#!/bin/bash

device=0
input_dir="input_files"
output_prefix="output"
structure_files="structure.*"
log_file="MD.log"

# First run (0-50)
acemd3 --device $device ${input_dir}/input_0-50_MD1 > $log_file
mkdir -p 1
cp ${output_prefix}.* $structure_files $log_file 1/

# Subsequent runs (51-500 in steps of 10)
for i in {51..500..10}; do
  step=$(( (i - 51) / 10 + 2 ))
  input_file="${input_dir}/input_${i}-$(($i + 9))_MD1"
  output_dir="$step"
  
  acemd3 --device $device $input_file > $log_file
  mkdir -p $output_dir
  cp ${output_prefix}.* $structure_files $log_file $output_dir/
done

echo "END OF SIMULATIONS"
