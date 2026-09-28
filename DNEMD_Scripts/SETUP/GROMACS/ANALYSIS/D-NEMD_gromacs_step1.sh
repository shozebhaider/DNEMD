grmsf="gmx rmsf"

echo "1- Determining the Calpha Deviations at: "

for rep in REP{01..20}
do
    echo "     $rep"

    perturbed_index=1  # Start with the first perturbed index

    for time in {500..4900..100}
    do
        echo "                               ${time} ns"

        # time_points={0:'0.1 ns', 4:'0.5 ns', 9:' 1 ns', 29: '3 ns' , 49: '5 ns' , 99:'10 ns'}
        for tframe in 0 4 9 29 49 99
        do
            mkdir -p  $rep/rms_ChA/${time}frame/

            # Calculate the RMSD between pairs of perturbed and unperturbed simulations at a given time
            ${grmsf} -f  $rep/lig/${time}frame/conf_${tframe}_ChA.pdb \
                     -s  $rep/apo/${perturbed_index}/conf_${tframe}_ChA.pdb \
                     -o  $rep/rms_ChA/${time}frame/rms_${tframe}_${perturbed_index}_ChA.xvg \
                     -od $rep/rms_ChA/${time}frame/rmsdev_${tframe}_${perturbed_index}_ChA.xvg \
		     -res \
                     -fit   <<EOF   
C-alpha
EOF

            # Housecleaning: Remove temporary files if needed
            rm -f  $rep/rms_ChA/${time}frame/rms_${tframe}_${perturbed_index}_ChA.xvg 
        done 

        # Increment the perturbed index for the next iteration
        ((perturbed_index++))
    done
done
