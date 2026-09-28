echo "3- Determining the average deviations at "

mkdir -p average_ChA/

sample_size=900  # Change the sample size here
# time_points={0:'0.1 ns', 4:'0.5 ns', 9:' 1 ns', 29: '3 ns' , 49: '5 ns' , 99:'10 ns'}

for tframe in  0 4 9 29 49 99
do
    echo "                                         ${tframe} ps "

    touch average_ChA/temp

    for rep in REP{01..20}  # Adjusted for your replicates
    do
        for time in {500..4900..100} 
        do
            for perturbed_index in {1..45}  # Loop over each perturbed simulation
            do
                rmsdev_file="${rep}/rms_ChA/${time}frame/rmsdev_${tframe}_${perturbed_index}_ChA.xvg"

                # Ensure the file exists before trying to process it
                if [[ -f "$rmsdev_file" ]]; then
                    egrep -v "#|@" $rmsdev_file | awk '{print $2}' > average_ChA/aux
                    paste average_ChA/temp average_ChA/aux > average_ChA/temp2
                    mv average_ChA/temp2 average_ChA/temp
                fi
            done
        done
    done

    # Calculating the global average
    awk -v sample_size=$sample_size '{for (i=1; i<=NF; i++) sum+=$i; printf "%10.3f%10i\n", sum/NF, sample_size; sum=0}' average_ChA/temp > average_ChA/avg_temp

    # Assuming residue number is in the first column of rmsdev files
    egrep -v "#|@" ${rep}/rms_ChA/500frame/rmsdev_${tframe}_0_ChA.xvg | awk '{printf "%10i\n", $1}' > average_ChA/res

    paste average_ChA/res average_ChA/avg_temp > average_ChA/average_${tframe}

    # Assuming 45 is the correct sample size
    paste average_ChA/temp average_ChA/average_${tframe} > average_ChA/stats

    echo "#       CA    Average     SD        Sample     SE     SE
          #                                    Size     (1SD)  (2SD)
###############################################" > average_ChA/SD_${tframe}

    awk -v sample_size=$sample_size '{for (i=1; i<=sample_size; i++) sum+=(($i-$48)*($i-$48)); printf "%10i%10.3f%10.3f%10i%10.3f%10.3f\n", $47,$48, sqrt((sum/sample_size)), $49, (sqrt(sum/sample_size))/(sqrt($49)), (2*sqrt(sum/sample_size))/(sqrt($49)) ; sum=0}' average_ChA/stats >> average_ChA/SD_${tframe}

    # Housecleaning
    rm -f average_ChA/res average_ChA/avg_temp average_ChA/temp average_ChA/stats
done

