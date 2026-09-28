echo "2- Determining the global average structure at "

mkdir -p global_average/apo/ global_average/lig/

for tframe in  0 100 5000 10000
do

    echo "                                      ${tframe}  ps"

    cat REP*/lig/*ns/conf_${tframe}_ChA.pdb > all_protein-lig_${tframe}.pdb
    cat REP*/apo/*/conf_${tframe}_ChA.pdb > all_protein_${tframe}.pdb


    gmx make_ndx -f all_protein-lig_${tframe}.pdb  \
                -o index_protein-lig.ndx  <<EOF
quit
EOF
	   
    echo "                                    Ligand-bound average"

    gmx rmsf   -f all_protein-lig_${tframe}.pdb  \
               -o rmsf_protein-lig_${tframe}ns.xvg \
               -s REP01/lig/50ns/conf_${tframe}_ChA.pdb  \
               -n  index_protein-lig.ndx \
               -ox global_average/lig/average_${tframe}ns.pdb \
               -oq global_average/lig/bfactor_${tframe}ns.pdb \
               -fit   <<EOF  
Protein         
EOF


    gmx make_ndx -f all_protein_${tframe}.pdb  \
                -o index_protein.ndx  <<EOF
quit
EOF
	   
    echo "                                    Apo average"

    gmx rmsf   -f all_protein_${tframe}.pdb  \
               -o rmsf_protein_${tframe}ns.xvg \
               -s REP01/apo/1/conf_${tframe}_ChA.pdb  \
               -n  index_protein.ndx \
               -ox global_average/apo/average_${tframe}ns.pdb \
               -oq global_average/apo/bfactor_${tframe}ns.pdb \
               -fit  <<EOF  
Protein         
EOF

    #housecleaning
    rm -f global_average/lig/bfactor_${tframe}ns.pdb rmsf_protein-lig_${tframe}ns.xvg  rmsf_protein_${tframe}ns.xvg global_average/apo/bfactor_${tframe}ns.pdb all_protein_${tframe}.pdb all_protein-lig_${tframe}.pdb index_protein.ndx index_protein-lig.ndx r*/lig/*ns/protein_${tframe}.pdb r*/apo/*ns/protein_${tframe}.pdb  

done
