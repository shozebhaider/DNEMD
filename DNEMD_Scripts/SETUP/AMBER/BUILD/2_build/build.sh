# Set constraints on protein
vmd -dispdev text -f x.prmtop x.pdb -e set_constraints.vmd
# Get cell dimensions and write xsc file
vmd -dispdev text -f x.prmtop x.pdb -e get_celldimension.vmd

cp x.pdb ../3_equil/
cp structure* ../3_equil/
cp x.prmtop ../3_equil/
cp input.xsc ../3_equil/
cp x.inpcrd ../3_equil/

cd ../3_equil/

nedit input_equil

