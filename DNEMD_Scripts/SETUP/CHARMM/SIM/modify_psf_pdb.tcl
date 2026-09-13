# Ensure that the VMD console is clear
package require psfgen

# Reset psfgen and load the psf and pdb files
resetpsf
readpsf structure.psf
coordpdb structure.pdb

# Load the psf and pdb files into VMD
mol load psf structure.psf pdb structure.pdb

# Delete the segid corresponding to the cholesterol molecules in the allosteric site
delatom P1
delatom P3
# Have psfgen write out the new psf and pdb files
writepsf x_out.psf
writepdb x_out.pdb

# Print a message to confirm completion
puts "PSF and PDB files have been written to x_out.psf and x_out.pdb, respectively."
