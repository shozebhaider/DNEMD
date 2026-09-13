# === D-NEMD ANALYSIS SCRIPT ===

# === LICENCE ===
# This project is licensed under the MIT License. 

# MIT License
#
# Copyright (c) 2025 University of Bristol
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


# === OVERVIEW ===
# This Python script allows the user to calculate the time evolution of a property—specifically, the displacement of alpha carbon atoms—across a series of time points. At each time point, the system’s response to the applied perturbation is averaged over the non-equilibrium ensemble, and the standard error of the 
# mean is calculated to ensure statistical significance. 
# The script produces the following outputs:
# - PDB files (one for each time point) in which the displacement of each alpha carbon atom is stored as the B-factor in a common reference structure. These files are useful for quantifying the magnitude of the displacement. 
# - PDB files (one for each time point), each showing a different structure that reflects the system’s positional changes over time. These are useful for visualizing the directionality of the displacements.

# === USAGE ===
# - The script requires an updated version of python and MDAnalysis

# === AUTHORS === 
# Main author: 
# - Dan Pike: ze23981@bristol.ac.uk
# Other author:
# - Lorenzo Tulli: lorenzo.tulli@bristol.ac.uk
# - Michael Beer: michael.beer@bristol.ac.uk
# - Sofia Oliveira: sofia.oliveira@bristol.ac.uk
  
# === CONTACT ===
# If you have any questions or feedback, feel free to reach out:
# - Dan Pike: ze23981@bristol.ac.uk
# - Lorenzo Tulli: lorenzo.tulli@bristol.ac.uk
# - Dr Michael Beer: michael.beer@bristol.ac.uk
# - Dr Sofia Oliveira: sofia.oliveira@bristol.ac.uk

# Thank you for using this script! We're looking forward to receiving your feedback.

import os
import MDAnalysis as mda
from MDAnalysis.analysis import align
import numpy as np
import glob

# === CONFIGURATION SECTION ===
# Define base paths and files - MODIFY THESE
eq_base_path = "/path/to/EQ"  # Replace with actual path
neq_base_path = "/path/to/NEQ"  # Replace with actual path
null_base_path = None  # Set to None to skip null simulations or add path to Null simulations

# Topology file paths (different for EQ and NEQ)
eq_topology_file = "/path/to/EQ_topology"  # Topology for equilibrium simulations
neq_topology_file = "/path/to/EQ_topology"  # Topology for non-equilibrium simulations
# Note: NULL uses the same topology as EQ

# Path to base reference structure (for C-alpha alignment)
base_ref_path = "path/to/base_reference.pdb"  # Replace with actual path

# List of time points to analyze
time_points = [50, 100, 150, 200]  # Example values, modify as needed

# Number of EQ runs and times to include in analysis
runs = list(range(1, 11))  # Runs 1-10
times = list(range(50, 246, 5))  # Time points 50, 55, ..., 245

# Analysis settings
significance_threshold = 2.0  # Significance threshold in standard errors

# Processing options
run_extraction = True  # Set to False to skip extraction and only run analysis
 
# Output directory
output_base_dir = "extracted_frames"
if not os.path.exists(output_base_dir):
    os.makedirs(output_base_dir)

# Analysis output directory (with threshold in the name)
analysis_output_dir = f"analysis_results_threshold_{significance_threshold}"
if not os.path.exists(analysis_output_dir):
    os.makedirs(analysis_output_dir)

# === PART 1: EXTRACT FRAMES FROM TRAJECTORIES ===
if run_extraction:
    print("=== EXTRACTING FRAMES ===")

    # Define simulation types to process with their corresponding topology files
    sim_paths = {
        "EQ": {"path": eq_base_path, "topology": eq_topology_file},
        "NEQ": {"path": neq_base_path, "topology": neq_topology_file}
    }

    # Add NULL path if provided (using EQ topology)
    if null_base_path:
        sim_paths["NULL"] = {"path": null_base_path, "topology": eq_topology_file}

    # Loop through directories 1-10
    for i in range(1, 11):
        # Loop through subdirectories 50, 55, ..., 245
        for x in range(50, 246, 5):
            for sim_type, sim_info in sim_paths.items():
                # Construct paths for this simulation type
                run_dir = os.path.join(sim_info["path"], str(i))
                time_dir = os.path.join(run_dir, str(x))
                
                # Create output directory for this combination
                output_dir = os.path.join(output_base_dir, f"{sim_type}_run_{i}_{x}")
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                
                # Define paths to topology and trajectory files
                topology_path = sim_info["topology"]            
                trajectory_pattern = os.path.join(time_dir, f"*{x}*.nc") # Adjust according to trajectory file type
                trajectory_path = glob.glob(trajectory_pattern)[0]  # Use first matching file
                
                try:
                    # Load the trajectory
                    u = mda.Universe(topology_path, trajectory_path)
                    
                    # Define atoms to keep (by inverting the strips)
                    keep_atoms = u.select_atoms("not (resname WAT or resname Na+)") # Adjust as needed
                    
                    # Extract specific frames that match time_points
                    for frame_number in time_points:
                        # Switch to this frame
                        u.trajectory[frame_number]
                        
                        # Create a new universe with only the atoms we want to keep
                        temp_universe = mda.Merge(keep_atoms)
                        
                        # Save the frame as PDB with a unique name
                        output_filename = os.path.join(output_dir, 
                                                     f"{sim_type}_run{i}_{x}_{frame_number}.pdb")
                        temp_universe.atoms.write(output_filename)
                        print(f"Wrote {output_filename}")
                
                except Exception as e:
                    print(f"Error processing {sim_type} run {i}, trajectory {x}: {e}")
                    continue
else:
    print("=== SKIPPING EXTRACTION PHASE ===")

# === PART 2: ANALYZE EXTRACTED FRAMES ===
print("\n=== ANALYZING FRAMES ===")
print(f"Using significance threshold: {significance_threshold}")

# Load the base C-alpha structure 
try:
    base_ref = mda.Universe(base_ref_path)
    original_base_ca = base_ref.select_atoms('name CA')
    
    # Create a new Universe object with only the C-alpha atoms
    ca_universe = mda.Merge(original_base_ca)
    
    # Initialize tempfactors attribute
    ca_universe.add_TopologyAttr('tempfactors')
    
    for tp in time_points:
        print(f"\nAnalyzing time point {tp}")
        vectors = []
        
        for run in runs:
            for time in times:
                try:
                    # Construct paths to extracted PDB files
                    eq_path = os.path.join(output_base_dir, f"EQ_run_{run}_{time}", 
                                          f"EQ_run{run}_{time}_{tp}.pdb")
                    neq_path = os.path.join(output_base_dir, f"NEQ_run_{run}_{time}", 
                                           f"NEQ_run{run}_{time}_{tp}.pdb")
                    
                    # Only include null path if it exists
                    if null_base_path:
                        null_path = os.path.join(output_base_dir, f"NULL_run_{run}_{time}", 
                                               f"NULL_run{run}_{time}_{tp}.pdb")
                    else:
                        null_path = None
                    
                    # Load structures
                    ref = mda.Universe(eq_path)
                    other = mda.Universe(neq_path)
                    
                    # Select C-alpha atoms
                    ca_ref = ref.select_atoms('name CA').positions
                    ca_other = other.select_atoms('name CA').positions
                    
                    # Align structures
                    align.alignto(other, ref, select='name CA')
                    
                    if null_path:
                        null_pert = mda.Universe(null_path)
                        align.alignto(null_pert, ref, select='name CA')
                        ca_null_pert = null_pert.select_atoms('name CA').positions
                        
                        # Calculate differences
                        diff2 = ca_null_pert - ca_ref
                        diff1 = ca_other - ca_ref
                        diff = diff1 - diff2
                    else:
                        # If no null perturbation, just use the direct difference
                        diff = ca_other - ca_ref
                    
                    vectors.append(diff)
                    print(f"  Processed run {run}, time {time}")
                    
                except Exception as e:
                    print(f"  Error processing {run}/{time}: {e}")
                    continue
        
        if len(vectors) == 0:
            print(f"  No valid data for time point {tp}, skipping analysis")
            continue
        
        # Convert to numpy array
        vectors_array = np.array(vectors)
        
        # Calculate average vectors
        avg_vectors = np.mean(vectors_array, axis=0)
        
        # Calculate statistics
        avg_disp = np.linalg.norm(avg_vectors, axis=1)
        sd = np.std(vectors_array, axis=0)
        se = sd / np.sqrt(len(vectors))
        avg_x = avg_vectors[:, 0]
        avg_y = avg_vectors[:, 1]
        avg_z = avg_vectors[:, 2]
        
        # SE of magnitude (added small epsilon to avoid division by zero)
        se_mag = (1 / np.maximum(avg_disp, 1e-10)) * np.sqrt(
            (avg_x**2)*(se[:, 0]**2) + 
            (avg_y**2)*(se[:, 1]**2) + 
            (avg_z**2)*(se[:, 2]**2)
        )
        
        # Adjust vectors based on significance threshold
        adjusted_avg_disp = np.where(np.abs(avg_disp) - significance_threshold*se_mag > 0, avg_disp, 0)
        adjusted_avg_vectors = np.where(np.abs(avg_vectors) - significance_threshold*se[:, np.newaxis] > 0, avg_vectors, 0)
        
        # Write summary stats 
        stats_filename = os.path.join(analysis_output_dir, f"vec_norm_stats_{tp}.txt")
        with open(stats_filename, 'w') as f:
            # Write header
            f.write("#" * 128 + "\n")
            f.write(f"# SIGNIFICANCE THRESHOLD: {significance_threshold} standard errors\n")
            f.write("#                 Average Displacement                Sample Size                                   SE (1SD)                SE (1SD) of           SE (std)                  SE (std) of       #\n")
            f.write("#              --------------------------       ---------------------      Average        --------------------------          Average        --------------------------          Average       #\n")
            f.write("#       CA     x-axis    y-axis    z-axis       x         y         z    Displacement     x-axis    y-axis    z-axis    Displacement vector  x-axis    y-axis    z-axis    Displacement vector #\n")
            f.write("#" * 128 + "\n")
            # Write data
            for i in range(len(avg_vectors)):
                f.write(f"{i+1:>9} {avg_vectors[i][0]:>8.3f} {avg_vectors[i][1]:>8.3f} {avg_vectors[i][2]:>8.3f} ")
                f.write(f"{len(vectors):>8} {len(vectors):>8} {len(vectors):>8} {avg_disp[i]:>12.3f} ")
                f.write(f"{se[i][0]:>10.3f} {se[i][1]:>10.3f} {se[i][2]:>10.3f} {se_mag[i]:>12.3f} ")
                f.write(f"{significance_threshold*se[i][0]:>10.3f} {significance_threshold*se[i][1]:>10.3f} {significance_threshold*se[i][2]:>10.3f} {significance_threshold*se_mag[i]:>12.3f}\n")
        
        # Map the norms of the C-alpha vectors onto the C-alpha atoms in the new Universe object
        for atom, bfactor in zip(ca_universe.atoms, adjusted_avg_disp):
            atom.tempfactor = bfactor
        
        # Write the C-alpha-only structure with mapped B-factors
        pdb_filename = os.path.join(analysis_output_dir, f"vec_norm_{tp}.pdb")
        ca_universe.atoms.write(pdb_filename)
        
        # Update CÎ± positions with average displacement vectors
        base_ca = original_base_ca.copy()
        base_ca.positions = original_base_ca.positions + avg_vectors
        
        # Write to PDB
        pdb_filename = os.path.join(analysis_output_dir, f"avg_vec_{tp}.pdb")
        with mda.Writer(pdb_filename, multiframe=False) as pdb_writer:
            pdb_writer.write(base_ca)
        
        print(f"  Wrote analysis results to {analysis_output_dir}")
    
except Exception as e:
    print(f"Error in analysis: {e}")

print("\nProcessing complete!")