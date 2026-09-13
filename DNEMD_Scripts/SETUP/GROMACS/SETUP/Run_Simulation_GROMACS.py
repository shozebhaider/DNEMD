# === GROMACS SIMULATIONS SET-UP ===

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

#!/usr/bin/env python3
"""
GROMACS batch runner with automatic input file generation.

Usage:
    python Run_Simulation_GROMACS.py --first-len 50 --interval 5 --n-intervals 10 --ligand-name XYZ
"""

import argparse
import glob
import os
import shutil
import subprocess
from pathlib import Path
from MDAnalysis import Universe

# ------------------------------
# Generate input ranges
# ------------------------------
def generate_input_ranges(first_len=50, interval=5, n_intervals=40):
    """
    Generate ranges like '0-50', '51-55', '56-60', ...
    """
    ranges = []
    start = 0
    end = first_len
    ranges.append(f"{start}-{end}")
    start = end + 1

    for _ in range(1, n_intervals):
        end = start + interval - 1
        ranges.append(f"{start}-{end}")
        start = end + 1

    return ranges

# ------------------------------
# Generate GROMACS .mdp input files
# ------------------------------

def generate_mdp_files(ranges, first_len=50, interval=5, base_dir="input_files", ligand_name="LIG"):
    """
    Generate .mdp files for each run length.
    First file uses first_len ns, subsequent ones increase by interval.
    """

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# Please change the .mdp parameters accordingly to your system. The provided example has been used to simulate a protein-ligand complex with the AMBER99SB-ILDN Force Fiels.
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

    mdp_template = """; GROMACS MD input parameters
integrator             = md
dt                     = 0.002
nsteps                 = {nsteps}       ; 1 ns = 500000 steps with 2 fs timestep
nstcomm                =  5000  ;  equal to nstenergy
nstxout-compressed     = 5000  ; writing every 5 ps
compressed-x-precision = 5000  ;  writing every 5 ps
nstenergy              = 5000  ;  writing every 5 ps
nstlog                 = 5000  ;  writing every 5 ps

;neighbor searching
nstlist             =  40 ; until 40 for gpu
ns_type             =  grid
cutoff-scheme       =  Verlet

;PME
rlist               = 1.2
coulombtype         = PME
rcoulomb            = 1.2
fourierspacing      = 0.16
pme_order           = 4
vdwtype             = cut-off
rvdw                = 1.2
DispCorr            = EnerPres

; temperature coupling is on in two groups
Tcoupl              =  V-rescale
tc-grps             =  Protein_{ligand_name}   water_and_ions
tau_t               =    0.1               0.1
ref_t               =    310               310

; Energy monitoring
energygrps          =  System

; Isotropic pressure coupling is now off
Pcoupl              =  Parrinello-Rahman
Pcoupltype          =  isotropic
tau_p               =  1.6
compressibility     =  4.5e-5
ref_p               =  1.0
refcoord-scaling    =  com

; Continue run without generating velocities
gen_vel             =  no
continuation        =  yes

constraint_algorithm=lincs
lincs_order         = 4
lincs-warnangle     = 90
constraints         =  h-bonds
"""

    os.makedirs(base_dir, exist_ok=True)

    def ns_to_steps(ns):
        # 1 ns = 500000 steps for 2 fs timestep
        return int(ns * 500000)

    filename = Path(base_dir) / f"input_{ranges[0]}.mdp"
    with open(filename, "w") as f:
        f.write(mdp_template.format(nsteps=ns_to_steps(first_len), ligand_name=ligand_name))
    print(f"[GEN] {filename} (run {first_len} ns)")

    runlen = first_len
    for r in ranges[1:]:
        runlen += interval
        filename = Path(base_dir) / f"input_{r}.mdp"
        with open(filename, "w") as f:
            f.write(mdp_template.format(nsteps=ns_to_steps(runlen), ligand_name=ligand_name))
        print(f"[GEN] {filename} (run {runlen} ns)")

# ------------------------------
# Manipulate coordinates (remove ligand)
# ------------------------------
def manipulate_structure(simdir, ligand_name):
    simdir = Path(simdir)
    gro_file = simdir / "confout.gro"
    tpr_file = simdir / "topol.tpr"

    u = Universe(str(tpr_file), str(gro_file))
    protein = u.select_atoms(f"not resname {ligand_name}")
    protein.write(str(simdir / "confout_nolig.gro"))
    print(f"[MDAnalysis] Removed residue {ligand_name} in {simdir}")

# ------------------------------
# Run commands with logging
# ------------------------------
def run(cmd, cwd=None, logfile=None, dry=False):
    msg = f"[RUN]{' (dry)' if dry else ''} cwd={cwd or os.getcwd()} :: " + " ".join(cmd)
    print(msg)
    if dry:
        return
    if logfile:
        with open(Path(cwd or ".") / logfile, "w") as fh:
            subprocess.run(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT, check=True)
    else:
        subprocess.run(cmd, cwd=cwd, check=True)

# ------------------------------
# Copy utility
# ------------------------------
def copy_globs(patterns, dest_dir, cwd=".", dry=False):
    dest = Path(dest_dir)
    if not dry:
        dest.mkdir(parents=True, exist_ok=True)
    for pat in patterns:
        matches = glob.glob(str(Path(cwd) / pat))
        if not matches:
            print(f"[WARN] No matches for pattern: {pat} in {cwd}")
        for src in matches:
            print(f"[CP]{' (dry)' if dry else ''} {src} -> {dest}/")
            if not dry:
                if os.path.isdir(src):
                    shutil.copytree(src, dest / os.path.basename(src), dirs_exist_ok=True)
                else:
                    shutil.copy2(src, dest)

# ------------------------------
# Main pipeline
# ------------------------------
def main():
    ap = argparse.ArgumentParser(description="Python wrapper for GROMACS batch runs")
    ap.add_argument("--first-len", type=int, default=50, help="First run length (ns)")
    ap.add_argument("--interval", type=int, default=5, help="Subsequent interval (ns)")
    ap.add_argument("--n-intervals", type=int, default=40, help="Total number of runs")
    ap.add_argument("--ligand-name", type=str, required=True, help="Unique residue name of the ligand to remove")
    ap.add_argument("--skip-equil", action="store_true", help="Skip Stage 1 (equilibrium runs)")
    ap.add_argument("--skip-noneq", action="store_true", help="Skip Stage 2 (nonequilibrium restarts)")
    ap.add_argument("--dry-run", action="store_true", help="Print actions only, do not execute")

    args = ap.parse_args()
    ranges = generate_input_ranges(args.first_len, args.interval, args.n_intervals)
    generate_mdp_files(ranges, args.first_len, args.interval, ligand_name=args.ligand_name)

    # === Stage 1: equilibrium GROMACS runs ===
    if not args.skip_equil:
        for idx, r in enumerate(ranges, start=1):
            print(f"\n=== Stage 1: Run {idx} ({r}) ===")
            mdp_file = f"input_files/input_{r}.mdp"

            # Preprocessing
            run(["gmx", "grompp", "-f", mdp_file, "-c", "structure.gro",
                 "-p", "topol.top", "-o", "md.tpr"], logfile="grompp.log", dry=args.dry_run)

            # MD run
            run(["gmx", "mdrun", "-deffnm", "md", "-v"], logfile="mdrun.log", dry=args.dry_run)

            # Save outputs
            copy_globs(["md.*", "confout.gro", "ener.edr", "md.log"], str(idx), dry=args.dry_run)

        print("\nEND OF EQUILIBRIUM SIMULATIONS")

    # === Stage 2: nonequilibrium restarts ===
    if not args.skip_noneq:
        for i in range(1, args.n_intervals + 1):
            simdir = Path(str(i))
            print(f"\n=== Stage 2: Nonequilibrium in {simdir} ===")

            # Copy parent files
            for fname in ["topol.top", "posre.itp"]:
                if Path(fname).exists():
                    shutil.copy2(fname, simdir / fname)

            # Remove ligand
            if not args.dry_run:
                manipulate_structure(simdir, args.ligand_name)

            # Prepare and run
            run(["gmx", "grompp", "-f", "input_files/input_0-50.mdp",
                 "-c", "confout_nolig.gro", "-p", "topol.top", "-o", "noneq.tpr"],
                cwd=simdir, logfile="grompp_noneq.log", dry=args.dry_run)

            run(["gmx", "mdrun", "-deffnm", "noneq", "-v"],
                cwd=simdir, logfile="mdrun_noneq.log", dry=args.dry_run)

            print(f"END OF NONEQ SIMULATION {i}")

if __name__ == "__main__":
    main()
