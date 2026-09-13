#!/usr/bin/env python3
"""
ACEMD batch runner with automatic input file generation.

Usage:
    python run_acemd_pipeline.py --first-len 50 --interval 5 --n-intervals 10
"""

import argparse
import glob
import os
import shutil
import subprocess
from pathlib import Path
from htmd.ui import Molecule

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
# Generate input files
# ------------------------------
def generate_input_files(ranges, first_len=50, interval=5, base_dir="input_files"):
    """
    Generate ACEMD input files for each range.
    First file uses 'run first_len ns', subsequent files increase by interval.
    """

    template = """bincoordinates                 input.coor
coordinates                    structure.pdb
cutoff                         9
extendedsystem                 input.xsc
langevin                       on
langevindamping                0.1
langevintemp                   300
parmfile                       x.prmtop
pme                            on
restart                        on
structure                      x.prmtop
switchdist                     7.5
switching                      on
temperature                    300
timestep                       4
xtcfile                        trajectory.xtc
xtcfreq                        25000

useflexiblecell                off
useconstantratio               off

run                            {run_length}ns
"""

    os.makedirs(base_dir, exist_ok=True)

    # first run
    filename = Path(base_dir) / f"input_{ranges[0]}_MD1"
    with open(filename, "w") as f:
        f.write(template.format(run_length=first_len))
    print(f"[GEN] {filename} (run {first_len} ns)")

    # subsequent runs
    runlen = first_len
    for r in ranges[1:]:
        runlen += interval
        filename = Path(base_dir) / f"input_{r}_MD1"
        with open(filename, "w") as f:
            f.write(template.format(run_length=runlen))
        print(f"[GEN] {filename} (run {runlen} ns)")

def manipulate_coor_and_vel(simdir, ligand_resid):
    simdir = Path(simdir)

    # coordinates
    m = Molecule(str(simdir / 'x.pdb'))
    m.read(str(simdir / "output.coor"))
    m.filter(f"not resid {ligand_resid}")
    m.write(str(simdir / "x_out.coor"), type="coor")
    
    #v = Molecule('x.pdb')
    #v.read('output.vel', type='coor')
    #v.filter(f"not resid {ligand_resid}")
    #v.write('x_out.vel', type='coor')

    print(f"[HTMD] Coordinate manipulation done in {simdir} (removed resid {ligand_resid})")

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
        subprocess.run(cmd, cwd=cwd, check=True, shell=True)

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
    ap = argparse.ArgumentParser(description="Python wrapper for ACEMD batch runs")
    ap.add_argument("--device", default="0", help="ACEMD GPU device id (default: 0)")
    ap.add_argument("--first-len", type=int, default=50, help="First run length (ns)")
    ap.add_argument("--interval", type=int, default=5, help="Subsequent interval (ns)")
    ap.add_argument("--n-intervals", type=int, default=40, help="Total number of runs (including first)")
    ap.add_argument("--skip-equil", action="store_true", help="Skip Stage 1 (initial ACEMD runs)")
    ap.add_argument("--skip-noneq", action="store_true", help="Skip Stage 2 (nonequilibrium restarts)")
    ap.add_argument("--dry-run", action="store_true", help="Print actions only, do not execute")
    ap.add_argument("--ligand-resid", type=int, required=True, help="Residue index of the ligand to remove")

    args = ap.parse_args()

    ranges = generate_input_ranges(args.first_len, args.interval, args.n_intervals)

    # generate input files automatically
    generate_input_files(ranges, args.first_len, args.interval, base_dir="input_files")

    # === Stage 1: initial ACEMD runs ===
    if not args.skip_equil:
        for idx, r in enumerate(ranges, start=1):
            print(f"\n=== Stage 1: Run {idx}  (input_files/input_{r}_MD1) ===")
            run(
                ["acemd", f"input_files/input_{r}_MD1"],
                cwd=".",
                logfile="MD.log",
                dry=args.dry_run
            )
            target = str(idx)
            copy_globs(["output.*", "x.*", "structure.*", "MD.log"], target, cwd=".", dry=args.dry_run)

        print("\nEND OF SIMULATIONS")

    # === Stage 2: nonequilibrium restarts ===
    if not args.skip_noneq:
        parent_files = ["x_out.prmtop", "x_out.pdb", "input-noneqb_MD1"]
        for i in range(1, args.n_intervals + 1):
            simdir = Path(str(i))
            print(f"\n=== Stage 2: Nonequilibrium in {simdir} ===")
            for fname in parent_files:
                src = Path(fname)
                dst = simdir / fname
                print(f"[CP]{' (dry)' if args.dry_run else ''} {src} -> {dst}")
                if not args.dry_run:
                    dst.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(src, dst)

            # Copy output.coor -> start_frame.coor
            src_coor = simdir / "output.coor"
            dst_coor = simdir / "start_frame.coor"
            print(f"[CP]{' (dry)' if args.dry_run else ''} {src_coor} -> {dst_coor}")
            if not args.dry_run:
                shutil.copy2(src_coor, dst_coor)

            if not args.dry_run:
                manipulate_coor_and_vel(simdir, args.ligand_resid)
            run(["acemd", "input-noneqb_MD1"],
                cwd=simdir, logfile="MD.log", dry=args.dry_run)

            print(f"END OF MD SIMULATION {i}")

if __name__ == "__main__":
    main()
