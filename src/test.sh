#!/bin/bash
set -Ceu

# Bcl-xl is a test data
python mdanalyzer.py -r ../../00_samples/sample2/ref.pdb \
                     -i ../../00_samples/sample2/traj_aligned.xtc \
                     --rmsd_ref /Users/siida/Dropbox/01code/software/myTools/00_samples/sample2/em2.gro 

python plot.py