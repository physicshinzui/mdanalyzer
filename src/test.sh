#!/bin/bash
set -Ceu

python mdanalyzer.py -r ../../00_samples/sample2/ref.pdb \
                     -i ../../00_samples/sample2/traj_aligned.xtc \