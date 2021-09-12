#!/usr/bin/env python3
import sys
from MDAnalysis import Universe
from tqdm import tqdm 

def parser():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("-f" , "--traj"    , required=True)
    p.add_argument("-s" , "--top"     , required=False)
    p.add_argument("-l", "--sele_list", required=True, nargs=2)
    args = p.parse_args()
    return args.traj, args.top, args.sele_list

def main():
    import numpy as np
    flog = open('distances.out','w')
    traj, top, selections = parser()
    sele1, sele2 = selections[0], selections[1]
    u = Universe(top, traj)
    sele1 = u.select_atoms(sele1)
    sele2 = u.select_atoms(sele2)

    if len(sele1) >=2: sys.exit("more than one atom")
    for frame in tqdm(u.trajectory):
        pos1, pos2 = sele1.positions[0], sele2.positions[0]
        distance   = np.linalg.norm(pos1 - pos2)
        flog.write(f'{distance}\n')

    flog.close()

if __name__ == '__main__':
    main()
