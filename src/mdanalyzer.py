#!/usr/bin/env python3
import MDAnalysis as mda
from MDAnalysis import Universe
from MDAnalysis.analysis.rms import RMSD
from MDAnalysis.analysis.rms import RMSF
from MDAnalysis.analysis import align
import numpy as np
import pandas as pd
import argparse

class MDAnalyzer():
    
    def __init__(self, ref, traj):
        self.__universe  = Universe(ref,traj)
        self.__reference = Universe(ref)

    def rmsd(self, selection = 'name CA'):
        print(f'Selection in RMSD calc.: {selection}')
       
        rmsd_obj = RMSD(atomgroup  = self.__universe, 
                        referenece = self.__reference, 
                        select     = selection, 
                        groupselections = None,
                        weights = None, 
                        weights_groupselections = False,
                        tol_mass = 0.1,
                        ref_frame = 0)
        rmsd_obj.run(start=None, stop=None, step=None, verbose=True)
        return rmsd_obj.rmsd
                
    def rmsf(self, selection = 'protein and name CA'):
        print(f'Selection in RMSF calc.: {selection}')
        aligner = align.AlignTraj(mobile    = self.__universe, 
                                  reference = self.__reference, 
                                  select=selection, 
                                  in_memory=True)
        aligner.run(start=None, stop=None, step=None, verbose=True)
        target = self.__universe.select_atoms(selection)
        rmsf_obj = RMSF(target, verbose = True).run()
        residue_names = [ resn+str(resi) for resn, resi in zip(target.resnames, target.resids) ]
        return list(zip(residue_names, rmsf_obj.rmsf))

    def cartesian_PCA(self, selection = 'name CA'):
        print(f'Selection in cPCA calc.: {selection}')
        
        return 0

    def distance_based_PCA(self, selection = 'name CA'):
        print(f'Selection in dPCA calc.: {selection}')
        
        return 0

    def sasa(self, type = 'residue'):
        pass

    def make_hist(self):
        pass

def main():
    p = argparse.ArgumentParser()
    p.add_argument('-r', '--ref' , required = True)
    p.add_argument('-i', '--traj', required = True)
    p.add_argument('-sx', '--suffix', required = False, default = 'mda')
    args = p.parse_args()
    ref = args.ref 
    traj = args.traj 
    outsuffix = args.suffix

    MDA = MDAnalyzer(ref, traj)

    #---RMSD calculation 
    rmsd = MDA.rmsd('(resid 103:109 and name CA) or (resid 142:149 and name CA)')
    np.savetxt(f'rmsd_{outsuffix}.csv', rmsd, delimiter = ',', header = 'frame, time (ps), RMSD (A)')

    #---RMSF calculation
    rmsf = MDA.rmsf()
    pd.DataFrame(rmsf).to_csv(f'rmsf_{outsuffix}.csv', header = ['residue', 'RMSF'], index = None)

    #---cPCA calculation
    cPC = MDA.cartesian_PCA() 

    #---dPCA calculation
    dPC = MDA.distance_based_PCA() 

    #---

if __name__ == '__main__':
    main()