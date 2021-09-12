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
        self.__traj_name = traj
        self.__ref_name  = ref

    def rmsd(self, selection = 'name CA', rmsd_ref_name = ''):
        # This function uses MDAnalysis.analysis.rms.RMSD.
        # See https://docs.mdanalysis.org/stable/documentation_pages/analysis/rms.html
        
        # if rmsd_ref is not give, the reference that is used to create a "Universe" object is also utilised for the one in RMSD calc.
        if rmsd_ref_name == '': 
            rmsd_ref_name = self.__ref_name
            rmsd_ref      = self.__reference 
        else:
            rmsd_ref = Universe(rmsd_ref_name)

        print('### RMSD calculation ###')
        print(f'Selection           : {selection}')
        print(f'Reference structure : {rmsd_ref_name}')      
        print(f'    - {rmsd_ref}')      

        rmsd_obj = RMSD(atomgroup = self.__universe, 
                        reference = rmsd_ref,
                        select     = selection, 
                        groupselections = None,
                        weights = None, 
                        weights_groupselections = False,
                        tol_mass = 0.1, 
                        ref_frame = 0, 
                        verbose = True)
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

def parser():
    p = argparse.ArgumentParser()
    p.add_argument('-r', '--ref' , required = True)
    p.add_argument('-i', '--traj', required = True)
    p.add_argument('--rmsd_ref'  , required = False, default = None)
    p.add_argument('--rmsf_ref'  , required = False, default = None)
    p.add_argument('--sele_rmsd' , required = False, default = 'protein and name CA')
    p.add_argument('--sele_rmsf' , required = False, default = 'protein and name CA')
    p.add_argument('-sx', '--suffix', required = False, default = 'mda')
    args = p.parse_args()
    return args

def main():
    args = parser()
    ref  = args.ref 
    traj = args.traj 
    sele_rmsd = args.sele_rmsd
    sele_rmsf = args.sele_rmsf
    outsuffix = args.suffix
    rmsd_ref_file = args.rmsd_ref
    MDA = MDAnalyzer(ref, traj)

    #---RMSD calculation
    # (1) if a reference file for RMSD calc. is not given, the first frame of a trajectory is used as a reference.
    if rmsd_ref_file == None:
        rmsd = MDA.rmsd(sele_rmsd)   
    # (2) if a reference file is given, it is used as a reference structure.
    else: 
        rmsd = MDA.rmsd(sele_rmsd, rmsd_ref_file)   

    np.savetxt(f'rmsd_{outsuffix}.csv', rmsd, delimiter = ',', header = 'frame, time (ps), RMSD (A)')

    #---RMSF calculation
    rmsf = MDA.rmsf(sele_rmsf)
    pd.DataFrame(rmsf).to_csv(f'rmsf_{outsuffix}.csv', header = ['residue', 'RMSF'], index = None)

    exit()
    #---cPCA calculation
    cPC = MDA.cartesian_PCA() 

    #---dPCA calculation
    dPC = MDA.distance_based_PCA() 

    #---

if __name__ == '__main__':
    main()