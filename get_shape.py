#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 11:53:09 2024

@author: seyma
"""

## imports
# general import
import sys
import argparse   
import os
import random
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt 
import seaborn as sns 

from rdkit import Chem 
from rdkit.Chem import AllChem, Descriptors3D, rdmolfiles

def get_arguments():
    parser = argparse.ArgumentParser(description = 
                                     "a PMI analysis tool as part of focAlyze")
    
    # Arguments:
    parser.add_argument("-p", "--pdb", nargs='+', dest="pdbfiles",
                        type=str, action="store", required=True,  
                        help="[Required] path(s) to pdb file(s)")
    
    parser.add_argument("-o", "--out", dest="outDir", type=str, default="./results",
                help="Path to an output directory that the analysis report(s) will be saved")
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
        
    return parser 


class myProtein(object):
    def __init__(self, pdbfile):
        self.pdbfile = pdbfile 
        self.mol = self._pdb2mol() 
        self.gyr = None
        self.pmi1 = None
        self.pmi2 = None
        self.pmi3 = None
    
    def _pdb2mol(self):
        try:
            mol = Chem.MolFromPDBFile(self.pdbfile) 
            Chem.AddHs(mol)
            Chem.SanitizeMol(mol) 
        except:
            mol = None 
        return mol
    
    def find_descriptors(self):
        if self.mol is None:
            print("Pdb file cannot be converted into object.")
            return 
        # Find properties:
        gyr = Descriptors3D.RadiusOfGyration(self.mol,confId=-1) 
        pmi1 = Descriptors3D.PMI1(self.mol,confId=-1) 
        pmi2 = Descriptors3D.PMI2(self.mol,confId=-1) 
        pmi3 = Descriptors3D.PMI3(self.mol,confId=-1)
        self.gyr  = gyr
        self.pmi1 = pmi1
        self.pmi2 = pmi2
        self.pmi3 = pmi3 
        self.npr1 = pmi1 / pmi3 
        self.npr2 = pmi2 / pmi3
        return gyr, pmi1, pmi2, pmi3

    def calculate_likeness(self): 
        if None in [self.pmi1, self.pmi2, self.pmi3, self.npr1, self.npr2]:
            print("One of PMI values is not properly calculated. Please run 'find_descriptors' method. ")
            return
        
        p_rod, p_disk, p_sphere = None, None, None
        rod    =  np.array([0.0, 1.0]) 
        disk   =  np.array([0.5, 0.5]) 
        sphere =  np.array([1.0, 1.0]) 
        
        a = np.array([self.npr1, self.npr2]) 
        d_rod = 1/np.linalg.norm(a-rod) 
        d_disk = 1/np.linalg.norm(a-disk) 
        d_sphere = 1/np.linalg.norm(a-sphere)
        total = d_rod+d_disk+d_sphere
        p_rod = round(d_rod/total, 3) 
        p_disk = round(d_disk/total, 3)
        p_sphere = round(d_sphere/total, 3) 
        
        self.p_rod     = p_rod 
        self.p_disk    = p_disk 
        self.p_sphere  = p_sphere 
        
        return p_rod, p_disk, p_sphere

def plot_PMI(npr1_list, npr2_list, labels=None, save=""): 
    fig, ax = plt.subplots(figsize=(8, 8))
    pmi_plot = sns.scatterplot(x=npr1_list, y=npr2_list, ax=ax) 
    sns.lineplot(x=[0.0, 0.5, 1.0], y=[1.0, 0.5, 1.0], color="violet", ax=ax) 
    sns.lineplot(x=[0.0, 1.0], y=[1.0, 1.0], color="violet", ax=ax) 
    pmi_plot.set(xlabel="NPR1", ylabel="NPR2", title="Shape space analysis by PMI plot", xlim=(-0.15, 1.15), ylim=(0.35, 1.15)) 
    pmi_plot.text(0.53,0.49,"Disk") 
    pmi_plot.text(0.96,1.01,"Sphere") 
    pmi_plot.text(-0.04,1.01,"Rod") 
    pmi_plot.text(0.4,1.01,"Oblate-shape", color="purple") 
    pmi_plot.text(0.67,0.65,"Prolate-shape", rotation=60, color="purple")
    if not labels is None:
        for i in range(len(npr1_list)):
            x, y = npr1_list[i], npr2_list[i] 
            lbl = labels[i] 
            pmi_plot.text(x+0.001,y+0.001,lbl) 
            
    plt.savefig(save, dpi=400)
    return pmi_plot 

def main(): 
    parser = get_arguments()
    args = parser.parse_args()
    
    outDir = args.outDir 
    if not os.path.exists(outDir):
        os.makedirs(outDir) 
    else:
        print(f"\nOutput to {outDir}")
        sys.exit("Output directory already exists! Aborting..\n")

    pdbs = args.pdbfiles
    if len(pdbs) == 0:
        sys.exit("No pdb file provided! Aborting..\n") 

    start_time = datetime.now() 
    print(f"{'-'*50}\nStarted: {start_time}\n{'-'*50}\n") 
    
    outDF = pd.DataFrame(None, columns=["Label", 
                                        "Rod-likeness", 
                                        "Disk-likeness",
                                        "Sphere-likeness", 
                                        "Gyration"],
                         index = [os.path.basename(p) for p in pdbs])
    
    npr1List = [] 
    npr2List = []
    for pdb in pdbs:
        mymol= myProtein(pdb)
        gyr, pmi1, pmi2, pmi3 = mymol.find_descriptors()
        gyr = round(gyr, 3) 
        
        npr1List.append(mymol.npr1)
        npr2List.append(mymol.npr2)
        
        p_rod, p_disk, p_sphere = mymol.calculate_likeness()
        max_val = max(p_rod, p_disk, p_sphere) 
        if max_val == p_rod:
            label = f"Rod({int(p_rod*100)}%)" 
        elif max_val == p_disk:
            label = f"Disk({int(p_disk*100)}%)" 
        elif max_val == p_sphere:
            label = f"Sphere({int(p_sphere*100)}%)" 
        
        wp_rod, wp_disk, wp_sphere, wgyr = f"{p_rod:.3f}", f"{p_disk:.3f}", f"{p_sphere:.3f}", f"{gyr:.3f}" 
        
        outDF["Label"][os.path.basename(pdb)] = label
        outDF["Rod-likeness"][os.path.basename(pdb)] = wp_rod
        outDF["Disk-likeness"][os.path.basename(pdb)] = wp_disk 
        outDF["Sphere-likeness"][os.path.basename(pdb)] = wp_sphere 
        outDF["Gyration"][os.path.basename(pdb)] = wgyr
    
    
    # Standard PMI Analysis Plot
    plot = plot_PMI(npr1List, npr2List, 
                    labels=[os.path.basename(p) for p in pdbs], 
                    save=os.path.join(outDir, "shape_space_plot.png"))
    
    outDF.to_csv(os.path.join(outDir, "shape_results.csv"))
    print(outDF)
    print() 
    
    end_time = datetime.now()  
    runtime_total = round((end_time-start_time).total_seconds(), 3) 
    print(f"\n{'-'*50}\nEnded: {end_time}\n{'-'*50}") 
    print(f"Total runtime: {runtime_total} sec\n" ) 

if __name__ == "__main__":
    main()









