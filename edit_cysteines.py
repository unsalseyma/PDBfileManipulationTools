#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 10:52:54 2022

@author: seyma
"""
import numpy as np
import sys

usage = "USAGE:\tpython edit_cysteines.py <path_to_the_pdbFile> <path_to_the_output_pdbFile> <varname_of_protein_in_leap>"
try:
    fin = open(sys.argv[1], 'r')
except:
    sys.exit(usage)
try:
    fout =open(sys.argv[2], 'w')
except:
    sys.exit(usage)
try:
    prot =str(sys.argv[3])
except:
    sys.exit(usage)

alllines = []
cut_s = 2.5
all_cys = []
res_prev = -1
nres = 0
nstart = 0
line = fin.readline()
while not line == '':
    if line.startswith('ATOM') or line.startswith('HETATM'):
        atname  = line[12:16].strip()
        resname = line[17:20].strip()
        ires = int(line[22:26])
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        if nstart == 0:
            nstart = ires
        if ires != res_prev: 
            res_prev = ires
            nres += 1
        if resname in ['CYS','CYM','CYX'] and atname == 'SG':
            all_cys.append([ires,nres,np.array([x,y,z])])
    line = line.replace(line[22:26],"{:>4}".format(nres))
    alllines.append(line)
    line = fin.readline()
fin.close()

bridges = []
bridged_res = []
cys_ires = []
ncys = len(all_cys)
for i in range(ncys):
    ii, ni, ci = all_cys[i]
    for j in range(i+1,ncys):
        ij, nj, cj = all_cys[j]
        d = np.linalg.norm(cj-ci)
        if d <= cut_s:
            bridges.append([ii,ij,ni,nj])
            if ii not in bridged_res:
                bridged_res.append(ni)
            if ij not in bridged_res:
                bridged_res.append(nj)

bridged_res = sorted(bridged_res)
newlines=[]
if nstart != 1:
    remark = f"REMARK\tShift of residue numbers\t({nstart} to 1)\n"
    newlines.append(remark)
for line in alllines:
    if line.startswith('ATOM') or line.startswith('HETATM'):
        atname  = line[12:16].strip()
        resname = line[17:20].strip()
        ires = int(line[22:26].strip())
        if ires in bridged_res:
            line=line.replace('CYS', 'CYX')
    newlines.append(line)

for line in newlines:
    fout.write(line)
fout.close()

sbonds = [f"bond {prot}.{x[2]}.SG {prot}.{x[3]}.SG\n" for x in bridges]
for line in sbonds:
    sys.stdout.write(line)
    