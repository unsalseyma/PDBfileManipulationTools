#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 13:32:09 2022

@author: seyma unsal beyge
@email: seymaunsalbeyge [at] meddenovo.com

"""

import sys

f = open(sys.argv[1])
poseNum=int(sys.argv[2])
fout = sys.stdout

res_dic={}
line = f.readline()
while not line == '':
    if line.startswith('MODEL'):
        n = int(line.split()[1])
        if n == poseNum:
            while not 'ENDMDL' in line:
                if line.startswith('ATOM') \
                    or line.startswith('HETATM'): 
                    resname=line.split()[3]
                    resid=line.split()[4]
                    if resname not in res_dic.keys():
                        res_dic[resname]=(resid, [line])
                    else:
                        res_dic[resname][1].append(line)
                line = f.readline()
            break
    line = f.readline()
f.close()

res_dic = dict(sorted(res_dic.items(), key=lambda item: item[1]))
n=1
for k, v in res_dic.items():
    for line in v[1]:
        # atomid=line[7:11]
        line = line[:7]+"{:>4}".format(n)+line[11:66]+'\n'
        sys.stdout.write(line)
        n+=1


