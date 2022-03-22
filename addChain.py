#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 19:15:56 2022

@author: seyma
"""

# py script to convert pdb file to a new version with chain information

import sys

f = sys.argv[1]
try:
    changeRes=sys.argv[2]  # should be restart or False
    # if residue numbers will be updated, It should be True. Otherwise, False.
except:
    changeRes=False
    # sys.exit("\n**** WARNING****\nYou should enter second input which is either 'True' or 'False'.\n"
    #          "True:\tIf you want to restart residue number for each chain.\n"
    #          "False:\tIf you don't want to change residue number.\n")
fout = sys.stdout

flines = []
data = open(f)
line = data.readline()
while not line == '':
    flines.append(line)
    line = data.readline()
data.close()

import string
alphabet = list(string.ascii_uppercase)

newlines = []
def addChainonly(lines):
    i = 0
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            newline = line[:21]+alphabet[i]+line[22:]
            newlines.append(newline)
        elif line.startswith('TER'):
            if len(line)>22:
                newline = line[:21]+alphabet[i]+line[22:]
                newlines.append(newline)
            else:
                newlines.append(line)
            i+=1
        else:
            newlines.append(line)
    return newlines

def addChainAndRestartRes(lines):
    i = 0
    current = 0
    prev = -1
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            resno = int(line[22:26].strip())
            if resno != prev:
                prev = resno
                current += 1 #the new residue number, changed based on the change of residue name
                newresno = "{:>4}".format(current)
            else:
                newresno = "{:>4}".format(current)
            newline = line[:21]+alphabet[i]+newresno+line[26:]
            newlines.append(newline)
        elif line.startswith('TER'):
            if len(line)>3:
                resno = int(line[22:26].strip())
                if resno != prev:
                    prev = resno
                    current += 1 #the new residue number, changed based on the change of residue name
                    newresno = "{:>4}".format(current)
                else:
                    newresno = "{:>4}".format(current)
                newline = line[:21]+alphabet[i]+newresno+line[26:]
                newlines.append(newline)
            else:
                newlines.append(line)
            i+=1
            current = 0
        else:
            newlines.append(line)
    return newlines


if changeRes == 'restart':
    newlines = addChainAndRestartRes(flines)
    for line in newlines:
        fout.write(f"{line}")
elif not changeRes:
    newlines = addChainonly(flines)
    for line in newlines:
        fout.write(f"{line}")
else:
    f.close()
    sys.exit("Second input should be either 'restart' or None.\n")

