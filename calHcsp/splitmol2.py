#!/usr/bin/env python


__description__ = \
"""
splitmol2.py

split the mol2 file containing multiple poses from glide docking
example "splitmol2.py poses.mol2 pose" will generate pose_001.mol2 pose_002.mol2 ... ...

"""
__author__ = "Zhuoqin Yu"


import sys, os
import itertools

def isplit(iterable,splitters):
    return [list(g) for k,g in itertools.groupby(iterable,lambda x:x in splitters) if not k]

def splitNMR(pdb):
    """
    Split each model in an NMR pdb file into its own pdb.
    """
    linerange=[]
    to_strip = "ROOT"
    for i in range(len(pdb)):
        if to_strip not in pdb[i]:
            linerange.append(i)
        elif to_strip in pdb[i]:
            linerange.append(i)
            linerange.append('NAN')
    newlinerange=isplit(linerange,('NAN',))
    return newlinerange

def main():
    """
    Function to call if run from commmand line.
    """

    pdb_file=sys.argv[1]
    newhead=sys.argv[2]

    f = open(pdb_file,'r')
    pdb = f.readlines()
    f.close()

    linerange= splitNMR(pdb)
        
		#short_pdb = os.path.split(pdb_file)[-1][:-5]
    for j in range(len(linerange)):
        g = open("%s_%s.mol2" % (newhead,str(int(j)+1).zfill(3)),"w")
        for line in linerange[j]:
            g.writelines(pdb[line])
        g.close()

if __name__ == "__main__":
    main()

