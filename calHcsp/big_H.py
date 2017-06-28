#!/usr/bin/env python
# -*- coding: utf-8 -*-
__description__ = \
"""
    big_H.py
    
    find the significantly perturbed protons
    example "big_H.py receptor.pdb pose_001.mol2 10 "
                       protein       ligand      distance cutoff  
    
    """
__author__ = "Zhuoqin Yu"


import numpy as np
import math
import sys, os
from operator import itemgetter
from moltools.get_csptype import CSPtype
receptor= sys.argv[1] # receptor pdb file
std_ligand=sys.argv[2] # ligand mol2 file used as a reference to filter residue within a certain distance
dis_cutoff=int(sys.argv[3]) # distance cutoff from ligand



#residrange=[150,30,36,166,167,168,169,42,171,172,173,174,175,176,177,178,179,39,40,170]# big residue and within 10A
#residrange=[36,166,167,168,169,42,171,174,175,178,39,40,170] #big residues and within 6 A


file1=open("exp_csp.txt") ##
lines1=file1.readlines()#lines=["line1","line2",...]
file1.close()
exp_csp_HA=[]
resid_list_HA=[]
atom_HA=[]
exp_csp_HN=[]
resid_list_HN=[]
atom_HN=[]
exp_csp_all=[]
exp_csp_HC=[]
resid_list_HC=[]
atom_HC=[]
exp_csp_PH=[]
resid_list_PH=[]
atom_PH=[]
for i in range(1,len(lines1)):
    linelist=lines1[i].split()
    csp_type=CSPtype(linelist[2],linelist[1])
    if csp_type=="HA":
        resid_list_HA.append(int(linelist[0]))
        atom_HA.append((int(linelist[0]),linelist[2]))
        exp_csp_HA.append(float(linelist[5])) #experimental complex cs - x2 calculated holo protein cs
    elif csp_type=="HN":
        resid_list_HN.append(int(linelist[0]))
        atom_HN.append((int(linelist[0]),linelist[2]))
        exp_csp_HN.append(float(linelist[5]))
    elif csp_type in ["HC","H1","H4","H5","HP","Har"]:
        resid_list_HC.append(int(linelist[0]))
        atom_HC.append((int(linelist[0]),linelist[2]))
        exp_csp_HC.append(float(linelist[5]))
    elif csp_type in ["HO","H"]:
        resid_list_PH.append(int(linelist[0]))
        atom_PH.append((int(linelist[0]),linelist[2]))
        exp_csp_PH.append(float(linelist[5]))


mean_HA=np.mean(exp_csp_HA)
std_HA=np.std(exp_csp_HA)
mean_HN=np.mean(exp_csp_HN)
std_HN=np.std(exp_csp_HN)
#print("***",mean_HN,std_HN)
mean_HC=np.mean(exp_csp_HC)
std_HC=np.std(exp_csp_HC)
mean_PH=np.mean(exp_csp_PH)
std_PH=np.std(exp_csp_PH)

######## 3sigma cycle for HN only:
#threesig=3*std_HN
#old_csp_HN=list(exp_csp_HN)
#new_csp_HN=[]
#for i in old_csp_HN:
#    if i < mean_HN+threesig and i > mean_HN-threesig:
#        new_csp_HN.append(i)
#while len(old_csp_HN)>len(new_csp_HN):
#    old_csp_HN=list(new_csp_HN)
#    new_csp_HN=[]
#    threesig=3*np.std(old_csp_HN)
#    mean_HN=np.mean(old_csp_HN)
#    for i in old_csp_HN:
#        if i < mean_HN+threesig and i > mean_HN-threesig:
#            new_csp_HN.append(i)
#std_HN=np.std(old_csp_HN)
#
#print("after cycle",mean_HN,std_HN)
##########

big_HA_resid=[]
big_HN_resid=[]
big_HC_resid=[]
big_PH_resid=[]
big_atoms=[]
for i in range(0,len(exp_csp_HA)):
    if exp_csp_HA[i]> mean_HA+std_HA or exp_csp_HA[i]<mean_HA-std_HA:
        big_HA_resid.append(resid_list_HA[i])
        big_atoms.append(atom_HA[i])
for i in range(0,len(exp_csp_HN)):
    if exp_csp_HN[i]> mean_HN+1*std_HN or exp_csp_HN[i]<mean_HN-1*std_HN:
    #if exp_csp_HN[i]>0.05 or exp_csp_HN<-0.05:
        big_HN_resid.append(resid_list_HN[i])
        big_atoms.append(atom_HN[i])
for i in range(0,len(exp_csp_HC)):
    if exp_csp_HC[i]> mean_HC+1*std_HC or exp_csp_HC[i]<mean_HC-1*std_HC:
        big_HC_resid.append(resid_list_HC[i])
        big_atoms.append(atom_HC[i])
for i in range(0,len(exp_csp_PH)):
    if exp_csp_PH[i]> mean_PH+1*std_PH or exp_csp_PH[i]<mean_PH-1*std_PH:
        big_PH_resid.append(resid_list_PH[i])
        big_atoms.append(atom_PH[i])

#big_resid=[val for val in big_HA_resid if val in big_HN_resid]
#print("Resids with significant CSP (exp_csp, HN): ")
#print(','.join(map(str, big_HN_resid)))
#print("Resids with significant CSP (exp_csp, HA): ")
#print(','.join(map(str, big_HA_resid)))
#print("Resids with significant CSP (exp_csp, HA or HN): ")
#print(','.join(map(str, set(big_HA_resid+big_HN_resid))))
#big_csp_res=list(set(big_HA_resid+big_HN_resid))
#print('big atoms')
#print(big_atoms)
reslist=[]
atomlist=[]
for atom in big_atoms:
    reslist.append(atom[0])
    atomlist.append(atom)
print("All resids with significant perturbed protons")
print(','.join(map(str,set(reslist))))
print("All significantly perturbed protons")
print(atomlist)

##read in protein file
profile=open(receptor)
lines=profile.readlines()
profile.close()
resids=[]
xs=[]
ys=[]
zs=[]
for i in range(0,len(lines)):
    linelist=lines[i].split()
    if linelist[0]=="ATOM" and linelist[3] not in ["ACE","NMA"]:
        resids.append(int(linelist[5]))
        xs.append(float(linelist[6]))
        ys.append(float(linelist[7]))
        zs.append(float(linelist[8]))
pro_matrx=np.array([xs,ys,zs])
##


currentresids=[]
## calculate the residues that are within distance cutoff of ligand
ligfile=open(std_ligand)
lines=ligfile.readlines()
ligfile.close()
for line in lines:
    if "@<TRIPOS>ATOM" in line: startline=lines.index(line)+1
    if "@<TRIPOS>BOND" in line: endline=lines.index(line)
lxs=[]
lys=[]
lzs=[]
resid_10a_dict=[]
for n in range(startline,endline):
    linelist=lines[n].split()
    lxs.append(float(linelist[2]))
    lys.append(float(linelist[3]))
    lzs.append(float(linelist[4]))
lig_mat=np.array([lxs,lys,lzs])
for pt in range(0,len(lxs)):
    lpt=np.array([[lxs[pt]],[lys[pt]],[lzs[pt]]])
    cal1=(pro_matrx-lpt)**2
    dis_mat_sq=cal1[0]+cal1[1]+cal1[2]
    for id in range(0,len(dis_mat_sq)):
        if dis_mat_sq[id]<=dis_cutoff**2:
            currentresids.append(resids[id])
resid_10a_dict=list(set(currentresids))
#residstr=','.join(resid_10a_dict[ligf])
print("Residues within "+str(dis_cutoff)+" A from the ligand:")
print(','.join(map(str,resid_10a_dict)))
print("Residues have significantly perturbed protons & within "+str(dis_cutoff)+" A from the ligand:")
finallist=[val for val in set(reslist) if val in resid_10a_dict]
print(','.join(map(str,finallist)))

print("Residues have significantly perturbed HA & within "+str(dis_cutoff)+" A from the ligand:")
finallist=[val for val in big_HA_resid if val in resid_10a_dict]
print(','.join(map(str,finallist)))

print("Residues have significantly perturbed HN & within "+str(dis_cutoff)+" A from the ligand:")
finallist=[val for val in big_HN_resid if val in resid_10a_dict]
print(','.join(map(str,finallist)))

print("Residues have significantly perturbed nonpolar sidechain protons & within "+str(dis_cutoff)+" A from the ligand:")
finallist=[val for val in big_HC_resid if val in resid_10a_dict]
print(','.join(map(str,finallist)))

print("Residues have significantly perturbed polar sidechain protons & within "+str(dis_cutoff)+" A from the ligand:")
finallist=[val for val in big_PH_resid if val in resid_10a_dict]
print(','.join(map(str,finallist)))

