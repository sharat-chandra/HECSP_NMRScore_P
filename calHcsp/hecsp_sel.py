#!/usr/bin/env python
# Filename: Hecsp.py
# -*- coding: utf-8 -*-

"""
    This is the Hecsp.py program written by Zhuoqin Yu in Merz Research Group
    at Michigan State University.
"""
#==============================================================================
# Load the module
#==============================================================================
from moltools.moltool import *
from moltools.classdata_mol2 import *
from moltools.classdata_csp import *
from math import *
from coordtrans.eigsys import *
import numpy as np
from coordtrans.coords import *
from copy import deepcopy
from ringpath.dfs import *
import warnings
import os
from optparse import OptionParser

parser = OptionParser("usage: -p protein_input_file -l ligand_input_file -s/--parm parm set \n"
                      "       [--cut    cutoff distance between residue and ligand] \n"
                      "       [-htp     proton types included in scoring] \n"
                      "       [-resid   residue ids included in scoring] \n"
                      "       [-v       verbose; print out CSP from each effect separately")
parser.add_option("-p", dest="pro_input", type='string',
                  help="Input file name for protein")
parser.add_option("-l", dest="lig_input", type='string',
                  help="Input file name for ligand")                  
parser.add_option("-s", "--parm", dest="set", type='string',
                  help="Parameter set selection: \n    gls: Gasteiger & least squares;\n    gloocv: Gasteiger & leave-one-out;\n    als: AM1-BCC & least squares")
parser.add_option("-v", "--verbose", dest="verbose", type='string',
                  help="verbose; print out CSP from each effect separately",default="False")
parser.add_option("--cut", dest="cutoff", type='float',
                  help="cutoff distance, default=0",default=0.0) 
parser.add_option("--htp", dest="proton_types", type='string',
                  help="proton types included, example: HA HN",default="HA HN H H1 H4 H5 HC HP Har HO")  
parser.add_option("--resid", dest="residue_IDs", type='string',
                  help="residue ids included, example: 20,21,22",default="ALL")
(options, args) = parser.parse_args()


#==============================================================================
# Useful constants and lists
#==============================================================================
Navo=6.022e23    
unsatbondtype=["2","3","am"]
scH_amber=["H","H1","H4","H5","HC","HP","Har"]
backbone_atomname=["H1","H2","H3","HA","HA2","HA3","HA1","H","N","CA","CB","C","HN","HN1","HN2","HN3"]
#csptypelist=["HA","HN","N","CA","CB","CO","H","H1","H4","H5","HC","HP","HO"]
efedisdict={"HA":10.0,"HN":10.0,"N":3.5,"CA":4.0,"CB":3.5,"CO":10.0,"H":10.0,\
"H1":10.0,"H1N":10.0,"H4":10.0,"H5":10.0,"HC":10.0,"HP":10.0,"Har":10.0,"HO":10.0,"HS":10.0}
hbdisdict={"HN":3.5,"HA":3.5,"H":3.5,"HO":3.5}
rcedisdict={"HA":6.0,"HN":6.0,"H":6.0,"H1":6.0,"H4":6.0,"H5":6.0,"HC":6.0,"HP":6.0,"Har":6.0,"HO":6.0,"HS":6.0}
#totringtypes=['CNNNN','CCNCN','CCCNN','CCCCS','CCNCO','CCCCCC','CCCCCN','CCCNCN','CNCNN'] #add ring type moved to the parameter file
Cref=191.775425
Href=31.75509167

# Print the title of the program
#version = '1.0'
#print_title('hecsp.py', version)
#==============================================================================
# parameters determined
#==============================================================================
if options.set=="gls":
    from parms.glsparm import *
elif options.set=="gloocv":
    from parms.gloocvparm import *
elif options.set=="als":
    from parms.alsparm import *

ring_I_dict={}
for i in range(len(totringtypes)):
    ring_I_dict[totringtypes[i]]=B[i]

def func(p,x): #calculate for HA HN # add ring type
    F,B0,B1,B2,B3,B4,B6,B7,B8,epi,a,b,t=p
    a0,a1,a2,a3,a4,a5,a6,a7,a8,e8,c9,c10,AS=x #B5 is for Benzene, keep it as 1.0
    return F*B0*a0+F*B1*a1+F*B2*a2+F*B3*a3 +F*B4*a4 +F*1*a5 +F*B6*a6+F*B7*a7+F*B8*a8+epi*e8+ c9*(a+b*c10)+t*AS
def sh_func(p,x): #calculate for side H
    F,B0,B1,B2,B3,B4,B6,B7,B8,epi,t=p
    a0,a1,a2,a3,a4,a5,a6,a7,a8,e8,AS=x #B5 is for Benzene, keep it as 1.0
    return F*B0*a0+F*B1*a1+F*B2*a2+F*B3*a3 +F*B4*a4 +F*1*a5 +F*B6*a6+F*B7*a7+F*B8*a8+epi*e8+t*AS
def func_efe(p,x):
    return p*x
def func_rce(p,x):
    F,B0,B1,B2,B3,B4,B6,B7,B8=p
    a0,a1,a2,a3,a4,a5,a6,a7,a8=x #B5 is for Benzene, keep it as 1.0
    return F*B0*a0+F*B1*a1+F*B2*a2+F*B3*a3 +F*B4*a4 +F*1*a5 +F*B6*a6+F*B7*a7+F*B8*a8
def func_ani(p,x):
    return p*x
def func_hb(p,x):
    a,b=p
    c9,c10=x #B5 is for Benzene, keep it as 1.0
    return c9*(a+b*c10)

#==============================================================================
# Functions to determine the effect that the nuclei is under
#==============================================================================
def efe_Nuc(atom): #whether consider the efe based on protein atomname
    if atom.predict:
        return True
    else:
        return False

def RCE_Nuc(atom): #judge whether under ringcurrent effet
    if atom.predict:
        return True
    else: return False
    
def HB_Nuc(atom): #judge whether can be under HB effet
    if atom.predict and atom.csptype in ["HA","HN"]:
        return True
    elif atom.predict and atom.csptype in ["H","HO"]:
        return True
    else: return False

def sybyl_hb_donor(sybyl):
    if sybyl in ["N.3","N.2","O.3","S.3","N.am","N.pl3","N.4"]:
        return True
    else:
        return False

def sybyl_hb_acceptor(sybyl):
    if sybyl in ["N.3","N.2","O.3","S.3","N.1","N.ar","O.2","O.co2","F","Br","Cl"]:
        return True
    else:
        return False

def Anyefe(cspterm): #whether the protein atom is under effective efe
    if cspterm!="NAN": 
        if cspterm.efe.efedepsi==0.0:
            return False
        else:
            return True
    else:
        return False
def Anyrce(cspterm):
    add=0.0
    if cspterm!="NAN":   
        for i in list(cspterm.rce.values()):
            add=add+i.G
        if add==0.0:
            return False
        else:
            return True
    else:
        return False
        
def Anyhb(cspterm):
    if cspterm!="NAN":        
        if cspterm.hb.lengthHA!=0.0:
            return True
        else:
            return False
    else:
        return False


def cal_A(rce): # X is the horizontal distance to the cm ,Y is the vertical distance
    Alist=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0] # add ring type
    for key in rce.keys():
        if rce[key].G!=0.0:
            for n in range(1,len(rce[key].ring)+1):
                newring=RotateMe(rce[key].ring,1,n)
                flipring=newring[::-1]
                if newring in totringtypes:
                    Alist[totringtypes.index(newring)]=Alist[totringtypes.index(newring)]+(1-(rce[key].costheta)**2*3)/(rce[key].distance**3)
                    break
                elif flipring in  totringtypes:
                    Alist[totringtypes.index(flipring)]=Alist[totringtypes.index(flipring)]+(1-(rce[key].costheta)**2*3)/(rce[key].distance**3)
                    break
    return Alist
#==============================================================================
# Read the local atom susceptibility from file and stored in sus_localatom_dict
#==============================================================================
hecsphome=os.getenv('HECSPHOME')
sus_localatom_dict={}
infile=open(hecsphome+"/lib/atom_susceptibility")
lines=infile.readlines()
infile.close()
for i in range(1,len(lines)):
    localatom,xkai,ykai,zkai=lines[i].split() #first line in the file are title
    sus_localatom_dict[localatom]=Susceptibility_atom(localatom,xkai,ykai,zkai)   
#==============================================================================
# Read in protein sybyl atom type from  atom_nm_amber_sybyl file and stored in pro_sybyl_dict
#==============================================================================
pro_sybyl_dict={}
infile=open(hecsphome+"/lib/atom_nm_amber_sybyl")
lines=infile.readlines()
infile.close()
okayresidue=[]
for i in range(1,len(lines)):
    number,residue,atomname,ambertype,sybyltype=lines[i].split()  
    okayresidue.append(residue)
    pro_sybyl_dict[(residue+atomname)]=Pro_sybyl(residue,atomname,ambertype,sybyltype)
    

def CSPtype(atomname,residue):
    if residue in okayresidue:
        if atomname in backbone_atomname:# CB included
            if atomname in ["H1","H2","H3","H","HN","HN1","HN2","HN3"]:
                return "HN"
            elif atomname in ["HA","HA2","HA3","HA1"]:
                return "HA"
            elif atomname=="C":
                return "CO"
            else:
                return atomname
        elif pro_sybyl_dict[(residue+atomname)].ambertype=="HA": # sidechain aromatic H is "HA" amber atomtype. we reassign as "Har"
            return "Har"
        else:
            return pro_sybyl_dict[(residue+atomname)].ambertype
    else:
        return "NAN" 
#==============================================================================
# Enter loops for different pdb structures
#==============================================================================
if options.pro_input[-4::]==".pdb" or options.pro_input[-6::]==".pdbqt":
    aa=options.pro_input.split(".")[0]
else:
    aa=options.pro_input
if "&" not in options.lig_input:
    if options.lig_input[-5::]==".mol2":
        bb=options.lig_input.split(".")[0]
    else:
        bb=options.pro_input
    pl_complex=[(aa,bb)]
#pl_complex=[("1e66","1e66")]#(p,l)

for code in pl_complex:
    
    print(code)
#==============================================================================
#    # Read ligand information from mol2 file and store in dictionary.
#==============================================================================    
    Latomdict={}
    Lbonddict={}

    if options.lig_input[-5::]==".mol2":
        fmol2=open(options.lig_input) ##
        lines=fmol2.readlines()#lines=["line1","line2",...]
        fmol2.close()
    else:
        fmol2=open(code[1]+".mol2")
        lines=fmol2.readlines()#lines=["line1","line2",...]
        fmol2.close()
# section the ligand mol2 file
    mol2secline=[]
    for i in range(0,len(lines)):
        if lines[i][0]=="@":
            mol2secline.append([lines[i][9:-1],i])
    mol2secline.append(["END",len(lines)-1])
    for i in range(0,len(mol2secline)):
        if mol2secline[i][0]=="ATOM":
            atomrange=range(mol2secline[i][1]+1,mol2secline[i+1][1])
        elif mol2secline[i][0]=="BOND":
            bondrange=range(mol2secline[i][1]+1,mol2secline[i+1][1])
# Put the atom and bond information in Latomdict and Lbonddict. 
    Lxyzlist=[]
    for i in atomrange:
        atnumber,atname,x,y,z,attype,resid,resname,charge=lines[i].split()
        element=get_element(attype)
        mass=get_mass(element)
        Latomdict[atnumber]=Atom(atnumber,atname,x,y,z,attype,resid,resname,charge,element,mass)
        Lxyzlist.append((float(x),float(y),float(z)))
    for i in bondrange:
        bnumber,atnumber1,atnumber2,btype=lines[i].split()
        Lbonddict[bnumber]=Bond(bnumber,atnumber1,atnumber2,btype)
#==============================================================================
#     calculate the center of the ligand
#==============================================================================
#    xsum=0
#    ysum=0
#    zsum=0
#    for i in Lxyzlist:
#        xsum=xsum+i[0]
#        ysum=ysum+i[1]
#        zsum=zsum+i[2]
#    Lcenter=(xsum/len(Lxyzlist),ysum/len(Lxyzlist),zsum/len(Lxyzlist))
        #==============================================================================
        #  Generate dictionary which contains the bounded atoms and bondtypes of each atom       
        #==============================================================================
    abdict=get_atombonds(Latomdict,Lbonddict)
        #==============================================================================
        #     Find polar atoms in ligand
        #==============================================================================
    polaratoms=[] #["1","3","5",...] items are atom numbers
    #print(abdict)
    for i in Latomdict:
        if Latomdict[i].element not in ["C","H"]:
            polaratoms.append(i)
        elif Latomdict[i].element in ["C","H"]:
            for j in abdict[i][0]: # j is [[bondatomnumbers],[bondtypes]]
                if Latomdict[j].element not in ["C","H"]:
                    polaratoms.append(i)
                    break
        #    print("polar atoms")
        #    print(polaratoms)
        #==============================================================================
        #     Find the aromatic rings in ligand (base line algorithm)
        #==============================================================================
        
    graph={}
    for i in Lbonddict:
        if Lbonddict[i].bondtype=="ar":
            if Lbonddict[i].atomnumber1 not in graph:
                graph[Lbonddict[i].atomnumber1]=[Lbonddict[i].atomnumber2]
                if Lbonddict[i].atomnumber2 not in graph:
                    graph[Lbonddict[i].atomnumber2]=[Lbonddict[i].atomnumber1]
                elif Lbonddict[i].atomnumber2 in graph:
                    graph[Lbonddict[i].atomnumber2].append(Lbonddict[i].atomnumber1)
            elif Lbonddict[i].atomnumber1 in graph:
                graph[Lbonddict[i].atomnumber1].append(Lbonddict[i].atomnumber2)
                if Lbonddict[i].atomnumber2 not in graph:
                    graph[Lbonddict[i].atomnumber2]=[Lbonddict[i].atomnumber1]
                elif Lbonddict[i].atomnumber2 in graph:
                    graph[Lbonddict[i].atomnumber2].append(Lbonddict[i].atomnumber1)
    
    cycles=[]
    cyclessort=[]
    for start in graph:
        ends=[]
        if len(graph[start])==2:
            ends.append(graph[start][0])
        elif len(graph[start])==1:
            ends.append(graph[start][0])
#        elif len(graph[start]>=3:
#            for i in graph[start]
#            ends.append
        for end in ends:
            path=find_all_paths(graph,start,end)
            for i in path:
                if len(i)>=4 and len(i)<=7:
                    if sorted(i) not in cyclessort:
                        cycles.append(i)
                        cyclessort.append(sorted(i))

    aroringdict={}
    if len(cycles)!=0:
        for c in cycles:
            for i in range(0,len(c)):
                c[i]=str(c[i])
            c=tuple(c)
            intensity=ring_intensity(c,Latomdict,ring_I_dict)
#            rcescale=rce_scale(c,Latomdict)
            R12=get_bondvector((Latomdict[c[0]].x,Latomdict[c[0]].y,Latomdict[c[0]].z),(Latomdict[c[1]].x,Latomdict[c[1]].y,Latomdict[c[1]].z))
            R1f=get_bondvector((Latomdict[c[0]].x,Latomdict[c[0]].y,Latomdict[c[0]].z),(Latomdict[c[-1]].x,Latomdict[c[-1]].y,Latomdict[c[-1]].z))           
            norm=np.cross(R12,R1f)/np.sqrt(vec_magnitude(R12)**2*vec_magnitude(R1f)**2-(np.dot(R12,R1f))**2) #np.array which is a unit vector
            ring_cm=ring_CM(c,Latomdict)
            aroringdict[c]=Aro_ring(c,norm,intensity,ring_cm)
#    print("aromatic ring",cycles)
#    print(cycles)
    
    cycleselement=deepcopy(cycles)
    for c in cycleselement:
        for i in range(len(c)):
            c[i]=Latomdict[c[i]].element        
    
        
                        
                    
        #==============================================================================
        #     group susceptibility effects
        #==============================================================================
    # 1erg/(G^2mole)=1cm^3/mol
    localatom_xyz_axis_dict={} #{'1':Local_xyz_axis,....}
    anisoatom=[]
    anisogroup=[]
    for bond in Lbonddict:
        if Lbonddict[bond].bondtype in unsatbondtype:
            anisoatom.append(Lbonddict[bond].atomnumber1)
            anisoatom.append(Lbonddict[bond].atomnumber2)
    anisoatom=list(set(anisoatom)) # remove duplicate in list anisoatom
#    print(anisoatom)

    for atnum in anisoatom[:]:# put the localatom infomation in anisoatom and generate localatom_xyz_axis_dict
        localatom=get_localatom(Latomdict[atnum],abdict[atnum][1])
        if localatom=="NONE":
            anisoatom.remove(atnum)
#            anisoatom[anisoatom.index(atnum)]=[atnum,localatom]        
        else:
            anisoatom[anisoatom.index(atnum)]=[atnum,localatom]# put the localatom infomation in anisoatom
            mass=Latomdict[atnum].mass
            atom1=(Latomdict[atnum].x,Latomdict[atnum].y,Latomdict[atnum].z)
            if localatom=="11C2":
                bdto=abdict[atnum][0][abdict[atnum][1].index("2")] #the atom number of the atom which is double bonded to the local atom
                atom2=(Latomdict[bdto].x,Latomdict[bdto].y,Latomdict[bdto].z)
                xaxis=get_bondvector(atom1,atom2)
                bdto=abdict[atnum][0][abdict[atnum][1].index("2")-1]                
#                bdto=abdict[atnum][0][abdict[atnum][1].index("1")]
                atom2=(Latomdict[bdto].x,Latomdict[bdto].y,Latomdict[bdto].z)
                singlebond=get_bondvector(atom1,atom2)
                zaxis=np.cross(xaxis,singlebond)
                yaxis=np.cross(zaxis,xaxis)
                localatom_xyz_axis_dict[atnum]=Local_xyz_axis(atnum,localatom,xaxis,yaxis,zaxis,mass)
            elif localatom=="O2":
                bdto=abdict[atnum][0][0]
                atom2=(Latomdict[bdto].x,Latomdict[bdto].y,Latomdict[bdto].z)
                xaxis=get_bondvector(atom1,atom2)
                bdtoto=abdict[bdto][0][abdict[bdto][1].index("2")-1]                
#                bdtoto=abdict[bdto][0][abdict[bdto][1].index("1")]
                atom3=(Latomdict[bdtoto].x,Latomdict[bdtoto].y,Latomdict[bdtoto].z)
                bdtovec=get_bondvector(atom2,atom3)
                zaxis=np.cross(xaxis,bdtovec)
                yaxis=np.cross(zaxis,xaxis)
                localatom_xyz_axis_dict[atnum]=Local_xyz_axis(atnum,localatom,xaxis,yaxis,zaxis,mass)
            elif localatom=="1C3":
                bdto=abdict[atnum][0][abdict[atnum][1].index("3")]
                atom2=(Latomdict[bdto].x,Latomdict[bdto].y,Latomdict[bdto].z)
                xaxis=get_bondvector(atom1,atom2)
                yaxis=Rz90(np.mat(xaxis).transpose())
                zaxis=np.cross(xaxis,yaxis)
                localatom_xyz_axis_dict[atnum]=Local_xyz_axis(atnum,localatom,xaxis,yaxis,zaxis,mass)
            elif localatom=="N3":
                bdto=abdict[atnum][0][0]
                atom2=(Latomdict[bdto].x,Latomdict[bdto].y,Latomdict[bdto].z)
                xaxis=get_bondvector(atom1,atom2)
                yaxis=Rz90(np.mat(xaxis).transpose())
                zaxis=np.cross(xaxis,yaxis)
                localatom_xyz_axis_dict[atnum]=Local_xyz_axis(atnum,localatom,xaxis,yaxis,zaxis,mass)
            elif localatom=="11N1":
                if Latomdict[atnum].atomtype=="N.am":
                    bdto=abdict[atnum][0][abdict[atnum][1].index("am")]
                    atom2=(Latomdict[bdto].x,Latomdict[bdto].y,Latomdict[bdto].z)
                    xaxis=get_bondvector(atom2,atom1) #xaxis is from atom2 to N.am
                    bdto=abdict[atnum][0][abdict[atnum][1].index("1")]
                    atom3=(Latomdict[bdto].x,Latomdict[bdto].y,Latomdict[bdto].z)
                    bdtovec=get_bondvector(atom1,atom3)
                    zaxis=np.cross(xaxis,bdtovec)
                    yaxis=np.cross(zaxis,xaxis)
                else:
                    bdto=abdict[0][0]
                    atom2=(Latomdict[bdto].x,Latomdict[bdto].y,Latomdict[bdto].z)
                    xaxis=get_bondvector(atom2,atom1)
                    bdto=abdict[0][1]
                    atom3=(Latomdict[bdto].x,Latomdict[bdto].y,Latomdict[bdto].z)
                    bdtovec=get_bondvector(atom1,atom3)
                    zaxis=np.cross(xaxis,bdtovec)
                    yaxis=np.cross(zaxis,xaxis)
                localatom_xyz_axis_dict[atnum]=Local_xyz_axis(atnum,localatom,xaxis,yaxis,zaxis,mass)                                    
            elif localatom=="S2":
                bdto=abdict[atnum][0][0]
                atom2=(Latomdict[bdto].x,Latomdict[bdto].y,Latomdict[bdto].z)
                xaxis=get_bondvector(atom1,atom2)
                yaxis=Rz90(np.mat(xaxis).transpose())
                zaxis=np.cross(xaxis,yaxis)
                localatom_xyz_axis_dict[atnum]=Local_xyz_axis(atnum,localatom,xaxis,yaxis,zaxis,mass) 
            elif localatom=="P3":
                bdto=abdict[atnum][0][0]
                atom2=(Latomdict[bdto].x,Latomdict[bdto].y,Latomdict[bdto].z)
                xaxis=get_bondvector(atom1,atom2)
                yaxis=Rz90(np.mat(xaxis).transpose())
                zaxis=np.cross(xaxis,yaxis)
                localatom_xyz_axis_dict[atnum]=Local_xyz_axis(atnum,localatom,xaxis,yaxis,zaxis,mass) 
    anisopair=[]
    anisogroup_pre=[[]]
    for i in anisoatom:
        for j in anisoatom[anisoatom.index(i)+1:]:
            if i[0] in abdict[j[0]][0]:
                anisopair.append([i,j])
                judge=0
                addcommon=[]
                for grp in anisogroup_pre:
                    if i in grp or j in grp:
                        grp.append(i)
                        grp.append(j)
                        addcommon.append(grp)
                        judge=judge+1
                if judge==0:
                    anisogroup_pre.append([i,j])
                if judge>=2:
                    anisogroup_pre.append(sum(addcommon,[]))
                    for g in addcommon:
                        anisogroup_pre.remove(g)
    anisogroup_pre.remove([])
    anisogroup=[]
    for grp in anisogroup_pre:
        anisogroup.append(make_unique(grp))
#    print("anisogroup")
#    print(anisogroup) #[[['4', '11C2'], ['3', '1N11'], ['2', '11C2'], ['1', 'O2']]]
   
    group_abc_axis_dict={}    
    for grp in anisogroup: # get the principal axis a b c and stored in group_abc_axis_dict
        n=len(grp)
        m=[0]
        x=[0]
        y=[0]
        z=[0]
        for atm in grp: #get m,x,y,z lists
            m.append(Latomdict[atm[0]].mass)
            x.append(Latomdict[atm[0]].x)
            y.append(Latomdict[atm[0]].y)            
            z.append(Latomdict[atm[0]].z)        
        x_old=list(x)
        y_old=list(y)
        z_old=list(z)
        cm_old=[0.0,0.0,0.0] # the coordinates of CM will be store in cm_old 
        MovetoCM(m,x,y,z,n,cm_old)
        pa_cab=PA_cab(m, x, y, z, n, -1)
        group_abc_axis_dict[anisogroup.index(grp)]=Group_abc_axis(grp,pa_cab[1],pa_cab[2],pa_cab[0],cm_old)
    
    sus_group_dict={}    
    for grp in anisogroup: # get the Kaiaa,Kaibb,Kaicc of the groups
        aakai,bbkai,cckai=0,0,0
        for atm in grp:
            aakai=aakai+cosine(localatom_xyz_axis_dict[atm[0]].xaxis,group_abc_axis_dict[anisogroup.index(grp)].aaxis)**2*sus_localatom_dict[atm[1]].xkai+\
            cosine(localatom_xyz_axis_dict[atm[0]].yaxis,group_abc_axis_dict[anisogroup.index(grp)].aaxis)**2*sus_localatom_dict[atm[1]].ykai+\
            cosine(localatom_xyz_axis_dict[atm[0]].zaxis,group_abc_axis_dict[anisogroup.index(grp)].aaxis)**2*sus_localatom_dict[atm[1]].zkai
            bbkai=bbkai+cosine(localatom_xyz_axis_dict[atm[0]].xaxis,group_abc_axis_dict[anisogroup.index(grp)].baxis)**2*sus_localatom_dict[atm[1]].xkai+\
            cosine(localatom_xyz_axis_dict[atm[0]].yaxis,group_abc_axis_dict[anisogroup.index(grp)].baxis)**2*sus_localatom_dict[atm[1]].ykai+\
            cosine(localatom_xyz_axis_dict[atm[0]].zaxis,group_abc_axis_dict[anisogroup.index(grp)].baxis)**2*sus_localatom_dict[atm[1]].zkai
            cckai=cckai+cosine(localatom_xyz_axis_dict[atm[0]].xaxis,group_abc_axis_dict[anisogroup.index(grp)].caxis)**2*sus_localatom_dict[atm[1]].xkai+\
            cosine(localatom_xyz_axis_dict[atm[0]].yaxis,group_abc_axis_dict[anisogroup.index(grp)].caxis)**2*sus_localatom_dict[atm[1]].ykai+\
            cosine(localatom_xyz_axis_dict[atm[0]].zaxis,group_abc_axis_dict[anisogroup.index(grp)].caxis)**2*sus_localatom_dict[atm[1]].zkai
        sus_group_dict[anisogroup.index(grp)]=Susceptibility_group(grp,aakai,bbkai,cckai)
##        
        #==============================================================================
        #   Read in protein information from pdb file
        #==============================================================================
    Patomdict={}
    csptermsdict={}
#    bb_dihedral={}
#    sc_dihedral={}
    resid_list=[]
    if options.pro_input[-4::]==".pdb" or options.pro_input[-6::]==".pdbqt":
        file=open(options.pro_input) ##
        lines=file.readlines()#lines=["line1","line2",...]
        file.close()
    else:
        file=open(code[0]+".pdb") ##
        lines=file.readlines()#lines=["line1","line2",...]  
        file.close()
    for i in range(0,len(lines)):
        if lines[i][0:4]=="ATOM" or lines[i][0:6]=="HETATM":
#            currentlinelist=lines[i].split()
#            posa,atomnum,atomname,residue,chain,resid,x,y,z,posb,posc,posd=lines[i].split()[0:12] # 6-31G** TMS C(ref)=191.775425 H(ref)=31.75509167
#            posa,atomnum,atomname,residue,resid,x,y,z,posb,posc,posd=lines[i].split()[0:11] # 6-31G** TMS C(ref)=191.775425 H(ref)=31.75509167
            atomnum=lines[i][6:11].replace(" ", "")
            atomname=lines[i][12:16].replace(" ","")
            residue=lines[i][17:20].replace(" ","")
            chain=lines[i][21]
            resid=lines[i][22:26].replace(" ","")
            x=lines[i][30:38].replace(" ","")
            y=lines[i][38:46].replace(" ","")
            z=lines[i][46:54].replace(" ","")
            csptype=CSPtype(atomname,residue)
            if residue in okayresidue:
                amber=pro_sybyl_dict[(residue+atomname)].ambertype
                sybyl=pro_sybyl_dict[(residue+atomname)].sybyltype
                csp=0.0
                if int(resid) not in resid_list:
                    resid_list.append(int(resid))
                if options.cutoff == 0 and options.residue_IDs=="ALL": # no cutoff nor resid id range
                    Patomdict[i]=Res_atom(atomname,atomnum,residue,resid,x,y,z,csp,amber,sybyl,csptype,True)
                elif options.cutoff !=0 or options.residue_IDs!="ALL":
                    Patomdict[i]=Res_atom(atomname,atomnum,residue,resid,x,y,z,csp,amber,sybyl,csptype,False)
                if Patomdict[i].predict:#give each predictable nuclear a start values
                    anishift=0.0
                    efe=Elecfieldeffect(0.0,epsdict[csptype])
                    hb=HBeffect(Patomdict[i],0.0,0.0,0.0,0.0)
                    rce={}
                    csptermsdict[i]=CSP_terms(efe,hb,rce,anishift)
                    patom=(Patomdict[i].x,Patomdict[i].y,Patomdict[i].z)
                    
    if options.cutoff !=0:
        resids_in_cutoff=[]
        for i in Patomdict:
            if Patomdict[i].resid not in resids_in_cutoff:
                for dot in Lxyzlist:
                    if get_distance((Patomdict[i].x,Patomdict[i].y,Patomdict[i].z),dot) <= options.cutoff:
                        resids_in_cutoff.append(Patomdict[i].resid)
                        break
                    else: continue
            else: continue
    else:
        resids_in_cutoff=resid_list
        
    if options.residue_IDs!="ALL":
        resids_in_range=map(int, options.residue_IDs.split(","))
    else:
        resids_in_range=resid_list
    if options.cutoff !=0 or options.residue_IDs!="ALL": 
        resid_met_both=list(set(resids_in_cutoff) & set(resids_in_range))
        for i in Patomdict:
            if Patomdict[i].resid in resid_met_both:
                Patomdict[i].metcutoff=True
            else:
                Patomdict[i].metcutoff=False
            #=============================================================================
            #           # calculate magnetic anisotropic shift of all predictable nuclei
            #==============================================================================
                        
    for i in csptermsdict:
        if Patomdict[i].metcutoff:                
            for grp in group_abc_axis_dict:         
                gcm=group_abc_axis_dict[grp].cm
                vec_p_gcm=get_bondvector(patom,gcm)
                r=vec_magnitude(vec_p_gcm)*1e-8 #convert from anstrong to centimeter
                csptermsdict[i].anishift+=(1/3)*(1/Navo)*pow(r,-3)*\
                (sus_group_dict[grp].akai*(3*cosine(vec_p_gcm,group_abc_axis_dict[grp].aaxis)**2-1)+\
                sus_group_dict[grp].bkai*(3*cosine(vec_p_gcm,group_abc_axis_dict[grp].baxis)**2-1)+\
                sus_group_dict[grp].ckai*(3*cosine(vec_p_gcm,group_abc_axis_dict[grp].caxis)**2-1))   
    #==============================================================================
    #   #   calculate the efe for each predictable nucli           
    #==============================================================================                
    for i in csptermsdict:
        if Patomdict[i].metcutoff:
            if efe_Nuc(Patomdict[i]): # the nuclei that is under efe effect.
                patom=(Patomdict[i].x,Patomdict[i].y,Patomdict[i].z)
                for plr in polaratoms:
                    plratom=(Latomdict[plr].x,Latomdict[plr].y,Latomdict[plr].z)
                    d=get_distance(patom,plratom)
                    if d<=efedisdict[Patomdict[i].csptype]:
                        q=Latomdict[plr].atomcharge
                        idx=partner_key_dict(Patomdict,i)
                        partatom=(Patomdict[idx].x,Patomdict[idx].y,Patomdict[idx].z)
                        ptoplr=np.array(plratom)-np.array(patom) 
                        ptopart=np.array(partatom)-np.array(patom)
                        costheta=cosine(ptoplr,ptopart)
                        csptermsdict[i].efe.efedepsi+=q*costheta/(d**2)
    #==============================================================================
    # # calculate aromatic ring geometrical factor G
    #==============================================================================
    if len(cycles)!=0:
        for c in cycles:
            c=tuple(c)
            normal=aroringdict[c].norm
            cm=aroringdict[c].cm
            ringelement="".join(cycleselement[cycles.index(list(c))])
            for i in csptermsdict: #calculate G for target i referring to ring c
                if Patomdict[i].metcutoff:
                    distancetocm=get_distance((Patomdict[i].x,Patomdict[i].y,Patomdict[i].z),cm)
                    costheta=np.dot(normal,get_bondvector((Patomdict[i].x,Patomdict[i].y,Patomdict[i].z),cm))/(vec_magnitude(normal)*distancetocm)
                    if RCE_Nuc(Patomdict[i]) and distancetocm<=rcedisdict[Patomdict[i].csptype]:
                        targetxyz=(Patomdict[i].x,Patomdict[i].y,Patomdict[i].z)                        
                        pro_O=project_on_plane((Patomdict[i].x,Patomdict[i].y,Patomdict[i].z),normal,(Latomdict[c[0]].x,Latomdict[c[0]].y,Latomdict[c[0]].z))
                        #array
                        geof=0.0
                        for k in range(0,len(c)): 
                            ROk=get_bondvector(pro_O,(Latomdict[c[k]].x,Latomdict[c[k]].y,Latomdict[c[k]].z))
                            Rtark=get_bondvector(targetxyz,(Latomdict[c[k]].x,Latomdict[c[k]].y,Latomdict[c[k]].z))
                            if k==len(c)-1:
                                ROk1=get_bondvector(pro_O,(Latomdict[c[0]].x,Latomdict[c[0]].y,Latomdict[c[0]].z))
                                Rtark1=get_bondvector(targetxyz,(Latomdict[c[0]].x,Latomdict[c[0]].y,Latomdict[c[0]].z))
                            else:   
                                ROk1=get_bondvector(pro_O,(Latomdict[c[k+1]].x,Latomdict[c[k+1]].y,Latomdict[c[k+1]].z))
                                Rtark1=get_bondvector(targetxyz,(Latomdict[c[k+1]].x,Latomdict[c[k+1]].y,Latomdict[c[k+1]].z))                                
                            dk=1/(vec_magnitude(Rtark)**3)+1/(vec_magnitude(Rtark1)**3)
                            if np.dot(np.cross(ROk,ROk1),normal)>0.0: #parallel
                                sign=-1
                            else:
                                sign=1
                            areaOkk1=sign*0.5*vec_magnitude(np.cross(ROk,ROk1))
                            geof += dk*areaOkk1
                        csptermsdict[i].rce[c]=Ringcurrenteffect(ringelement,geof,Targetfactordict[Patomdict[i].csptype],costheta,distancetocm,float(distancetocm)*np.sqrt(1-float(costheta)**2),abs(float(distancetocm)*float(costheta)))
                    else:
                        csptermsdict[i].rce[c]=Ringcurrenteffect(ringelement,0.0,Targetfactordict[Patomdict[i].csptype],0.0,0.0,0.0,0.0)
#                    
    #                    
    #==============================================================================
    # Find all Hydrogen bonds                
    #==============================================================================
    HBdict={}
##### HN and HA HB effects, HA determine by itself
    for i in csptermsdict:
        if Patomdict[i].metcutoff:
            if HB_Nuc(Patomdict[i]) and get_element(Patomdict[i].sybyl)=="H": #target nuclear is donor H 
                targetnuc=Patomdict[i]
                for j in Latomdict:
                    if sybyl_hb_acceptor(Latomdict[j].atomtype):
                        lengthHA=get_distance((Patomdict[i].x,Patomdict[i].y,Patomdict[i].z),(Latomdict[j].x,Latomdict[j].y,Latomdict[j].z))                
                        if lengthHA<hbdisdict[Patomdict[i].csptype]: #H...A<3.5
                            pairatom=Latomdict[j]
                            donor_idx=partner_key_dict(Patomdict,i)
                            lengthDA=get_distance((Patomdict[donor_idx].x,Patomdict[donor_idx].y,Patomdict[donor_idx].z),(Latomdict[j].x,Latomdict[j].y,Latomdict[j].z))
                            lengthDH=get_distance((Patomdict[donor_idx].x,Patomdict[donor_idx].y,Patomdict[donor_idx].z),(Patomdict[i].x,Patomdict[i].y,Patomdict[i].z))                      
                            cosangDHA=(lengthHA**2+lengthDH**2-lengthDA**2)/(2*lengthHA*lengthDH)
                            if lengthDA<=3.5 and cosangDHA<0.0: #D...A<3.6 and angle DHA larger than 90 degrees.
    #                            print("HB",i,Patomdict[i].csptype,lengthHA)
                                if i not in HBdict:
                                    HBdict[i]=HBeffect(targetnuc,pairatom,lengthHA,lengthDA,cosangDHA)                             
                                elif i in HBdict:
                                    if HBdict[i].lengthHA >lengthHA: #the current length is shorter than previous one, we need to replace it.
                                        HBdict[i]=HBeffect(targetnuc,pairatom,lengthHA,lengthDA,cosangDHA)
                    
    if len(HBdict)!=0:
        for i in HBdict:
            csptermsdict[i].hb=HBdict[i]
#==============================================================================
# bounded leastsqbound fit all parameters    collect parameters
##==============================================================================        
    cal_nuc=options.proton_types.split(" ")
#    cal_nuc=["HA","HN","H","H1","H4","H5","HC","HP","Har","HO"] 
    lines=[]
    with open(code[0]+"_"+code[1]+"_csp","w") as cspfile:
        width=(15,8,6,15,25,2)
        line=["Conditions:","cutoff:",str(options.cutoff),"  proton_types:",options.proton_types,"\n"]
        cspfile.write("".join("%*s" % i for i in zip(width, line)))        
        width=(20,30,2)
        line=["     Selected residue IDs: ",options.residue_IDs,"\n"]
        cspfile.write("".join("%*s" % i for i in zip(width, line)))
        if options.verbose=="False":
            width=(10,10,10,10,15,2)
            line=["ATOM","ATOM_ID","RES_ID","RESIDUE","CSP(ppm)","\n"]
        elif options.verbose=="True":
            width=(10,10,10,10,15,15,15,15,15,2)
            line=["ATOM","ATOM_ID","RES_ID","RESIDUE","CSP(ppm)","RCE","EFE","HB","ANISO","\n"]
        cspfile.write("".join("%*s" % i for i in zip(width, line)))
#        sumofcsp2=0
        for i in csptermsdict:
            if Patomdict[i].csptype in cal_nuc and Patomdict[i].metcutoff:
                if Anyefe(csptermsdict[i]):
                    E=[csptermsdict[i].efe.efedepsi]
                else:E=[0.0]
                if Anyhb(csptermsdict[i]):
                    C=(1/(csptermsdict[i].hb.lengthHA**3),csptermsdict[i].hb.lengthHA**3)                  
                else: C=(0.0,0.0)
                if Anyrce(csptermsdict[i]):  #totringtypes=['CNNNN','CCNCN','CCCNN','CCCCS','CCNCO','CCCCCC','CCCCCN','CCCNCN']
                    A=cal_A(csptermsdict[i].rce)
                else: A=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]# add ring type
                newy=0
                if Patomdict[i].csptype=="HA":             
                    newy=func((Targetfactordict["HA"],B[0],B[1],B[2],B[3],B[4],B[6],B[7],B[8],epsdict["HA"],HB_dict["HA"][0],HB_dict["HA"][1],AnishiftCoeffdict["HA"]),\
                    (A[0],A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8],E[0],C[0],C[1],csptermsdict[i].anishift))
                    effects=[func_rce((Targetfactordict["HA"],B[0],B[1],B[2],B[3],B[4],B[6],B[7],B[8]),(A[0],A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8])),func_efe(epsdict["HA"],E[0]),func_hb((HB_dict["HA"][0],HB_dict["HA"][1]),(C[0],C[1])),func_ani(AnishiftCoeffdict["HA"],csptermsdict[i].anishift)]
                elif Patomdict[i].csptype=="HN":             
                    newy=func((Targetfactordict["HN"],B[0],B[1],B[2],B[3],B[4],B[6],B[7],B[8],epsdict["HN"],HB_dict["HN"][0],HB_dict["HN"][1],AnishiftCoeffdict["HN"]),\
                    (A[0],A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8],E[0],C[0],C[1],csptermsdict[i].anishift))
                    effects=[func_rce((Targetfactordict["HN"],B[0],B[1],B[2],B[3],B[4],B[6],B[7],B[8]),(A[0],A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8])),func_efe(epsdict["HN"],E[0]),func_hb((HB_dict["HN"][0],HB_dict["HN"][1]),(C[0],C[1])),func_ani(AnishiftCoeffdict["HN"],csptermsdict[i].anishift)]
                elif Patomdict[i].csptype=="H":             
                    newy=func((Targetfactordict["H"],B[0],B[1],B[2],B[3],B[4],B[6],B[7],B[8],epsdict["H"],HB_dict["H"][0],HB_dict["H"][1],AnishiftCoeffdict["H"]),\
                    (A[0],A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8],E[0],C[0],C[1],csptermsdict[i].anishift))
                    effects=[func_rce((Targetfactordict["H"],B[0],B[1],B[2],B[3],B[4],B[6],B[7],B[8]),(A[0],A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8])),func_efe(epsdict["H"],E[0]),func_hb((HB_dict["H"][0],HB_dict["H"][1]),(C[0],C[1])),func_ani(AnishiftCoeffdict["H"],csptermsdict[i].anishift)]
                elif Patomdict[i].csptype=="HO":             
                    newy=func((Targetfactordict["HO"],B[0],B[1],B[2],B[3],B[4],B[6],B[7],B[8],epsdict["HO"],HB_dict["HO"][0],HB_dict["HO"][1],AnishiftCoeffdict["HO"]),\
                    (A[0],A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8],E[0],C[0],C[1],csptermsdict[i].anishift))
                    effects=[func_rce((Targetfactordict["HO"],B[0],B[1],B[2],B[3],B[4],B[6],B[7],B[8]),(A[0],A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8])),func_efe(epsdict["HO"],E[0]),func_hb((HB_dict["HO"][0],HB_dict["HO"][1]),(C[0],C[1])),func_ani(AnishiftCoeffdict["HO"],csptermsdict[i].anishift)]
                elif Patomdict[i].csptype in ["Har","H1","H4","H5","HC","HP"]:
                    newy=sh_func((Targetfactordict[Patomdict[i].csptype],B[0],B[1],B[2],B[3],B[4],B[6],B[7],B[8],epsdict[Patomdict[i].csptype],AnishiftCoeffdict[Patomdict[i].csptype]),\
                    (A[0],A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8],E[0],csptermsdict[i].anishift))
                    effects=[func_rce((Targetfactordict[Patomdict[i].csptype],B[0],B[1],B[2],B[3],B[4],B[6],B[7],B[8]),(A[0],A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8])),func_efe(epsdict[Patomdict[i].csptype],E[0]),0.000,func_ani(AnishiftCoeffdict[Patomdict[i].csptype],csptermsdict[i].anishift)]
                if options.verbose=="False":
                    line=[Patomdict[i].atomname,int(Patomdict[i].atomnum),Patomdict[i].resid,Patomdict[i].residue,("%.3f" % round(newy,3)),"\n"]
                elif options.verbose=="True":
                    line=[Patomdict[i].atomname,int(Patomdict[i].atomnum),Patomdict[i].resid,Patomdict[i].residue,("%.3f" % round(newy,3)),("%.3f" % round(effects[0],3)),("%.3f" % round(effects[1],3)),("%.3f" % round(effects[2],3)),("%.3f" % round(effects[3],3)),"\n"]
                lines.append(line)
#               sumofcsp2=sumofcsp2+newy**2
        lines.sort(key=lambda x:x[1])
        for line in lines:
            line[1]=str(line[1])
            cspfile.write("".join("%*s" % i for i in zip(width, line)))
#print("sum of csp^2:",str(sumofcsp2))
quit()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
