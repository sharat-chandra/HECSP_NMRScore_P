# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 10:50:35 2015

@author: zhuoqinyu
"""
import numpy as np
import math

def find(string, ch):
    for i, ltr in enumerate(string):
        if ltr == ch:
            yield i
 #============================================================================            
def contained(candidate, container):
    temp = container[:]
    try:
        for v in candidate:
            temp.remove(v)
        return True
    except ValueError:
        return False
 #============================================================================
def get_element(sybyltype):#get the element from sybyl atomtype
    a=sybyltype.split('.')[0]
    return a.title() #only the first letter is capitalized
 #============================================================================
def get_mass(element):
    table={ 'H'  :   1.0079 , 'He' :   4.0026 , 'Li' :   6.941  ,
         'Be' :   9.0122 , 'B'  :  10.811  , 'C'  :  12.0107 ,
         'N'  :  14.0067 , 'O'  :  15.9994 , 'F'  :  18.9984 ,
         'Ne' :  20.1797 , 'Na' :  22.9898 , 'Mg' :  24.3050 ,
         'Al' :  26.9815 , 'Si' :  28.0855 , 'P'  :  30.9738 , 
         'S'  :  32.065  , 'Cl' :  35.453  , 'Ar' :  39.948  ,
         'K'  :  39.0983 , 'Ca' :  40.078  , 'Sc' :  44.9559 , 
         'Ti' :  47.867  , 'V'  :  50.9415 , 'Cr' :  51.9961 ,
         'Mn' :  54.9380 , 'Fe' :  55.845  , 'Co' :  58.9331 , 
         'Ni' :  58.6934 , 'Cu' :  63.546  , 'Zn' :  65.409  ,
         'Ga' :  69.723  , 'Ge' :  72.64   , 'As' :  74.9216 , 
         'Se' :  78.96   , 'Br' :  79.904  , 'Kr' :  83.798  ,
         'Rb' :  85.4678 , 'Sr' :  87.62   , 'Y'  :  88.9059 , 
         'Zr' :  91.224  , 'Nb' :  92.9064 , 'Mo' :  95.94   ,
         'Tc' :  98.     , 'Ru' : 101.07   , 'Rh' : 102.9055 , 
         'Pd' : 106.42   , 'Ag' : 107.8682 , 'Cd' : 112.411  ,
         'In' : 114.818  , 'Sn' : 118.710  , 'Sb' : 121.760  , 
         'Te' : 127.60   , 'I'  : 126.9045 , 'Xe' : 131.293  ,
         'Cs' : 132.9055 , 'Ba' : 137.327  , 'La' : 138.9055 , 
         'Ce' : 140.116  , 'Pr' : 140.9077 , 'Nd' : 144.242  ,
         'Pm' : 145.     , 'Sm' : 150.36   , 'Eu' : 151.964  , 
         'Gd' : 157.25   , 'Tb' : 158.9254 , 'Dy' : 162.500  ,
         'Ho' : 164.9303 , 'Er' : 167.259  , 'Tm' : 168.9342 , 
         'Yb' : 173.04   , 'Lu' : 174.967  , 'Hf' : 178.49   ,
         'Ta' : 180.9479 , 'W'  : 183.84   , 'Re' : 186.207  , 
         'Os' : 190.23   , 'Ir' : 192.217  , 'Pt' : 195.084  ,
         'Au' : 196.9666 , 'Hg' : 200.59   , 'Tl' : 204.3833 , 
         'Pb' : 207.2    , 'Bi' : 208.9804 , 'Po' : 209.     ,
         'At' : 210.     , 'Rn' : 222.     , 'Fr' : 223.     , 
         'Ra' : 226.     , 'Ac' : 227.     , 'Th' : 232.0381 ,
         'Pa' : 231.0359 , 'U'  : 238.0289 , 'Np' : 237.     , 
         'Pu' : 244.     , 'Am' : 243.     , 'Cm' : 247.     ,
         'Bk' : 247.     , 'Cf' : 251.     , 'Es' : 252.     , 
         'Fm' : 257.     , 'Md' : 258.     , 'No' : 259.     ,
         'Lr' : 262.     , 'Rf' : 261.     , 'Db' : 262.     , 
         'Sg' : 266.     , 'Bh' : 264.     , 'Hs' : 277.     ,
         'Mt' : 268.     , 'Ds' : 281.     , 'Rg' : 272.     , 
         'Cn' : 285.     , 'Uut': 284.     , 'Uuq': 289.     ,
         'Uup': 288.     , 'Uuh': 292.     , 'Uus': 291.     , 
         'Uuo': 294.     , 'EP' : 0.000000 }     
    return table[element]                                  
 #============================================================================    
def get_atombonds(atoms,bonds):#atoms and bonds are Atom and Bond class 
    Atom_bonds={}
    for i in atoms:
        blist=[[],[]]
        for j in bonds:
            if bonds[j].atomnumber1==atoms[i].atomnumber:
                blist[0].append(bonds[j].atomnumber2)
                blist[1].append(bonds[j].bondtype)
            elif bonds[j].atomnumber2==atoms[i].atomnumber:
                blist[0].append(bonds[j].atomnumber1)
                blist[1].append(bonds[j].bondtype)
#        blist=list(set(blist))
        Atom_bonds[i]=blist
    return Atom_bonds
 #============================================================================
def get_localatom(atom,bondlist):
    if atom.element=='C' and bondlist.count("2")==1 and len(bondlist)==3:
        return "11C2"
    elif atom.element=="C" and bondlist.count("3")==1 and len(bondlist)==2:
        return "1C3"
    elif atom.element=="O" and bondlist.count("2")==1 and len(bondlist)==1:
        return "O2"
    elif atom.element=="N" and bondlist.count("3")==1 and len(bondlist)==1:
        return "N3"
    elif atom.atomtype in ["N.am","N.pl3"] and len(bondlist)==3:
        return "11N1"
    elif atom.element=="S" and len(bondlist)==1:
        return "S2"
    elif atom.element=="P" and len(bondlist)==1:
        return "P3"
    elif atom.element=="F" and len(bondlist)==1:
        return "F"
    else:
        return "NONE"
#============================================================================
def get_distance(atom1,atom2):
    dis=0
    for i in range(0,len(atom1)):
        dis=dis+(atom1[i]-atom2[i])**2
    return np.sqrt(dis)
#============================================================================
def vec_magnitude(vec):
    add=0
    for i in vec:
        add=add+i**2
    return np.sqrt(add)
#============================================================================      
def get_bondvector(atom1,atom2):
    vec=np.array(atom1) #return np.array
    for i in range(0,len(atom1)):
        vec[i]=atom2[i]-atom1[i]
    return vec #array type
#============================================================================        
def Rz90(vec): # vec is colum vector
    R=np.matrix([[0,-1,0],[1,0,0],[0,0,1]])
    return np.array(R*vec).T[0] # result is array

#============================================================================    
def cosine(x,y): #x and y are two vectors(list)
    return np.dot(x,y)/(np.sqrt(np.dot(x,x))*np.sqrt(np.dot(y,y)))
#============================================================================    
def common_elements(list1, list2):
    result = []
    for element in list1:
        if element in list2:
            result.append(element)
    return result    
#============================================================================
def make_unique(original_list):
    unique_list = []
    for obj in original_list:
        if obj not in unique_list:
            unique_list.append(obj)
    return unique_list
#============================================================================
def partner_key_dict(Patomdict,atm):
    for i in Patomdict:
        if Patomdict[i].atomname==Patomdict[atm].partner[0]:
            if Patomdict[i].resid==Patomdict[atm].resid+int(Patomdict[atm].partner[1]):
                return i
                break
            
def reverse_partner_key_dict(Patomdict,atm):
    for i in Patomdict:
        if Patomdict[atm].atomname==Patomdict[i].partner[0]:
            if Patomdict[atm].resid==Patomdict[i].resid+int(Patomdict[i].partner[1]):
                return i 
                break
def idx_patomdict(Patomdict,reid,atname):
    for i in Patomdict:
        if Patomdict[i].atomname==atname and Patomdict[i].resid==reid:
            return i
#============================================================================
def ring_CM(cycle,Latomdict):
    n=len(cycle)
    mCM = 0e0; xCM = 0e0; yCM = 0e0; zCM = 0e0
    for i in range(0,n):
        mCM += Latomdict[cycle[i]].mass
        xCM += Latomdict[cycle[i]].mass*Latomdict[cycle[i]].x
        yCM += Latomdict[cycle[i]].mass*Latomdict[cycle[i]].y
        zCM += Latomdict[cycle[i]].mass*Latomdict[cycle[i]].z
    xCM /= mCM; yCM /= mCM; zCM /= mCM
    return (xCM,yCM,zCM)
#============================================================================
def project_on_plane(point,norm,ptonplane): #norm is unit vector
    vec1=np.array(point)-np.array(ptonplane)
    pro_norm=np.dot(vec1,norm)*norm # pro_norm is the vec1 project onto norm vector
    pro_point=np.array(point)
    for i in range(0,len(point)):
        pro_point[i]=point[i]-pro_norm[i]
    return pro_point #np.array   
    
    
def dihedral_ang(p1,p2,p3,p4):
    b1=np.array(p2)-np.array(p1)
    b2=np.array(p3)-np.array(p2)
    b3=np.array(p4)-np.array(p3)
    n1=np.cross(b1,b2)/np.linalg.norm(np.cross(b1,b2))
    n2=np.cross(b2,b3)/np.linalg.norm(np.cross(b2,b3))
    m1=np.cross(n1,b2/np.linalg.norm(b2))
    x=np.dot(n1,n2)
    y=np.dot(m1,n2)
    return math.atan2(y,x)*180/math.pi #from -180 to 180
#===========================================================================
def locateHAformHN(patomdict,i):
    for j in patomdict:
        if patomdict[j].csptype=="HA" and patomdict[j].resid== patomdict[i].resid-1:
            return j
        if j==len(patomdict)-1:
            return "NAN"
    