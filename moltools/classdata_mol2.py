# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 10:58:31 2015

@author: zhuoqinyu
"""
class Atom:
    def __init__(self,atnumber,atname,x,y,z,attype,resid,resname,charge,element,mass):
        self.atomnumber=atnumber
        self.atomname=atname
        self.x=float(x)
        self.y=float(y)
        self.z=float(z)
        self.atomtype=attype
        self.resid=resid
        self.resname=resname
        self.atomcharge=float(charge)
        self.element=element
        self.mass=float(mass)
        
class Bond:
    def __init__(self,bnumber,atnumber1,atnumber2,btype):
        self.bondnumber=bnumber
        self.atomnumber1=atnumber1
        self.atomnumber2=atnumber2
        self.bondtype=btype
        
class Molecule:
    def __init__(self,atom,bond):
        self.atom=atom
        self.bond=bond
        
class Aro_ring:
    def __init__(self,ring,norm,intensity,cm):
        self.ring=ring
        self.norm=norm
        self.intensity=intensity
        self.cm=cm
#        self.rcescale=rcescale


#==============================================================================
# Determine ring type and ring intensity
#==============================================================================
totringtypes=['CNNNN','CCNCN','CCCNN','CCCCS','CCNCO','CCCCCC','CCCCCN','CCCNCN']  

def RotateMe(text_ori,mode=0,steps=1):
  # Takes a text string and rotates
  # the characters by the number of steps.
  # mode=0 rotate right
  # mode=1 rotate left
  text=text_ori
  length=len(text)
 
  for step in range(steps):
  # repeat for required steps
 
    if mode==0:
      # rotate right
      text=text[length-1] + text[0:length-1]
    else:
      # rotate left
      text=text[1:length] + text[0]
 
  return text  
  

def ringI(ring,ring_I_dict):
    for n in range(1,len(ring)+1):
        newring=RotateMe(ring,1,n)
        flipring=newring[::-1]
        if newring in ring_I_dict:
            return ring_I_dict[newring]
            break
        elif flipring in ring_I_dict:
            return ring_I_dict[flipring]
            break
    else:
        return 0.0
    
def ring_intensity(cyc,Latomdict,ring_I_dict):
    ring='' #string

    for i in cyc:
        ring=ring+Latomdict[i].element # generate the ring string
    for n in range(1,len(ring)+1):
        newring=RotateMe(ring,1,n)
        flipring=newring[::-1]
        if newring in ring_I_dict:
            return ring_I_dict[newring]
            break
        elif flipring in ring_I_dict:
            return ring_I_dict[flipring]
            break
    else:
        return 0.0
        

    
#==============================================================================
# objects used in susceptibility calculation
#==============================================================================
class Susceptibility_atom:
    def __init__(self,localatom,xkai,ykai,zkai):
        self.localatom=localatom
        self.xkai=float(xkai)
        self.ykai=float(ykai)
        self.zkai=float(zkai)
        
class Local_xyz_axis:
    def __init__(self,atomnumber,localatom,xaxis,yaxis,zaxis,mass):
        self.atomnumber=atomnumber
        self.localatom=localatom
        self.xaxis=xaxis
        self.yaxis=yaxis
        self.zaxis=zaxis
        self.mass=mass
        
class Group_abc_axis:
    def __init__(self,group,aaxis,baxis,caxis,cm):
        self.group=group
        self.aaxis=aaxis
        self.baxis=baxis
        self.caxis=caxis
        self.cm=cm
        
class Susceptibility_group:
    def __init__(self,group,akai,bkai,ckai):
        self.group=group
        self.akai=akai
        self.bkai=bkai
        self.ckai=ckai

        