# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 13:50:52 2015

@author: zhuoqinyu
"""
from classdata_mol2 import *
class Res_atom:
    def __init__(self,atomname,atomnum,residue,resid,x,y,z,csp,amber,sybyl,csptype,metcutoff):
        self.atomname=atomname
        self.atomnum=atomnum
        self.residue=residue
        self.resid=int(resid)
        self.x=float(x)
        self.y=float(y)
        self.z=float(z)
        self.csp=float(csp)
        self.csptype=csptype
        self.predict=canbepredict(csptype)
#        self.predict=canbepredict(atomname,residue)
        self.partner=partner_atom(atomname,residue,resid)
        self.amber=amber
        self.sybyl=sybyl
        self.metcutoff=metcutoff

        
class CSP_terms:
    def __init__(self,efe,hb,rce,anishift):
        self.efe=efe #Elecfieldeffect
        self.hb=hb   #HBeffect (direct hydrogen effect)
        self.rce=rce #Ringcurrenteffect
        self.anishift=anishift #just a number

class Elecfieldeffect:
    def __init__(self,efedepsi,epsilon):
        self.efedepsi=efedepsi
        self.epsilon=epsilon

class HBeffect:
    def __init__(self,targetnuc,pairatom,lengthHA,lengthDA,cosangDHA):
        self.targetnuc=targetnuc
        self.pairatom=pairatom
        self.lengthHA=lengthHA
        self.lengthDA=lengthDA
        self.cosangDHA=cosangDHA
#        self.parm=parm
#        self.parm=HB_parm(targetnuc,pairatom,lengthHA,lengthDA,cosangDHA)
        
class Ringcurrenteffect:
    def __init__(self,ring,G,F,costheta,distance,paradistocm,verdistocm):
        self.ring=ring
        self.G=G
        self.F=F
        self.costheta=costheta
        self.distance=distance
        self.paradistocm=paradistocm
        self.verdistocm=verdistocm

                
class Pro_sybyl:
    def __init__(self,residue,atomname,ambertype,sybyltype):
        self.residue=residue
        self.atomname=atomname
        self.ambertype=ambertype
        self.sybyltype=sybyltype
#        self.predict=canbepredict(atomname,residue)
        
#==============================================================================
#         
#==============================================================================
##AnishiftCoeffdict={"HA":1.63296705,
##         "HN":1.28122878,
##         "N":1.0,
##         "CA":11.92792469,
##         "CO":45.23406152,
##         "CB":34.65732757,
##         "H":1.5,
##         "H1":0.69844375,
##         "H4":1.5,
##         "H5":1.5,
##         "HC":1.01431831,
##         "HP":0.5,
##         "Har":1.09706597,
##         "HS":0.93516782,
##         "HO":0.93516782}
#def Cal_AS_csp(anishift,csptype):
#    global AnishiftCoeffdict
#    return anishift*AnishiftCoeffdict[csptype]
##Targetfactordict={"HA":19.66426867,
##         "HN":17.60884178 ,
##         "N":1.0,
##         "CA":11.92792469,
##         "CO":45.23406152,
##         "CB":34.65732757,
##         "H":17.82552361 ,
##         "H1S2":18.91144413 ,
##         "H1S1":14.6807748,
##         "H1O1":38.45644582,
##         "H1O2":47.92793665 ,
##         "H1N":7.93065427,
##         "H4":22.4584929,
##         "H5":45.77939255,
##         "HC":11.46835127,
##         "HP":13.70722603,
##         "Har":17.98674432,
##         "HS":15.05901746 ,
##         "HO":15.05901746 }
#def F_target_factor(res_atom):
#    global Targetfactordict
#    if res_atom.csptype in Targetfactordict:
#        return Targetfactordict[res_atom.csptype]
#    else:
#        return 0.0  
        
global totringtypes
#=['CNNNN','CCNCN','CCCNN','CCCCS','CCNCO','CCCCCC','CCCCCN','CCCNCN']
def RotateMe(text,mode=0,steps=1):
  # Takes a text string and rotates
  # the characters by the number of steps.
  # mode=0 rotate right
  # mode=1 rotate left
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
          
#def rce_surface_parm(res_atom,ring):
#    if res_atom.csptype in ["HA","HN","N","CA","CB","CO","H","H1","H4","H5","HC","HP"]:
#        for n in range(1,len(ring)+1):
#            newring=RotateMe(ring,1,n)
#            flipring=newring[::-1]
#            if newring in totringtypes:
#                ring=newring
#                break
#            elif flipring in totringtypes:
#                ring=flipring
#                break
#        return surfaceparmdict[ring][res_atom.csptype]
#    else: 
#        return [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
#        
        
#def Cal_RCE_csp(rce): # X is the horizontal distance to the cm ,Y is the vertical distance
#    opt=2
#    if opt==1:# 3d surface fitting parm method
#        rce_sum=0
#        rcelist=[]
#        for key in rce.keys(): #maybe multiple rings
#            params=rce[key].surfaceparm
#            X=rce[key].paradistocm
#            Y=rce[key].verdistocm
#        #    if len(params)==10:    #order 3    
#            rcelist.append(params[0]+params[1]*X+params[2]*Y+params[3]*X**2+params[4]*X*Y+params[5]*Y**2+params[6]*X**3+params[7]*X**2*Y+params[8]*X*Y**2+params[9]*Y**3)           
#            rce_sum=rce_sum+rcelist[-1]
#        #        if len(params)==6:    #order 2
#        #            return X**2*params[0]+Y**2*params[1]+params[2]*X*Y+params[3]*X+params[4]*Y+params[5]
#    #    print(rcelist)
#        return rce_sum
#    elif opt==2:#(1-3costheta**3)/(r**3)
#        rce_sum=0
#        rcelist=[]
#        for key in rce.keys():
#            if rce[key].G!=0.0:
#                I=ringI(rce[key].ring)
#                F=rce[key].F
#                rcelist.append(I*F*(1-(rce[key].costheta)**2*3)/(rce[key].distance**3))
#                rce_sum=rce_sum+rcelist[-1]
#        return rce_sum
#        
            
#epsdict={"HA":10.77194867 ,
#         "HN":14.35949083,
#         "N":62.57000814,
#         "CA":-28.24780097,#-18.90593421,
#         "CO":55.88736235,#32.65422755,
#         "H":27.91951052,
#         "H1S2":15.86024939,
#         "H1S1":3.93873698,
#         "H1O1":10.56349744,
#         "H1O2":13.53684692,
#         "H1N":16.04401555,
#         "H4":6.82131575,
#         "H5":16.63081061,
#         "HC":14.42885469,
#         "HP":13.32538414 ,
#         "Har":11.16700558,
#         "HS":27.54639405,
#         "HO":27.54639405}

#def epsilon(res_atom):
#    global epsdict
#    if res_atom.csptype in epsdict:
#        return epsdict[res_atom.csptype]
#    else:
#        return 0.0  
#def Cal_efe_csp(Elecfieldeffect):
#    return Elecfieldeffect.efedepsi*Elecfieldeffect.epsilon
    
    
#def Cal_HB_csp(HBeffect):
#    if HBeffect.parm["function"]=="2d":
#        return HBeffect.parm["a"]/(HBeffect.lengthHA**3)+HBeffect.parm["b"] #a/r**3 +b
#    elif HBeffect.parm["function"]=="3d":
#        C=HBeffect.parm["Coeff"]
#        x=HBeffect.lengthDA
#        y=HBeffect.cosangDHA
#        return C[4]*x**2+C[5]*y**2+C[3]*x*y+C[1]*x+C[2]*y+C[0]
#    elif HBeffect.parm==0.0:
#        return 0.0
    
def partner_atom(atom,res,resid):
    if atom in ["H","H1","H2","H3","HN1","HN2","HN3","HN"]:
        return ("N",0)
    elif atom[0:2]=="HA":  ## first two letters of the atomname
        return ("CA",0)
    elif atom[0:2]=="HB":
        return ("CB",0)
    elif res in ["LYS","LYN"]:
        if atom[0:2]=="HG":return ("CG",0)
        elif atom[0:2]=="HD": return ("CD",0)
        elif atom[0:2]=="HE": return ("CE",0)
        elif atom[0:2]=="HZ": return ("NZ",0)
    elif res in ["MET"]:
        if atom[0:2]=="HG":return ("CG",0)
        elif atom[0:2]=="HE": return ("CE",0)
    elif res in ["PHE"]:
        if atom[0:2]=="HD":return ("CD"+atom[2],0)
        elif atom[0:2]=="HE": return ("CE"+atom[2],0)
        elif atom[0:2]=="HZ": return ("CZ",0)
    elif res in ["PRO"]:
        if atom[0:2]=="HG":return ("CG",0)
        elif atom[0:2]=="HD": return ("CD",0)
    elif res in ["SER"]:
        if atom[0:2]=="HG":return ("OG",0)
    elif res in ["THR"]:
        if atom[0:2]=="HG" and len(atom)==3:return ("OG"+atom[2],0)
        elif len(atom)==4: return ("CG2",0)
    elif res in ["TRP"]:
        if atom=="HD1": return ("CD1",0)
        elif atom=="HE1": return ("NE1",0)
        elif atom=="HE3": return ("CE3",0)
        elif atom[0:2]=="HZ": return ("CZ"+atom[2],0)
        elif atom=="HH2": return("CH2",0)
    elif res in ["TYR"]:
        if atom[0:2]=="HD":return ("CD"+atom[2],0)
        elif atom[0:2]=="HE": return ("CE"+atom[2],0)
        elif atom=="HH": return ("OH",0)
    elif res in ["VAL"]:
        if "HG" in atom: return ("CG"+atom[atom.index("G")+1],0)
    elif res in ["ARG"]:
        if atom[0:2]=="HG":return ("CG",0)
        elif atom[0:2]=="HD": return ("CD",0)
        elif atom[0:2]=="HE": return ("NE",0)
        elif "HH" in atom:
            if atom[0]=="H": return ("NH"+atom[2],0)
            else: return ("NH"+atom[3],0)
    elif res in ["ASN"]:
        if "HD" in atom:
            if atom[0]=="H": return ("ND"+atom[2],0)
            else: return ("ND"+atom[3],0)
    elif res in ["CYS","CYX","CYM"]:
        if atom=="HG": return ("SG",0)
    elif res in ["ASP","ASH"]:
        if "HD" in atom: return ("OD"+atom[2],0)
    elif res in ["GLU","GLH"]:
        if atom[0:2]=="HG":return ("CG",0)
        elif atom[0:2]=="HE":return ("OE"+atom[2],0)
    elif res in ["GLN"]:
        if atom[0:2]=="HG":return ("CG",0)
        elif "HE" in atom:
            if atom[0]=="H": return ("NE"+atom[2],0)
            else: return ("NE"+atom[3],0)
    elif res in ["HID","HIE","HIP","HIS"]:
        if atom=="HD1": return ("ND1",0)
        elif atom=="HD2": return ("CD2",0)
        elif atom=="HE1": return ("CE1",0)
        elif atom=="HE2": return ("NE2",0)
    elif res in ["ILE"]:
        if "HG" in atom: return ("CG"+atom[atom.index("G")+1],0)
        elif "HD" in atom: return ("CD"+atom[atom.index("D")+1],0)
    elif res in ["LEU"]:
        if atom[0:2]=="HG":return ("CG",0)
        elif "HD" in atom: return ("CD"+atom[atom.index("D")+1],0)
    elif atom=="C":
        return ("O",0)
    elif atom=="N":
        return ("C",-1)
    else: return "NAN"
#==============================================================================
def canbepredict(csptype):
    if csptype in ["HA","HN","H","H1","H4","H5","HC","HP","Har","HO"]:
        return True
    else:
        return False


