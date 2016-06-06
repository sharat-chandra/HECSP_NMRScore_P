# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 16:33:08 2015

@author: zhuoqinyu
"""

#-------------------------------- coords.py ---------------------------------
#  Contains routines for transforming the coordinates of systems of particles
#  Author: Titus Beu, 2013
#http://phys.ubbcluj.ro/~tbeu/INP/libraries.html
#----------------------------------------------------------------------------
from eigsys import *

#============================================================================
def MovetoCM(m, x, y, z, n, cm_old): 
#----------------------------------------------------------------------------
#  Moves a system of n particles to the CM system
# cm_old is the CM of the system before moving
# CM after moving is (0,0,0)
#----------------------------------------------------------------------------
   mCM = 0e0; xCM = 0e0; yCM = 0e0; zCM = 0e0
   for i in range(1,n+1):
      mCM += m[i]
      xCM += m[i] * x[i]
      yCM += m[i] * y[i]
      zCM += m[i] * z[i]
   xCM /= mCM; yCM /= mCM; zCM /= mCM

   for i in range(1,n+1):
      x[i] -= xCM; y[i] -= yCM; z[i] -= zCM
   cm_old[0]=xCM
   cm_old[1]=yCM
   cm_old[2]=zCM
#============================================================================
def PrincipalAxes(m, x, y, z, n, isort):
#----------------------------------------------------------------------------
#  Rotates a set of n particles to the system of principal axes
#  isort =  1 - highest symmetry axis along x-axis
#  isort = -1 - highest symmetry axis along z-axis
#https://books.google.com/books?id=T13SBQAAQBAJ&pg=PA256&lpg=PA256&dq=principal+inertia+axis++python&source=bl&ots=2KUwGppm8A&sig=xELnf1467C8P8kyTUECTkysQBdc&hl=en&sa=X&ei=yYzrVPuQDo2dyQSD74BA&ved=0CF0Q6AEwCQ#v=onepage&q=principal%20inertia%20axis%20%20python&f=false
#----------------------------------------------------------------------------
   Inert = [[0]*4 for i in range(4)]
   Rot   = [[0]*4 for i in range(4)]
   MomInert = [0]*4

   for i in range(1,4):
      for j in range(1,i+1): Inert[i][j] = 0e0

   for i in range(1,n+1):
      mi = m[i]; xi = x[i]; yi = y[i]; zi = z[i]
      Inert[1][1] += mi * (yi*yi + zi*zi)                    # inertia tensor
      Inert[2][2] += mi * (zi*zi + xi*xi)
      Inert[3][3] += mi * (xi*xi + yi*yi)
      Inert[2][1] -= mi * xi*yi
      Inert[3][1] -= mi * zi*xi
      Inert[3][2] -= mi * yi*zi

   Jacobi(Inert,Rot,MomInert,3)                  # diagonalize inertia tensor
   EigSort(Rot,MomInert,3,isort)          # sort eigenvalues and eigenvectors

   for i in range(1,n+1):                   # rotate system to principal axes
      xi = x[i]; yi = y[i]; zi = z[i]
      x[i] = Rot[1][1] * xi + Rot[2][1] * yi + Rot[3][1] * zi
      y[i] = Rot[1][2] * xi + Rot[2][2] * yi + Rot[3][2] * zi
      z[i] = Rot[1][3] * xi + Rot[2][3] * yi + Rot[3][3] * zi
#   print(Inert)
#   print(Rot)
#   print(MomInert) 

   return MomInert
 #============================================================================
def PA_cab(m, x, y, z, n, isort):
#----------------------------------------------------------------------------
#  The vector of principal axes of a set of n particles c a b represented in atomic xyz coordinates
#  isort =  1 - highest symmetry axis along x-axis
#  isort = -1 - highest symmetry axis along z-axis
#----------------------------------------------------------------------------
   Inert = [[0]*4 for i in range(4)]
   Rot   = [[0]*4 for i in range(4)]
   MomInert = [0]*4

   for i in range(1,4): # 3 D coordinates 
      for j in range(1,i+1): Inert[i][j] = 0e0

   for i in range(1,n+1):
      mi = m[i]; xi = x[i]; yi = y[i]; zi = z[i]
      Inert[1][1] += mi * (yi*yi + zi*zi)                    # inertia tensor
      Inert[2][2] += mi * (zi*zi + xi*xi)
      Inert[3][3] += mi * (xi*xi + yi*yi)
      Inert[2][1] -= mi * xi*yi
      Inert[3][1] -= mi * zi*xi
      Inert[3][2] -= mi * yi*zi

   Jacobi(Inert,Rot,MomInert,3)                  # diagonalize inertia tensor
   EigSort(Rot,MomInert,3,isort)          # sort eigenvalues and eigenvectors
   
   pa_cab=[[0e0,0e0,1e0],[1e0,0e0,0e0],[0e0,1e0,0e0]]#principal axis system and will be transformed back to atomic coordinate in the following loop
   
   for i in range(0,3):                   #  principal axes vectors in atomic xyz coordinates
      xi=pa_cab[i][0]; yi=pa_cab[i][1]; zi=pa_cab[i][2]
      pa_cab[i][0] = Rot[1][1] * xi + Rot[1][2] * yi + Rot[1][3] * zi
      pa_cab[i][1] = Rot[2][1] * xi + Rot[2][2] * yi + Rot[2][3] * zi
      pa_cab[i][2] = Rot[3][1] * xi + Rot[3][2] * yi + Rot[3][3] * zi
#   print(Inert)
#   print(Rot)
#   print(MomInert) 

   return pa_cab