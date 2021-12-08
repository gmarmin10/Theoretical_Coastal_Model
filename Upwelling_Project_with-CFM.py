#!/usr/bin/env python
# coding: utf-8

# #### Modeling the elemental stoichiometry of phytoplankton and surrounding surface waters in and upwelling or estuarine system
# >Steps to complete project:
# >1. Translate matlab physical model into python
# >2. Substitute Dynamic CFM into model for eco component
# >3. Analyze the results
# 

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math as m
from math import pi
import time
import pylab as lab
from Solver_3D import *
from Solver_2D import *
from af001_energy_calculation import *


# In[2]:


#create the time grid
sperd=60*60*24   #seconds per day
spery=365*sperd  #seconds per year

nyr=1            #number of years to simulate. This can be chagned by the user
T=spery*nyr      #length of time in seconds to simulate
dt=spery/720      #time step
t=np.arange(0,T+dt,dt) #create a time array 
Nt=T/dt+1         #number of time grid points
Nt=int(Nt)


# In[3]:


#create the spatial grid
Lx=5e5           #length in meters
dx=1e3           #grid spacing
Nx=(Lx/dx)+1      #number of grid points
Nx=int(Nx)
xx=np.arange(0,Lx+dx,dx)  #create a space array


# In[4]:


#wind forcing
U0=0.1     #velocity (m/s) subject to change


# In[5]:


#initial conditions of nutrient variables NO3 (Nnut) and PO4 (Pnut)
Pnut=2*np.ones(Nx)  #creating a ones array the size of the number of spatial grid points
Pnut_i=Pnut[0]      #initial value of phosphorus at the coastal boundary
Rup_Po=10*Pnut_i/spery #baseline P uptake rate (will be changed with CFM-Phyto)
Nnut=Pnut*15        #assume NO3 concentration higher than PO4. Change based on field observations
Nnut_i=Pnut_i*15     #intial value of NO3 available at the coastal boundary
Rup_No=10*Nnut_i/spery  #baseline N uptake rate (will be changed with CFM-Phyto)


# In[6]:


# #initial condition of biomass variables- to be replaced with CFM later
# Pbio=0.01*Pnut
# Pbio_i=0.01*Pnut_i    #phytoplankton P at coast (boundary condition)
# Nbio=0.01*Nnut
# Nbio_i=0.01*Nnut_i    #phytoplankton N at coast (boundary condition)


# In[7]:


#initial biological parameters
Kp=0.1   #half-saturation constant for Pnut
Kn=1     #half-saturation constant for Nnut
mu= 1/sperd  #growth rate per sec
phi=0.5     #fraction of uptake remineralized locally
Snp=16     #redfield N:P ratio of phytoplankton
deathd=1e-1 #(d-1) death rate per day
death=deathd/86400  #(s-1) death rate
m2=death   


# In[8]:


#Dynamic CFM initial conditions
x0=1e9
Nbio_i=2.5e-3         #(mol N m-3) original N concentration
Pbio_i=2.5e-2/16      #(mol P m-3) original P concentration


# In[9]:


#Light energy computations
E3=evalue()
E=E3.E
Qc=1.015*10**(-12)/12      #(molC/cell) biomass C per cell (196-18)(average of N and P limited cases from Healey 1985)
YchlN_C=4/55


# In[10]:


#Dynamic CFM starting parameters
Nconst_protein=5.4e-15*0.8
Nstore_max=5.4e-15*0.6             #(molN cell-1) Maximum nitrogen storage (193-25)
Qp_max=25/(3.097e16)               #(mol P cell-1) maximum phosphorus quota
Pconst_other=4.5e-17               #(molP cell-1) Constant part of phosphorus (193-26) * This includes ATP ADP, Phospholipid and DNA RNA
#================================
Molar_mass_DNA_AT_average=307.47        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
Molar_mass_DNA_CG_average=307.97        #(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")

Molar_mass_RNA_AT_average=316.47        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
Molar_mass_RNA_CG_average=323.97        #(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
#================================
#E coli
#================================
CG_Ecoli=0.506          #(dimensionless) from [http://www.ncbi.nlm.nih.gov/genome/167 (accessed 06/18/2016)]
AT_Ecoli=1-CG_Ecoli     #(dimensionless) 

Molar_mass_DNA_Ecoli=Molar_mass_DNA_AT_average*CG_Ecoli+Molar_mass_DNA_CG_average*AT_Ecoli     #(g mol-1) Molar mass of DNA unit
Molar_mass_RNA_Ecoli=Molar_mass_RNA_AT_average*CG_Ecoli+Molar_mass_RNA_CG_average*AT_Ecoli     #(g mol-1) Molar mass of RNA unit

RNA_DNA_mass_ratio=20/7.6   #(ug/ug) Bremer and Dennis 1996
RNA_DNA_mass_ratio=17.844/6.5239  #(ug/ug) from values ad D=0 "07 Bremer and Dennis 1996 data plot.xlsx"

RNA_DNA_molar_ratio=RNA_DNA_mass_ratio/Molar_mass_RNA_Ecoli*Molar_mass_DNA_Ecoli    #(mol mol-1)
#================================
#Stoichiometric parameters
#================================
YcyanoC_N=2                             #(molC molN) C/N molar ratio of cyanophycin
YpgC_P=40                           #(molC molP) C/P molar ratio of PG: Phosphatidyl glycerol (assuming C16 fatty acids (since actually mostly C16 (Huflejt et al., 1990)

CG=0.5755                   #GC%    [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]
YnucacidP_N=1/(3.5*(1-CG)+4*CG)               #(molP molN-1) P/N molar ratio of RNA (193-26) values (193-28) excel file "08 N to P ratio in DNA and RNA.xlsx"

YdnaC_N=3.5*(1-CG)+2.5*CG       #(molC molN-1) C/N molar ratio of dna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)
YrnaC_N=3.25*(1-CG)+2.5*CG      #(molC molN-1) C/N molar ratio of rna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)

DNAmb=2.65354                   #(Mb) Megabase pair of synechococcus DNA in mega (million) base pairs [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]
Avogadro=6.022*10**23           #(molecules mol-1) Avogadro constant
Pdna_const=DNAmb*2*10**6/Avogadro                #(molP cell-1) Constant part of DNA in phosphorus 
Prna_const=Pdna_const*RNA_DNA_molar_ratio       #(molP cell-1) Constant part of RNA in phosphorus
#* Make sure to multiply by 2 as they are base PAIRs"
Ndna_const=Pdna_const/YnucacidP_N      #(molN cell-1) Constant part of DNA in nitrogen
Nrna_const=Ndna_const*RNA_DNA_molar_ratio   #(molN cell-1) Constatn part of RNA in nitrogen
Ndna=Ndna_const    #(molN cell-1) DNA in nitrogen (here assuming constant)
Pdna=Ndna*YnucacidP_N   #(mol P cell-1) DNA in phosphorus
#=======================================
#Calculation of carbon usage (195-16)
#=======================================
CNprotein=4.22945   #(molC molN) the ratio of C to N in protein (derived from Brown 1991) calculation in "13 Amino acid composition of different phytoplankton.xlsx"
#ooooooooooooooooooooooooooooooooooo
#Photosynthetic parameters
#ooooooooooooooooooooooooooooooooooo
I=64    #(umolE m-2 s-1) light intensity
m=5*10**(-19)           #(mol C s-1 cell-1) maintenance carbonhydrate consumption (idea from 172-7)
Pmax0=7                     #(g C /(g Chl h) Maximum production rate per chlorophyll (around 6 by Cullen 1990)
Mchl=893.49             #(g / mol chlorophyll) mollar mass of chlorophyll
Pmax=Pmax0*Mchl/12/3600/55       #(mol C s-1 mol chl-1) carbon fixing rate (156-10) (156-15) for unit conversion)
O=0.0025
T=5
Pchl=Pmax*(1-exp(-O*I*T)) #(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)(193-25)
#================================
#Constant parameters
#================================
Ynphoto_chl=3.2           #((molN cell-1)/(molC chl cell-1)) the stoichiometric ratio for cell photosynthetic enzyme (Rubisco etc.) nitrogen to chlorophyll (193-25)
Ypthylakoid_chl=0.03          #((molP cell-1)/(molC chl cell-1)) the shoichiometric ratio for cell phosphorus in thylakoid membrane to chlorophyll (193-26)
Cnbiosynth=5e-10          #(molN cell-1 s) Constant for varible part of biosynthesis protein nitrogen (193-37)
Cnrna_variable=480*10*1.3         #(s) Constant for Variable part of RNA (193-26)
Qp_max=25/(3.097e16)
Cessential=8.33e-14/10*0.8          #(molC cell-1) essential carbon (lipid membrane, etc.) *8.33e-14/10 is 10%


# In[11]:


#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#Setting arrays for each variable
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
zero=zeros(size(xx))  
D=copy(zero)
Nbio=copy(zero)
Vn=copy(zero)
Pbio=copy(zero)
Vp=copy(zero)

A=copy(zero)
B=copy(zero)
G=copy(zero)
H=copy(zero)
I=copy(zero)
J=copy(zero)
K=copy(zero)
L=copy(zero)
M=copy(zero)
N=copy(zero)
O=copy(zero)
P=copy(zero)
R=copy(zero)
S=copy(zero)

aN=copy(zero)
bN=copy(zero)
cN=copy(zero)
dN=copy(zero)
Dn=copy(zero)

aP=copy(zero)
bP=copy(zero)
cP=copy(zero)
dP=copy(zero)
Dp=copy(zero)

Qn=copy(zero)
Qn_max=copy(zero)
Qp=copy(zero)
Qp_test=copy(zero)
Qn_test=copy(zero)


ls=copy(zero)
Chl=copy(zero)
Nchl=copy(zero)
Nphoto=copy(zero)
Nbiosynth=copy(zero)
Nprotein=copy(zero)
Nrna_variable=copy(zero)
Nrna=copy(zero)
Pthylakoid=copy(zero)
Prna=copy(zero)

Nstore=copy(zero)
Pstore=copy(zero)
x=copy(zero)

Ntot=copy(zero)
Ptot=copy(zero)
Which=copy(zero)


# In[12]:


#initial condition
x[0]=x0
Nbio[0]=Nbio_i         
Pbio[0]=Pbio_i


# In[13]:


#==============================
#Nutrient related
#==============================
Vnmax=3.1467e-19      #(mol N cell-1 s-1) maximum nitrogen uptake rate   (estimated from Healey 1985 based on maximum growth rate of 2 per day) (value from "01 Vmax estimation.xlsx")
Kno3=2.25e-4      #(mol N m-3) Half saturation constant for NO3 uptake (average of oceanic species from table 2 in Eppley 1969)
 
Vpmax=Vnmax/16     #(mol P cell-1 s-1) maximum phosphorus uptake rate (estimated from Vnmax based on the Redfield ratio) 
Kpo4=Kno3/16         #(mol P m-3) Half saturation constant for PO4 uptake (estimated from Kno3 based on the Redfield ratio)


# In[14]:


Loop_array=arange(0,np.size(xx),1)     #array for loops
for i in Loop_array:
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Main calculation 199-21~
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    Vn[i]=Vnmax*(Nbio[i])/(Nbio[i]+Kno3)       #(mol N cell-1 s-1) nitrogen uptake per cell  
    Vp[i]=Vpmax*(Pbio[i])/(Pbio[i]+Kpo4)       #(mol P cell-1 s-1) phosphorus uptake per cell  
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For obtaining D
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #=====================
    #N limited 199-21
    #=====================
    A[i]=((1+E)*Qc*Ynphoto_chl)/Pchl+Cnbiosynth
    B[i]=Nconst_protein+(m*Ynphoto_chl)/Pchl
    G[i]=((1+E)*Qc*YchlN_C)/Pchl
    H[i]=(m*YchlN_C)/Pchl
    I[i]=A[i]*Cnrna_variable
    J[i]=A[i]+B[i]*Cnrna_variable+G[i]
    K[i]=B[i]+Nrna_const+Ndna+H[i]
    
    #if i>-1:
    if i==0:
        aN[i]=I[i]
        bN[i]=J[i]
        cN[i]=K[i]
        dN[i]=-Vn[i]
    else:
        aN[i]=I[i]
        bN[i]=J[i]+I[i]/dt
        cN[i]=K[i]+J[i]/dt
        dN[i]=K[i]/dt-Qn[i-1]/dt-Vn[i]
    
    #-----------------------
    #Solving 3D equation
    #-----------------------
    aNf=float(aN[i])
    bNf=float(bN[i])
    cNf=float(cN[i])
    dNf=float(dN[i])   
    
    Dn0=solver_3D(aNf,bNf,cNf,dNf)
    Dn[i]=Dn0.X

    #======================
    #P limited 199-22
    #======================
    L[i]=((1+E)*Qc*Ypthylakoid_chl)/Pchl
    M[i]=(m*Ypthylakoid_chl)/Pchl
    N[i]=A[i]*Cnrna_variable*YnucacidP_N
    O[i]=L[i]+B[i]*Cnrna_variable*YnucacidP_N
    P[i]=M[i]+Nrna_const*YnucacidP_N+Ndna*YnucacidP_N+Pconst_other
    
    if i==0:
        aP[i]=N[i]
        bP[i]=O[i]
        cP[i]=P[i]
        dP[i]=-Vp[i]
    
    else:
        aP[i]=N[i]
        bP[i]=O[i]+N[i]/dt
        cP[i]=P[i]+O[i]/dt
        dP[i]=(P[i]-Qp[i-1])/dt-Vp[i]        
    
    #-----------------------
    #Solving 3D equation
    #-----------------------
    aPf=float(aP[i])
    bPf=float(bP[i])
    cPf=float(cP[i])
    dPf=float(dP[i])   
    
    Dp0=solver_3D(aPf,bPf,cPf,dPf)
    Dp[i]=Dp0.X
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Which D to apply
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    if Dp[i]>Dn[i]:
        D[i]=Dn[i]
    elif Dn[i]>Dp[i]: 
        D[i]=Dp[i]
    else:
        print("check")      #In case something is not right
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Obtaining output values (Similar to the previous steady state one)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #=========================
    #Chlorophyll related
    #dN[i]========================
    ls[i]=D[i]*Qc   #(molC s-1) Biomass synthesis rate (193-25)
    Chl[i]=((1+E)*ls[i]+m)/Pchl       #(molC chl cell-1) cN[i]hlrophyll concentration (193-25)
    #=========================
    #Nitrogen related
    #========================= 
    Nchl[i]=Chl[i]*YchlN_C         #(molN chl cell-1) Chlorophyll N concentration
    Nphoto[i]=Chl[i]*Ynphoto_chl  #(molN cell-1) Photosynthesis related protein nitrogen (193-25)
    Nbiosynth[i]=D[i]*Cnbiosynth             #(molN cell-1) various part of biosynthesis related protein in N (193-37)
    Nprotein[i]=Nphoto[i]+Nconst_protein+Nbiosynth[i]    #(molN cell-1) All the proteins in N (193-26)
    Nrna_variable[i]=Nprotein[i]*D[i]*Cnrna_variable        #(molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
    Nrna[i]=Nrna_const+Nrna_variable[i]                 #(molN cell-1) nitrogen in RNA (193-26)(193-37)
    
    Qn[i]=Nchl[i]+Nconst_protein+Nphoto[i]+Nbiosynth[i]+Nrna[i]+Ndna+Nstore[i]          #(mol N cell-1) cellular nitrogen quota
    Qn_max[i]=Qn[i]+Nstore_max      #(mol N cell-1) maximum Qn given store amount is maximum
    #=========================
    #Phosphorus related
    #=========================
    Pthylakoid[i]=Chl[i]*Ypthylakoid_chl          #(molP cell-1) Phosphorus in thylakoid membranes: phospholipid, etc. (193-26)
    Prna[i]=Nrna[i]*YnucacidP_N     #(molP cell-1) Phosphorus in RNA 

    Qp[i]=Pthylakoid[i]+Prna[i]+Pdna+Pconst_other    #(molP cell-1)
    
    #==================================
    #Storage calculation 199-5,6,37
    #==================================
    #---------------------
    #Preparation
    #---------------------
    if i==0:
        Qn_test[i]=Vn[i]/D[i]
        Qp_test[i]=Vp[i]/D[i]
    else:
        Qn_test[i]=(Vn[i]+Qn[i-1]/dt)/(D[i]+1/dt)
        Qp_test[i]=(Vp[i]+Qp[i-1]/dt)/(D[i]+1/dt)
    Qn_other=Qn[i]
    Qp_other=Qp[i]
    #---------------------
    #Main part
    #---------------------
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #N limitation -> P store and Fe store  (199-5, 37, a800_04_27)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if Dp[i]>Dn[i]:   
        Which[i]=1      
        #.........................
        #About P store
        #.........................
        if Qp_test[i]<Qp_max:
            Qp[i]=Qp_test[i]
            Pstore[i]=Qp_test[i]-Qp_other
        else:
            Qp[i]=Qp_max
            Pstore[i]=Qp_max-Qp_other

            
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #P limitation -> N store and Fe store  (199-6,37, a800_04_27)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    elif Dn[i]>Dp[i]:
        Which[i]=2       
        #.........................
        #About N store
        #.........................
        if Qn_test[i]<Qn_max[i]:
            Qn[i]=Qn_test[i]
            Nstore[i]=Qn_test[i]-Qn_other  #(molN cell-1) Nitrogen storage in the cell
        else:
            Qn[i]=Qn_max[i]
            Nstore[i]=Nstore_max
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #In case something is not right
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    else:
        print("check")
    #=============================
    #Updating nutrient uptake            #This is necessary not only for updating storiage, but also for making sure the conservation after 3D computation 
    #=============================
    
    if i==0:
        Vn[i]=Qn[i]*D[i]
        Vp[i]=Qp[i]*D[i]
    else:
        Vn[i]=Qn[i]*D[i]+(Qn[i]-Qn[i-1])/dt
        Vp[i]=Qp[i]*D[i]+(Qp[i]-Qp[i-1])/dt    
    
    #==============================
    #Total values
    #==============================
    
    if i==0:
        Ntot[i]=Qn[0]*x[i]+Nbio[i]
        Ptot[i]=Qp[0]*x[i]+Pbio[i]
    else:
        Ntot[i]=Qn[i-1]*x[i]+Nbio[i]
        Ptot[i]=Qp[i-1]*x[i]+Pbio[i]
    #==============================
    #For the next time step
    #==============================
    if i<Loop_array[-1]:
        x[i+1]=x[i]+dt*D[i]*x[i]-dt*m2*x[i]        #(cell m-3) computation of the next step cell density 
        Nbio[i+1]=Nbio[i]-x[i]*Vn[i]*dt+dt*m2*x[i]*Qn[i]      #(mol N m-3) computation of next step NO3 concentration
        Pbio[i+1]=Pbio[i]-x[i]*Vp[i]*dt+dt*m2*x[i]*Qp[i]      #(mol P m-3) computation of next step PO4 concentration
print(Nbio)


# In[15]:


#starting variability with both phytoplanton and environmental plankton
Period=1   #period of oscillation in forcing (velocity) (yr)
w=(2*pi)/Period  #frequency of oscillation
A0=0.5   #amplitude of oscillation
nn=0    #year counter


# In[16]:


it_Nx=np.arange(0,Nx-1,1)
it_Nt=np.arange(0,Nt+1,1)
f=[]
Ua=[]
for n in it_Nt:
    #vary the circulation rates
    for y in t:
        f=A0*(sin(w*y/spery))
        #fn.append(f)
        U0_array=np.full_like(f,U0)
        Ua=U0*f
        U=U0+Ua
    #calculate the biological rates-to be replaced by CFM
    RgrowN=mu*Nbio*(Nnut/(Nnut+Kn))
    RmortN=m2*Nbio**2
    RbioN=RgrowN-RmortN
    RnutN=-RgrowN+phi*RmortN
    RbioP=RbioN/Snp
    RnutP=RnutN/Snp
    
    #update the distribution: Advection scheme
    for i in it_Nx:
        Pnut[i+1]=((dt/dx)*U*Pnut[i]+Pnut[i+1]+RnutP[i]*dt)/(1+dt/dx*U)
        Nnut[i+1]=((dt/dx)*U*Nnut[i]+Nnut[i+1]+RnutN[i]*dt)/(1+dt/dx*U)
        
        Pbio[i+1]=((dt/dx)*U*Pbio[i]+Pbio[i+1]+RbioP[i]*dt)/(1+dt/dx*U)
        Nbio[i+1]=((dt/dx)*U*Nbio[i]+Nbio[i+1]+RbioN[i]*dt)/(1+dt/dx*U)

print((Nbio))


# In[17]:


#some plotting
ax=plt.figure(1)
plt.subplot(2,2,1)
x=np.arange(0,Lx+dx,dx)
x=x*1e-3
plt.plot(x,Pnut,marker='o',color='orange')
plt.xlabel('horizontal distance (km)')
plt.ylabel('PO4 (uM)')
plt.subplot(2,2,2)
plt.plot(x,Nnut,marker='o',color='green')
plt.xlabel('horizontal distance (km)')
plt.ylabel('NO3 (uM)')
plt.subplot(2,2,3)
plt.plot(x,Pbio,marker='o',color='red')
plt.xlabel('horizontal distance (km)')
plt.ylabel('Phyto P (uM)')
plt.subplot(2,2,4)
plt.plot(x,Nbio,marker='o',color='blue')
plt.xlabel('horizontal distance (km)')
plt.ylabel('Phyto N (uM)') 
plt.tight_layout()
plt.savefig('Nutrient_Concentrations.png')
plt.show()


# In[ ]:




