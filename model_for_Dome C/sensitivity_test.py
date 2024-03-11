# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 04:47:01 2022

@author: lenovo
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from model.model import main
from model.function import sza_repartition, ozone_index


##############################################################################
# configure the model
# 1. set the localtion of ice core drilling site
localtion = 'DomeC, Antarctica'
lat = -75.1
lon = 123.35
tz = 8 # time zone
theta_min, sza_rep = sza_repartition(lat, lon, tz)
sza_rep = np.vstack((sza_rep[26:,:],sza_rep[:26,:]))
# 2. read input file and set model parameters
dfs = pd.read_excel('input/input_snowpack.xlsx',sheet_name='DomeC_weekly',index_col=0)
# # # LAI contents in ng g-1 in snow grain
cBC = 10e-9 # ppb
cdust = 0e-9
cHULIS = 0e-9
# # snow property
ssa = 55 # 25 to 90 m2 kg-1, Libos et al., 2015 TC
density = 300
# # nitrate photolysis quantum yield
phi = 0.015
file_path = 'inputs_tuv/DomeC/DC'
# # reformed nitrate export fraction
fexp = 0.2
fc = 0.15
# # varitions of TCO
ozone = ozone_index(dfs['ozone'].values)
# # accumulated snow thichness in m
accumulation = dfs['accumulation'].values
cfactor = accumulation[0]*1e-9*density*14/62 # convert ppb to kgN in each snow layer
# # 3. set initial conditions (the archival snowpack depth profile)
conc0 = 21.2
d15N0 = 273.6
D17O0 = 26

def func(t,t0=26,a=0.50551384,b=2.98887728,sig=4.85371715):
    # shifted Gaussian
    return a+b*np.exp(-(t-t0)**2/sig**2)

temp = func(np.arange(52),t0=35)*21.2*(np.random.random(52)+0.5)
SN_conc = temp/temp.mean()*21.2
SN_d15N = dfs['d15N'].values
SN_D17O = dfs['D17O'].values
FP_D17O = dfs['D17O_FP'].values

import pandas as pd
FP_D17O = pd.Series(FP_D17O).interpolate().values
FP_D17O[np.where(FP_D17O>=35)] = 35

# # 4. run the model
# df = main(theta_min,sza_rep,accumulation,SN_conc,SN_d15N,SN_D17O,
#           cBC,cdust,cHULIS,ssa,density,ozone,file_path,fexp,phi, 
#           FP_D17O,fc=fc)


fexpa = np.linspace(0,1,11)
fca = np.linspace(0,1,11)

FR = np.zeros((11,11))
d15N = np.zeros((11,11))
D17O = np.zeros((11,11))

for i in range(len(fexpa)):
    for j in range(len(fca)):
        fexp = fexpa[i]
        fc = fca[j]
        theta_min, sza_rep = sza_repartition(lat, lon, tz)
        df = main(theta_min,sza_rep,accumulation,SN_conc,SN_d15N,SN_D17O,
                  cBC,cdust,cHULIS,ssa,density,ozone,file_path,fexp,phi, 
                  FP_D17O,fc=fc,printf=False)
        d15N[i,j] = df['FA_d15N']-df['FPRI_d15N_annual']
        D17O[i,j] = df['FA_D17O']-df['FPRI_D17O_annual']
        FR[i,j] = (1-df['FA_mass']/df['FPRI_mass_annual'])*100
        print('(%s, %s) finisded! fc=%.1f, fe=%.1f'%(i,j,fc,fexp))


import pickle
 
res = {'flux_ratio':FR,
        'd15N':d15N,
        'D17O':D17O}
    
with open('Sensi_fexp_fc_DomeC.pickle','wb') as f:
    pickle.dump(res,f)


































