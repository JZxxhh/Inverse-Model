# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 05:57:56 2022

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
# invert the order of the sza_rep if the location is in south hemisphere
sza_rep = np.vstack((sza_rep[26:,:],sza_rep[:26,:])) 

# 2. read input file and set model parameters
dfs = pd.read_excel('input/input_snowpack.xlsx',sheet_name='DomeC_weekly',index_col=0)
# # # LAI contents in ng g-1 in snow grain
cBC = 10e-9 # ppb, adjust to ensure a same e-folding depth as in Erbland et al., 2015 ACP
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
# # varitions of TCO
ozone = ozone_index(dfs['ozone'].values)
# # accumulated snow thichness in m
accumulation = dfs['accumulation'].values
cfactor = accumulation[0]*1e-9*density*14/62 # convert ppb to kgN in each snow layer
# # 3. set initial conditions (the archival snowpack depth profile)
conc0 = 21.2
d15N0 = 273.6
D17O0 = 26

c0 = conc0/70.3
conc = (c0*(0.5+3.5*np.exp(-(np.arange(52)-26)**2/50)))
SN_conc = conc/conc.mean()*21.2 # adjust the archived nitrate profile to mimic the skin layer nitrate varitions

SN_d15N = dfs['d15N'].values
SN_D17O = dfs['D17O'].values
FP_D17O = dfs['D17O_FP'].values

####################################################################################
## the following are the skin layer isotopes, we also tried to adjust the archived 
## nitrate isotopic profile according the SK observations
## and also adjusted the D17O amplitude because it's not as large in snowpack as in SK
## Note it only impact the modeled snowpack nitrate profile (see figure )
d15N_sk = np.array([15.83004485, 14.14813025, 14.37307575, 17.84182805, 19.94740999,
                    34.53661718, 21.07587447, 23.5445841 , 21.44271093, 20.27096172,
                    20.39011093,  9.978935  , 23.63179121, 19.03247211, 14.433153  ,
                    9.8338339 ,  2.05872762, -8.90549566, 11.86067082,  8.64721797,
                    7.34770021,  9.75085098, 12.24917543, 10.35493719, 13.02684763,
                    20.16553186, 23.29838889, 30.10811676, 34.32772736, 37.08270501,
                    40.3997878 , 42.35700784, 11.00711052,  5.19653215,  6.98641179,
                    9.40068196, 14.73208997, 16.00460149, 19.47512879, 22.9456561 ,
                    11.29830708, 20.06466416, 29.25523619, 24.28601039, 19.31678458,
                    27.4979777 , 25.29271788, 18.72517835, 17.92643364, 16.29932632,
                    14.67221901, 14.67221901])

D17O_sk = np.array([34.63309247, 36.39567048, 34.77962323, 42.25174787, 40.57262003,
                    39.5240913 , 42.23202733, 42.94344154, 42.52254749, 39.11241698,
                    40.32379094, 32.74427458, 34.92303164, 35.51508929, 36.10714695,
                    36.6992046 , 34.28978466, 33.04269321, 34.97768862, 31.79045108,
                    30.5788485 , 30.71231837, 28.99848903, 29.90485114, 30.79706323,
                    30.50391055, 30.41277623, 28.86117682, 30.05802659, 30.07896606,
                    29.73016514, 29.92731567, 24.45777062, 28.64451005, 26.93890915,
                    27.19205279, 31.62903678, 28.63394188, 28.9017193 , 29.16949672,
                    32.22148023, 33.46957843, 33.30970218, 33.61645755, 33.92321291,
                    29.95034305, 33.90150166, 34.00501891, 33.94042356, 34.87538666,
                    35.81034977, 35.81034977])

SN_d15N1 = SN_d15N + d15N_sk-d15N_sk.mean() # add the seasonality
SN_D17O1 = SN_D17O + (D17O_sk-D17O_sk.mean())*1/4

import pandas as pd
FP_D17O = pd.Series(FP_D17O).interpolate().values
FP_D17O[np.where(FP_D17O>=35)] = 35

# # 4. run the model
df = main(theta_min,sza_rep,accumulation,SN_conc,SN_d15N1,SN_D17O1,
          cBC,cdust,cHULIS,ssa,density,ozone,file_path,fexp,phi, 
          FP_D17O,fc=0.15)

# convert into monthly results
import datetime as dt
tt = np.array([(dt.datetime(2002,1,1)+dt.timedelta(7*i+3)).month for i in range(52)])
# find the time index for each month
# tidx = []
res = np.zeros((12,6))
FD_res = np.zeros((12,4))
epsd = +10
for i in range(12):
    res[i,0] = df['FPRI_mass'][np.where(tt==i+1)].mean()*52 # 52 is 1/dt
    res[i,1] = df['FPRI_mass'][np.where(tt==i+1)].std()
    res[i,2] = df['FPRI_d15N'][np.where(tt==i+1)].mean()
    res[i,3] = df['FPRI_d15N'][np.where(tt==i+1)].std()
    res[i,4] = df['FPRI_D17O'][np.where(tt==i+1)].mean()
    res[i,5] = df['FPRI_D17O'][np.where(tt==i+1)].std()
    # calculate atmosphere varitions
    FD_res[i,0] = df['FD_d15N'][np.where(tt==i+1)].mean()
    FD_res[i,1] = df['FD_d15N'][np.where(tt==i+1)].std()
    FD_res[i,2] = df['FD_D17O'][np.where(tt==i+1)].mean()
    FD_res[i,3] = df['FD_D17O'][np.where(tt==i+1)].std()    

# # import data
AT_ob = np.genfromtxt('data/atmosphere_monthly.txt')
SK_ob = np.genfromtxt('data/skin_layer_monthly.txt')
SN1 = np.genfromtxt('data/snow pit1.txt',skip_header=1)
SN2 = np.genfromtxt('data/snow pit2.txt',skip_header=1)
SN3 = np.genfromtxt('data/snow pit3.txt',skip_header=1)
SN4 = np.genfromtxt('data/snow pit4.txt',skip_header=1)

# # compared with Joesph's results (Erbland et al., 2015 ACP, the original TRANSITS model output)
dataJ = np.genfromtxt('input/Joseph_snowpack_results.txt')
depthJ = np.linspace(0,1,1000)

########################################################################
# PLOT SNOW PROFILE
file_name = 'output/DomeC_SN_adjusted.png'
fig, (ax1,ax2,ax3) = plt.subplots(1,3,sharey=True,figsize=(9, 5))
plt.subplots_adjust(hspace=0.1)
plt.subplots_adjust(wspace=0.1)

ax11 = ax1.twiny()
ax11.plot(SN1[:,1],-SN1[:,0]/100,'ro--',linewidth=1,markersize=3)
ax11.plot(SN2[:,1],-SN2[:,0]/100,'ro--',linewidth=1,markersize=3)
ax11.plot(SN3[:,1],-SN3[:,0]/100,'ro--',linewidth=1,markersize=3)
ax11.plot(SN4[:,1],-SN4[:,0]/100,'ro--',linewidth=1,markersize=3,label='Observations')
ax11.plot(df['SN_mass'][:286,26]/cfactor,-np.arange(286)*accumulation[0],
          color='k',linewidth=2,label='Inverse model')
ax11.plot(dataJ[0,:500],-depthJ[:500],
          color='b',linewidth=2,label='Forward model')
# ax11.set_xlim([-100,1500])
ax11.legend(loc=('lower right'),fontsize=12)
ax11.set_xlabel(r'[NO$_{3}^{-}$] (ppb)', fontsize=14, color='k', ha="center", va="center", labelpad=14) #set xlabel and ylabel

# #d15N
ax21 = ax2.twiny()
ax21.plot(SN1[:,2],-SN1[:,0]/100,'ro--',linewidth=1,markersize=3)
ax21.plot(SN2[:,2],-SN2[:,0]/100,'ro--',linewidth=1,markersize=3)
ax21.plot(SN3[:,2],-SN3[:,0]/100,'ro--',linewidth=1,markersize=3)
ax21.plot(SN4[:,2],-SN4[:,0]/100,'ro--',linewidth=1,markersize=3)
ax21.plot(df['SN_d15N'][:286,26],-np.arange(286)*accumulation[0],
          color='k',linewidth=2)
ax21.plot(dataJ[1,:500],-depthJ[:500],
          color='b',linewidth=2)
ax21.set_xlabel(r'10$^{3} \times$ $\delta^{15}$N', fontsize=14, color='k', ha="center", va="center", labelpad=14) #set xlabel and ylabel

ax31 = ax3.twiny()
ax31.plot(SN1[:,3],-SN1[:,0]/100,'ro--',linewidth=1,markersize=3)
ax31.plot(SN2[:,3],-SN2[:,0]/100,'ro--',linewidth=1,markersize=3)
ax31.plot(SN3[:,3],-SN3[:,0]/100,'ro--',linewidth=1,markersize=3)
ax31.plot(SN4[:,3],-SN4[:,0]/100,'ro--',linewidth=1,markersize=3)
ax31.plot(df['SN_D17O'][:286,26],-np.arange(286)*accumulation[0],
          color='k',linewidth=2)
ax31.plot(dataJ[2,:500]-3,-depthJ[:500],
          color='b',linewidth=2)
ax31.set_xlabel(r'10$^{3} \times$ $\Delta^{17}$O', fontsize=14, color='k', ha="center", va="center", labelpad=14) #set xlabel and ylabel

ax1.axes.get_xaxis().set_visible(False)
ax2.axes.get_xaxis().set_visible(False)
ax3.axes.get_xaxis().set_visible(False)
ax1.set_ylabel('Depth (m)',fontsize=14)
ax1.text(0.05,0.95,'(a)',fontsize=16,transform=ax1.transAxes)
ax2.text(0.05,0.95,'(b)',fontsize=16,transform=ax2.transAxes)
ax3.text(0.05,0.95,'(c)',fontsize=16,transform=ax3.transAxes)

ax11.axhspan(-0.3, 0, alpha=0.5, color='yellow')
ax21.axhspan(-0.3, 0, alpha=0.5, color='yellow')
ax31.axhspan(-0.3, 0, alpha=0.5, color='yellow')

plt.savefig(file_name,dpi=300,bbox_inches='tight')





