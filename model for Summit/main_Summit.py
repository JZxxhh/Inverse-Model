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
# 1. set the localtion of ice core drilling site
lat = 72.56
lon = -38.45
tz = -3 # time zone
theta_min, sza_rep = sza_repartition(lat, lon, tz)

# read input file
dfs = pd.read_excel('input/input_snowpack.xlsx',sheet_name='Summit',index_col=0)

# set snowpack properties
cBC = 1.4e-9 # ppb
cdust = 137e-9  # ppb
cHULIS = 31e-9  # ppb
ssa = 47 # m2/kg
density = 380 # kg/m3

# set other parameters
phi = 0.002 # quantum yield for nitrate photolysis
file_path = 'inputs_tuv/Summit/sum' # input tuv file path
fexp = 0.35 # export fraction factor
ozone = ozone_index(dfs['ozone'].values) # total column ozone
accumulation = dfs['accumulation'].values # snow accumulation rate
epsd = +10 # depositional fractionation factor for d15N

# set initial snowpack conditions
SN_conc = dfs['conc'].values
SN_d15N = dfs['d15N'].values
SN_D17O = dfs['D17O'].values
FP_D17O = dfs['D17O_FP'].values

# run the simulations
df = main(theta_min,sza_rep,accumulation,SN_conc,SN_d15N,SN_D17O,
          cBC,cdust,cHULIS,ssa,density,ozone,file_path,fexp,phi, 
          FP_D17O,fc=0.15)

# convert into monthly result
import datetime as dt
tt = np.array([(dt.datetime(2002,1,1)+dt.timedelta(7*i+3)).month for i in range(52)])
# find the time index for each month
# tidx = []
res = np.zeros((12,6))
FD_res = np.zeros((12,4))
for i in range(12):
    res[i,0] = df['FPRI_mass'][np.where(tt==i+1)].mean()*52 # 52 is 1/dt
    res[i,1] = df['FPRI_mass'][np.where(tt==i+1)].std()
    res[i,2] = df['FPRI_d15N'][np.where(tt==i+1)].mean()
    res[i,3] = df['FPRI_d15N'][np.where(tt==i+1)].std()
    res[i,4] = df['FPRI_D17O'][np.where(tt==i+1)].mean()
    res[i,5] = df['FPRI_D17O'][np.where(tt==i+1)].std()
    # calculate atmosphere varitions
    FD_res[i,0] = df['FD_d15N'][np.where(tt==i+1)].mean()-epsd
    FD_res[i,1] = df['FD_d15N'][np.where(tt==i+1)].std()
    FD_res[i,2] = df['FD_D17O'][np.where(tt==i+1)].mean()
    FD_res[i,3] = df['FD_D17O'][np.where(tt==i+1)].std()    

AT = pd.read_excel('input/input_snowpack.xlsx',sheet_name='Summit_AT',index_col=0)
tt = np.array([i.month for i in AT.index])
at_m = np.zeros((12,4))
for i in range(12):
    at_m[i,0] = AT['d15Ntrue'].values[np.where(tt==i+1)].mean()
    at_m[i,1] = AT['d15Ntrue'].values[np.where(tt==i+1)].std()
    at_m[i,2] = AT['D17Otrue'].values[np.where(tt==i+1)].mean()
    at_m[i,3] = AT['D17Otrue'].values[np.where(tt==i+1)].std() 
    
SN_m = pd.read_excel('input/input_snowpack.xlsx',sheet_name='SN_monthly',index_col=0)
SN_conc = np.genfromtxt('input/SN_conc.txt',skip_header=1)

# the final snowpack nitrate isotopes were loaded separately
with open('input/final_snowpack_results.txt','r') as f:
    snowpack_res = np.genfromtxt(f,skip_header=0)

########################################################################
# now draw primary nitrate flux
fig, (ax1,ax2,ax3) = plt.subplots(nrows=3,sharex=True,figsize=(9,12))
plt.subplots_adjust(hspace=0.)

# subplot1: d15N
# Fpri
ax1.errorbar(np.arange(1,13),res[:,2],res[:,3],fmt='o-',color='r',markersize=8,
             capsize=2,label=r'$\delta^{15}$N_F$_{pri}$')
# snowpack
ax1.errorbar(np.arange(1,13),snowpack_res[1,:],snowpack_res[3,:],fmt='^--',
             color='k',alpha=0.3,markersize=8,
             capsize=2,label=r'$\delta^{15}$N_FA')
ax1.set_ylabel(r'10$^{3} \times$ $\delta^{15}$N', fontsize=14, 
                color='r', ha="center", va="center", labelpad=14)
ax1.legend(loc='upper right',fontsize=14)

# subplot2: D17O
# Fpri
ax2.get_yaxis().set_visible(False)
ax21 = ax2.twinx()
ax21.errorbar(np.arange(1,13),res[:,4],res[:,5],fmt='o-',color='r',markersize=8,
              capsize=2,label=r'$\Delta^{17}$O_F$_{pri}$')
# snowpack
ax21.errorbar(np.arange(1,13),snowpack_res[0,:],snowpack_res[2,:],fmt='^--',
              color='k',alpha=0.3,markersize=8,
              capsize=2,label=r'$\Delta^{17}$O_FA')
ax21.set_ylabel(r'10$^{3} \times$ $\Delta^{17}$O', fontsize=14, 
                color='r', ha="center", va="center", labelpad=14)
ax21.legend(loc='lower right',fontsize=14)

# subplot3: Fpri
ls1 = ax3.bar(np.arange(1,13),res[:,0]/12*1e7,width=0.6,color='coral',
              label=r'F$_{pri}$')
ax3.set_ylim([0,12])
ax3.set_yticks([0,3,6,9,12])
ax3.set_ylabel(r'$\times 10^{-7}$ KgN m$^{-2}$ month$^{-1}$',fontsize=14,color='coral')
ax3.tick_params(axis='y', colors='coral')

xlabels = ['Jan','Feb','Mar','Apr','May','Jun','Jul',
            'Aug','Sep','Oct','Nov','Dec']
ax31 = ax3.twinx()
ls2 = ax31.errorbar(np.arange(1,13),SN_conc[:,1],SN_conc[:,2],fmt='^--',
                    color='k',alpha=0.3,markersize=8,
                    capsize=2,label=r'[NO$_{3}^{-}$]')

lines, labels = ax3.get_legend_handles_labels()
lines2, labels2 = ax31.get_legend_handles_labels()
ax31.legend(lines + lines2, labels + labels2, loc='upper right',fontsize=14)
ax31.set_ylabel(r'ng g$^{-1}$', fontsize=14, 
                color='k', ha="center", va="center", labelpad=14)

ax1.tick_params(top=True,which='both', pad=4, direction='in')
ax2.tick_params(top=True,which='both', pad=4, direction='in')
ax3.tick_params(top=True,which='both', pad=4, direction='in')

ax1.text(0.02,0.9,'(a)',fontsize=16,transform=ax1.transAxes)
ax2.text(0.02,0.9,'(b)',fontsize=16,transform=ax2.transAxes)
ax3.text(0.02,0.9,'(c)',fontsize=16,transform=ax3.transAxes)

ax3.set_xticks(np.arange(1,13))
ax3.set_xticklabels(xlabels,fontsize=14)

plt.savefig('output/Summit_Fpri_results.png',dpi=300,bbox_inches='tight')

##############################################################
# plot2: comparing with the local atmospheric observations
fig, (ax1,ax2) = plt.subplots(nrows=2,sharex=True,figsize=(9,8))
plt.subplots_adjust(hspace=0.)
# d15N plot
ax1.errorbar(np.arange(1,13),at_m[:,0],at_m[:,1],
             # linestyle='-',linewidth=2,
             fmt='*',color='b',markersize=10,
             capsize=2,label='Observed')

ax1.errorbar(np.arange(1,13),FD_res[:,0],FD_res[:,1],fmt='o',color='r',markersize=8,
             capsize=2,label='Modelled')

ax1.errorbar(np.arange(1,13),snowpack_res[1,:],snowpack_res[3,:],fmt='^--',
             color='k',alpha=0.3,markersize=8,
             capsize=2,label=r'Initial snowpack')
ax1.set_ylabel(r'10$^{3} \times$ $\delta^{15}$N(NO$_{3}$$^{-}$)', fontsize=14, 
               color='k', ha="center", va="center", labelpad=14)

ax2.get_yaxis().set_visible(False)
ax21 = ax2.twinx()
ls1 = ax21.errorbar(np.arange(1,13),at_m[:,2],at_m[:,3],fmt='*',color='b',markersize=10,
                    capsize=2,label='Observed atmosphere')

ls2 = ax21.errorbar(np.arange(1,13),FD_res[:,2],FD_res[:,3],fmt='o',color='r',markersize=8,
                    capsize=2,label='Modelled atmosphere')

ls3 = ax21.errorbar(np.arange(1,13),snowpack_res[0,:],snowpack_res[2,:],fmt='^--',
                     color='k',alpha=0.3,markersize=8,
                     capsize=2,label=r'Initial snowpack')
ax21.set_ylabel(r'10$^{3} \times$ $\Delta^{17}$O(NO$_{3}$$^{-}$)', fontsize=14, 
               color='k', ha="center", va="center", labelpad=18)

ax21.legend(loc='lower left',fontsize=12)
ax21.set_ylim([21,37])

# set xlabel
xlabels = ['Jan','Feb','Mar','Apr','May','Jun','Jul',
           'Aug','Sep','Oct','Nov','Dec']

ax1.tick_params(top=True,which='both', pad=4, direction='in')
ax2.tick_params(top=True,which='both', pad=4, direction='in')

ax1.text(0.07,0.9,'(a)',fontsize=16,transform=ax1.transAxes)
ax2.text(0.07,0.9,'(b)',fontsize=16,transform=ax2.transAxes)

ax2.set_xticks(np.arange(1,13))
ax2.set_xticklabels(xlabels,fontsize=14)

plt.savefig('output/Summit_atmospheric_results.png',dpi=300,bbox_inches='tight')
########################################################################









