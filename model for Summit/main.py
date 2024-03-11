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

cBC = 1.4e-9 # ppb
cdust = 137e-9
cHULIS = 31e-9
ssa = 47
density = 380

phi = 0.001
file_path = 'inputs_tuv/Summit/sum'
fexp = 0.35

# ze, Jz14, eps = cal_Jz_depth_profile(cBC,cdust,chulis,AC_df,AC_dd,phi,ssa,rho_snow,theta_min)
ozone = ozone_index(dfs['ozone'].values)
accumulation = dfs['accumulation'].values

SN_conc = dfs['conc'].values
SN_d15N = dfs['d15N'].values
SN_D17O = dfs['D17O'].values
FP_D17O = dfs['D17O_FP'].values

df = main(theta_min,sza_rep,accumulation,SN_conc,SN_d15N,SN_D17O,
          cBC,cdust,cHULIS,ssa,density,ozone,file_path,fexp,phi, 
          FP_D17O,fc=0.15)

# plt.plot(1-df['frem'])
# print(df['frem_net'])

# convert into monthly result
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

# at_m[0,1] = 0.2 # since only 1 data point for Jan, use measured atandard error to represent
SN_m = pd.read_excel('input/input_snowpack.xlsx',sheet_name='SN_monthly',index_col=0)

fig, (ax1,ax2,ax3) = plt.subplots(nrows=3,sharex=True,figsize=(10,12))
# ax1.set_xticks(np.arange(13))
###############################################################################
# d15N plot
ax1.errorbar(np.arange(1,13),at_m[:,0],at_m[:,1],fmt='*',color='b',markersize=10,
             capsize=2,label='Observed')

ax1.errorbar(np.arange(1,13),FD_res[:,0],FD_res[:,1],fmt='o',color='r',markersize=8,
             capsize=2,label='Modelled')

ax1.plot(np.arange(1,13),SN_m['d15N'],'k--',linewidth=1,alpha=0.5)
ax1.set_ylabel(r'10$^{3} \times$ $\delta^{15}$N(NO$_{3}$$^{-}$)', fontsize=14, 
               color='k', ha="center", va="center", labelpad=14)

# D17O plot
ax2.get_yaxis().set_visible(False)
ax21 = ax2.twinx()
ls1 = ax21.errorbar(np.arange(1,13),at_m[:,2],at_m[:,3],fmt='*',color='b',markersize=10,
                    capsize=2,label='Observed')

ls2 = ax21.errorbar(np.arange(1,13),FD_res[:,2],FD_res[:,3],fmt='o',color='r',markersize=8,
                    capsize=2,label='Modelled')

ax21.plot(np.arange(1,13),SN_m['D17O'],'k--',linewidth=1,alpha=0.3)
ax21.set_ylabel(r'10$^{3} \times$ $\Delta^{17}$O(NO$_{3}$$^{-}$)', fontsize=14, 
               color='k', ha="center", va="center", labelpad=18)
labels = []
lines = []
for ax in [ax2,ax21]:
    axLine, axLabel = ax.get_legend_handles_labels()
    lines.extend(axLine)
    labels.extend(axLabel)

ax2.legend(lines,labels,loc='upper center',fontsize=14)

# plot FPRI and its isotopes
ax3.bar(np.arange(1,13),res[:,0]*1e6,width=0.6,color='coral')
ax3.set_ylim([0,25])
# plt.yticks(np.arange(0,15,3))
ax3.set_yticks([0,4,8,12,16])
ax3.set_ylabel(r'$\times 10^{-6}$ KgN m$^{-2}$ yr$^{-1}$',fontsize=14,color='coral')
ax3.tick_params(axis='y', colors='coral')

ax31 = ax3.twinx()
ls1 = ax31.plot(np.arange(1,13),res[:,2],'ro-',linewidth=1,markersize=8,label=r'$\delta^{15}$N')
ax31.set_ylim(-20,12)
ax31.set_yticks([-12,-6,0,6])
ax31.set_ylabel(r'10$^{3} \times$ $\delta^{15}$N(NO$_{3}$$^{-}$)', fontsize=14, 
                color='r', ha="center", va="center", labelpad=14)

ax32 = ax3.twinx()
ax32.spines['right'].set_position(("axes", 1.13))
ls2 = ax32.plot(np.arange(1,13),res[:,4],'b*--',linewidth=1,markersize=10,label=r'$\Delta^{17}$O')
ax32.set_ylim(24,35)
ax32.set_yticks([24,28,32,36])
ax32.set_ylabel(r'10$^{3} \times$ $\Delta^{17}$O(NO$_{3}$$^{-}$)', fontsize=14, 
                color='b', ha="center", va="center", labelpad=14)
ax32.tick_params(axis='y', colors='b')

labels = []
lines = []
for ax in [ax31,ax32]:
    axLine, axLabel = ax.get_legend_handles_labels()
    lines.extend(axLine)
    labels.extend(axLabel)

ax32.legend(lines,labels,loc='upper right')

# set xlabel
xlabels = ['Jan','Feb','Mar','Apr','May','Jun','Jul',
           'Aug','Sep','Oct','Nov','Dec']

ax1.tick_params(top=True,which='both', pad=4, direction='in')
ax2.tick_params(top=True,which='both', pad=4, direction='in')
ax3.tick_params(top=True,which='both', pad=4, direction='in')

ax1.text(0.02,0.9,'(a)',fontsize=16,transform=ax1.transAxes)
ax2.text(0.02,0.9,'(b)',fontsize=16,transform=ax2.transAxes)
ax3.text(0.02,0.9,'(c)',fontsize=16,transform=ax3.transAxes)

ax32.set_xticks(np.arange(1,13))
ax32.set_xticklabels(xlabels)
plt.subplots_adjust(hspace=0)

# plt.savefig('output/Summit_results.png',dpi=300,bbox_inches='tight')












