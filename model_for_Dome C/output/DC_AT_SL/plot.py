# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 01:36:02 2023

@author: lenovo
"""

import matplotlib.pyplot as plt
import numpy as np


DC_AT_ob = np.genfromtxt('DC_AT_ob.txt',skip_header=1)
DC_SK_ob = np.genfromtxt('DC_SK_ob.txt',skip_header=1)
datam = np.genfromtxt('DC_AT_SK.txt',skip_header=0).T

# we use box plot
# two subplots: d15N and D17O
# data order: SK ob, SK model, AT ob, AT model
data1 = []
temp = DC_SK_ob[:,0]
temp = temp[~np.isnan(temp)]
data1.append(temp) # ob sk d15N
data1.append(datam[:,0]) # model sk d15N
temp = DC_AT_ob[:,0]
temp = temp[~np.isnan(temp)]
data1.append(temp) # ob AT d15N
data1.append(datam[:,0]-10)

data2 = []
temp = DC_SK_ob[:,1]
temp = temp[~np.isnan(temp)]
data2.append(temp) # ob sk d15N
data2.append(datam[:,1]) # model sk d15N
temp = DC_AT_ob[:,1]
temp = temp[~np.isnan(temp)]
data2.append(temp) # ob AT d15N
data2.append(datam[:,1])

# plot
label = ['Observed',
         'Modeled',
         'Observed',
         'Modeled']


fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(9, 4))
fig.subplots_adjust(hspace=0,wspace=0.3)

bplot1 = ax[0].boxplot(data1,patch_artist=True, 
              widths = 0.6, meanline=True, showmeans=True)
ax[0].set_xticks([1,2,3,4],label,fontsize=10)
ax[0].set_ylabel(r'$\delta^{15}$N(NO$_{3}^{-}$)(‰)', fontsize = 14)
ax[0].axvline(2.5,linestyle='--',linewidth=1)

bplot2 = ax[1].boxplot(data2,patch_artist=True, 
              widths = 0.6, meanline=True, showmeans=True)
ax[1].set_xticks([1,2,3,4],label,fontsize=10)
ax[1].set_ylabel(r'$\Delta^{17}$O(NO$_{3}^{-}$)(‰)', fontsize = 14)

ax[1].axvline(2.5,linestyle='--',linewidth=1)

colors = ['pink', 'lightgreen']*2
for bplot in (bplot1, bplot2):
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
        
ax[0].text(0.02,0.92,'(a)',fontsize=14,transform=ax[0].transAxes)
ax[1].text(0.02,0.92,'(b)',fontsize=14,transform=ax[1].transAxes)

ax[0].text(0.07,1.03,'Skin layer',fontsize=14,transform=ax[0].transAxes)
ax[0].text(0.57,1.03,'Atmosphere',fontsize=14,transform=ax[0].transAxes)
ax[1].text(0.07,1.03,'Skin layer',fontsize=14,transform=ax[1].transAxes)
ax[1].text(0.57,1.03,'Atmosphere',fontsize=14,transform=ax[1].transAxes)
# ax[1].text(0.02,0.92,'(b)',fontsize=14,transform=ax[1].transAxes)

plt.savefig('sensi.pdf',bbox_inches='tight',dpi=300)
'''
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(9, 4))

# Fixing random state for reproducibility
np.random.seed(19680801)


# generate some random test data
all_data = [np.random.normal(0, std, 100) for std in range(6, 10)]

# plot violin plot
axs[0].violinplot(all_data,
                  showmeans=False,
                  showmedians=True)
axs[0].set_title('Violin plot')

# plot box plot
axs[1].boxplot(all_data)
axs[1].set_title('Box plot')

# adding horizontal grid lines
for ax in axs:
    ax.yaxis.grid(True)
    ax.set_xticks([y + 1 for y in range(len(all_data))],
                  labels=['x1', 'x2', 'x3', 'x4'])
    ax.set_xlabel('Four separate samples')
    ax.set_ylabel('Observed values')

plt.show()

'''










