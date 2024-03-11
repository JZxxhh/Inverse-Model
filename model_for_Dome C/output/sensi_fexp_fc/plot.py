# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 21:49:19 2022

@author: lenovo
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import FormatStrFormatter
import pickle


with open('Sensi_fexp_fc_Summit.pickle','rb') as f:
    res1 = pickle.load(f)
FR_S = res1['flux_ratio']
d15N_S = res1['d15N']
d15N_S[np.where(abs(d15N_S)<0.1)] = 0
D17O_S = res1['D17O']

with open('Sensi_fexp_fc_DomeC.pickle','rb') as f:
    res2 = pickle.load(f)
FR_D = res2['flux_ratio']
d15N_D = res2['d15N']
d15N_D[np.where(abs(d15N_D)<0.1)] = 0
D17O_D = res2['D17O']
D17O_D[np.where(abs(D17O_D)>600)] = 600

X, Y = np.meshgrid(np.linspace(0,1,11),np.linspace(0,1,11))

fig, ax = plt.subplots(nrows=2,ncols=3,figsize=(10,7))
plt.subplots_adjust(hspace=0.4,wspace=0.35)

RF_S = 1-FR_S
RF_S[np.where(abs(RF_S)<0.1)] = 0
CS = ax[0,0].contourf(X,Y,RF_S,cmap='cool')
ax[0,0].grid(visible=None, which='major', axis='both',linestyle='--')
cb = fig.colorbar(ScalarMappable(norm=CS.norm, cmap='cool'), ax=ax[0][0], 
                  ticks=np.linspace(0., 0.15, 5),location='top', pad=0.02)
cb.ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
cb.set_label(label='annual loss fraction',fontsize=12)
ax[0][0].tick_params(axis='both', pad=4, direction='in')
ax[0][0].locator_params(axis='x', nbins=6)
ax[0][0].locator_params(axis='y', nbins=6)
ax[0][0].set_xticklabels(['0','0.2','0.4','0.6','0.8','1'])
ax[0][0].set_yticklabels(['0','0.2','0.4','0.6','0.8','1'])
ax[0][0].scatter(0.15,0.35,marker='*',s=40,color='r')

ax[0][0].set_xlabel(r'$f_{c}$',fontsize=12,labelpad=-0.5)
ax[0][0].set_ylabel(r'$f_{exp}$',fontsize=12,labelpad=-0.5)

CS = ax[0,1].contourf(X,Y,d15N_S,cmap='cool')
ax[0,1].grid(visible=None, which='major', axis='both',linestyle='--')
cb = fig.colorbar(ScalarMappable(norm=CS.norm, cmap='cool'), ax=ax[0][1], 
                  ticks=np.linspace(0, 8, 5),location='top', pad=0.02)
cb.ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
cb.set_label(label=r'$\delta^{15}$N (FA-F$_{pri}$), per mil',fontsize=12)
ax[0][1].tick_params(axis='both', pad=8, direction='in')
ax[0][1].locator_params(axis='x', nbins=6)
ax[0][1].locator_params(axis='y', nbins=6)
ax[0][1].scatter(0.15,0.35,marker='*',s=40,color='r')
ax[0][1].set_xlabel(r'$f_{c}$',fontsize=12,labelpad=-0.5)
ax[0][1].set_ylabel(r'$f_{exp}$',fontsize=12,labelpad=-0.5)


CS = ax[0,2].contourf(X,Y,D17O_S,cmap='cool')
ax[0,2].grid(visible=None, which='major', axis='both',linestyle='--')
cb = fig.colorbar(ScalarMappable(norm=CS.norm, cmap='cool'), ax=ax[0][2], 
                  ticks=np.linspace(-1.1, 0.1, 5),location='top', pad=0.02)
cb.ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
cb.set_label(label=r'$\Delta^{17}$O (FA-F$_{pri}$), per mil',fontsize=12)
ax[0][2].tick_params(axis='both', pad=8, direction='in')
ax[0][2].locator_params(axis='x', nbins=6)
ax[0][2].locator_params(axis='y', nbins=6)
ax[0][2].scatter(0.15,0.35,marker='*',s=40,color='r')
ax[0][2].set_xlabel(r'$f_{c}$',fontsize=12,labelpad=-0.5)
ax[0][2].set_ylabel(r'$f_{exp}$',fontsize=12,labelpad=-0.5)


CS = ax[1,0].contourf(X,Y,FR_D/100,cmap='cool')
ax[1,0].grid(visible=None, which='major', axis='both',linestyle='--')
cb = fig.colorbar(ScalarMappable(norm=CS.norm, cmap='cool'), ax=ax[1][0], 
                  ticks=np.linspace(1, 100, 5)/100,location='top', pad=0.02)
cb.ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
cb.set_label(label='annual loss fraction',fontsize=12)
ax[1][0].tick_params(axis='both', pad=8, direction='in')
ax[1][0].locator_params(axis='x', nbins=6)
ax[1][0].locator_params(axis='y', nbins=6)
# ax[1][0].set_xticklabels([])
ax[1][0].scatter(0.15,0.2,marker='*',s=40,color='r')
# ax[1][0].set_yticklabels([])
ax[1][0].set_xlabel(r'$f_{c}$',fontsize=12,labelpad=-0.5)
ax[1][0].set_ylabel(r'$f_{exp}$',fontsize=12,labelpad=-0.5)


CS = ax[1,1].contourf(X,Y,d15N_D,cmap='cool')

ax[1,1].grid(visible=None, which='major', axis='both',linestyle='--')
cb = fig.colorbar(ScalarMappable(norm=CS.norm, cmap='cool'), ax=ax[1][1],
                  ticks=np.linspace(0, 400, 6),location='top', pad=0.02)
cb.ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
cb.set_label(label=r'$\delta^{15}$N (FA-F$_{pri}$), per mil',fontsize=12)
ax[1][1].tick_params(axis='both', pad=8, direction='in')
ax[1][1].locator_params(axis='x', nbins=6)
ax[1][1].locator_params(axis='y', nbins=6)
# ax[1][1].set_xticklabels([])
ax[1][1].scatter(0.15,0.2,marker='*',s=40,color='r')
ax[1][1].set_xlabel(r'$f_{c}$',fontsize=12,labelpad=-0.5)
ax[1][1].set_ylabel(r'$f_{exp}$',fontsize=12,labelpad=-0.5)
# ax[1][1].set_yticklabels([])

# D17O_D[D17O_D<-200] = -200
CS = ax[1,2].contourf(X,Y,D17O_D,cmap='cool')
ax[1,2].grid(visible=None, which='major', axis='both',linestyle='--')
cb = fig.colorbar(ScalarMappable(norm=CS.norm, cmap='cool'), ax=ax[1][2],
                  ticks=np.linspace(-600, 600, 6),location='top', pad=0.02)
cb.ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
cb.set_label(label=r'$\Delta^{17}$O (FA-F$_{pri}$), per mil',fontsize=12)
ax[1][2].tick_params(axis='both', pad=8, direction='in')
ax[1][2].locator_params(axis='x', nbins=6)
ax[1][2].locator_params(axis='y', nbins=6)
# ax[1][2].set_xticklabels([])
ax[1][2].scatter(0.15,0.2,marker='*',s=40,color='r')
ax[1][2].set_xlabel(r'$f_{c}$',fontsize=12,labelpad=-0.5)
ax[1][2].set_ylabel(r'$f_{exp}$',fontsize=12,labelpad=-0.5)

ax[0][0].text(0.05,0.9,'(a)',fontsize=14,transform=ax[0][0].transAxes)
ax[0][1].text(0.05,0.9,'(b)',fontsize=14,transform=ax[0][1].transAxes)
ax[0][2].text(0.05,0.9,'(c)',fontsize=14,transform=ax[0][2].transAxes)
ax[1][0].text(0.05,0.9,'(d)',fontsize=14,transform=ax[1][0].transAxes)
ax[1][1].text(0.05,0.9,'(e)',fontsize=14,transform=ax[1][1].transAxes)
ax[1][2].text(0.05,0.9,'(f)',fontsize=14,transform=ax[1][2].transAxes)

ax[0][0].text(-0.2, 1.2, "Summit", size=15, rotation=0.,
         ha="right", va="top",
         bbox=dict(boxstyle="square",
                   ec=(1., 0.5, 0.5),
                   fc=(1., 0.8, 0.8),
                   )
         )

ax[1][0].text(-0.2, 1.2, "Dome C", size=15, rotation=0.,
         ha="right", va="top",
         bbox=dict(boxstyle="square",
                   ec=(1., 0.5, 0.5),
                   fc=(1., 0.8, 0.8),
                   )
         )

plt.savefig('sensi_fexp_fc.png',dpi=300,bbox_inches='tight')













