# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 11:15:12 2022

@author: lenovo
"""

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from .Impurities import cal_Jz_depth_profile
import warnings
warnings.filterwarnings("ignore") 

###############################################################################
def cal_sza(ttime,currentDayOfYear,latitude,longitude,timezone=0):
#    ttime = get_second_of_day(ttime)
#    cosx_threshold = 1.0E-30  
    GMT = (ttime/60-timezone+24)%24 # can use numpy array as input
    theta = 2.0 * np.pi * currentDayOfYear / 365.0 
    # The sun declination is the angle between the center of the Sun 
    # and Earths equatorial plane. 
    b0 =  0.006918 
    b1 = -0.399912 
    b2 =  0.070257 
    b3 = -0.006758 
    b4 =  0.000907 
    b5 = -0.002697 
    b6 =  0.001480 
    # solar declination
    dec = b0 + b1 * np.cos(theta) + b2 * np.sin(theta) + b3 * np.cos(2.0 * theta) + b4 * np.sin(2.0 * theta) + b5 * np.cos(3.0 * theta) + b6 * np.sin(3.0 * theta) 
    # The equation of time accounts for the discrepancy between the 
    # apparent and the mean solar time at a given location. 
    c0 = 0.000075 
    c1 = 0.001868 
    c2 = -0.032077 
    c3 = -0.014615 
    c4 = -0.040849 
    eqtime = c0 + c1 * np.cos(theta) + c2 * np.sin(theta) + c3 * np.cos(2.0 * theta) + c4 * np.sin(2.0 * theta)     
    # this is an error in MCM program for calculating lha
    lha = np.pi * (GMT / 12.0 - 1.0 + longitude / 180.0) + eqtime 
    lat = latitude * np.pi / 180.0 
    sinld = np.sin( lat ) * np.sin( dec ) 
    cosld = np.cos( lat ) * np.cos( dec ) 
    # Calculate the cosine of the solar zenith angle. 
    cosx = np.cos( lha ) * cosld + sinld 
    return np.arccos(cosx)*180/np.pi

def sza_repartition(lat, lon, tz, resolution='Weekly'):
    '''
    Make solar zenith angle repartition file for a given site
    Parameters
    ----------
    lat : float, latitude
    lon : float, longtitude
    tz  : time zone of the site
    Returns: None
    -------
    None.

    '''
    sza_lisd = []
    # calculate the sza varitaions in one day
    ttime = np.arange(1440)
    for i in range(365):
        sza_lisd.append(np.round(cal_sza(ttime,i,lat,lon,timezone=tz)))
    szad = np.array(sza_lisd)
    theta_min = int(szad.min())
    # print('minimal solar zenith angle is %s degree'%theta_min)
    theta_max = 90
    sza_lisd = []
    # we choose the resolution according to the accumulation rate
    # the best accumulat velocity is 1-2 cn per time step
    if resolution == 'Weekly': # for normal accumulation rate sites such as Summit, Greenand
        for i in range(52):
            temp = []
            for j in range(theta_min,theta_max+1):
                temp.append(np.where(szad[i*7:(i+1)*7,:]==j)[0].shape[0])
            sza_lisd.append(temp)
    elif resolution == 'Monthly': # for accumulation rate is lower such as DML
        yr = 2002 # any year is good.
        dof = []
        for i in range(1,13):
            dof.append((dt.datetime(yr,i,1)-dt.datetime(yr,1,1)).days)
        dof.append(-1)
        for i in range(12):
            temp = []
            for j in range(theta_min,theta_max+1):
                temp.append(np.where(szad[dof[i]:dof[i+1],:]==j)[0].shape[0])
            sza_lisd.append(temp)    
    elif resolution == 'Bi-monthly': # for accumulation rate is extremly low such as DC Antarctica
        # not recommand
        yr = 2002 # any year is good.
        dof = []
        for i in range(1,12,2):
            dof.append((dt.datetime(yr,i,1)-dt.datetime(yr,1,1)).days)
        dof.append(-1)
        for i in range(6):
            temp = []
            for j in range(theta_min,theta_max+1):
                temp.append(np.where(szad[dof[i]:dof[i+1],:]==j)[0].shape[0])
            sza_lisd.append(temp)     
    SZA_rep = np.array(sza_lisd)
    return theta_min, SZA_rep 

##############################################################################
def calculating_f_ary(nb_yr,nb_dt,ozone,cBC,cdust,chulis,accumulation,
                      phi,ssa,rho_snow,sza_rep,theta_min,file_path):
    # 2 arrays used to save the f and eps
    num = nb_dt*nb_yr
    fary = np.zeros((num,num))
    eps = np.zeros((nb_dt,)) # use surface eps for the whole snowpack
    Js_week = np.zeros((nb_dt,)) # surface J value in each time step
    sza_rep *= 60 # convert to second
    for ts in range(nb_dt):
        # open the tuv file
        wl = np.arange(280,351,1)
        # conversion factor from W m-2 nm-1 to photon cm-2 s-1 nm-1
        cfactor = 5.031e11*wl
        with open(file_path+'%s.txt'%ozone[ts],'r') as f: # toal downward irradiance
            IR_dt = np.genfromtxt(f, skip_header=4,skip_footer=1)[:,1:]
        AC_dt = IR_dt*cfactor[:,None]
        with open(file_path+'0%s.txt'%ozone[ts],'r') as f: # direct downward irradiance
            IR_dd = np.genfromtxt(f, skip_header=4,skip_footer=1)[:,1:]
        AC_dd = IR_dd*cfactor[:,None]
        AC_df = AC_dt-AC_dd
        AC_dd[AC_dd<0] = 0
        _, J_14, J_15 = cal_Jz_depth_profile(cBC,cdust,chulis,AC_df,
                                                AC_dd,phi,ssa,rho_snow,theta_min)
        accu_thickness = np.append(accumulation[:ts+1][::-1],np.tile(accumulation[::-1],nb_yr-1))
        # find the J value of the middle depth for each layer
        mid_depth = (accu_thickness.cumsum()+np.insert(accu_thickness,0,0).cumsum()[:-1])/2
        J_ary14 = np.zeros((J_14.shape[0],mid_depth.shape[0]))   
        for i in range(J_14.shape[0]):
            J_ary14[i,:] = np.interp(mid_depth,np.arange(101)*0.01,J_14[i,:],right=0)
        fary[:len(accu_thickness),ts] = np.exp(-(J_ary14*sza_rep[ts,:][:,None]).sum(axis=0))
        eps[ts] = (np.nansum((J_15[:,0]*sza_rep[ts,:]))/np.nansum((J_14[:,0]*sza_rep[ts,:]))-1)*1000
        # eps_d = J_15/J_14-1
        # eps[ts] = (eps_d[:,0]*sza_rep[ts,:]).sum()/sza_rep[ts,:].sum()*1000
        # print(eps[ts])
        for year in range(1,nb_yr):
            accu_thickness = np.tile(accumulation[::-1],nb_yr)
            mid_depth = (accu_thickness.cumsum()+np.insert(accu_thickness,0,0).cumsum()[:-1])/2+accumulation.sum()*(year-1)+accumulation[:ts+1].sum()
            J_ary14 = np.zeros((J_14.shape[0],mid_depth.shape[0]))   
            for i in range(J_14.shape[0]):
                J_ary14[i,:] = np.interp(mid_depth,np.arange(101)*0.01,J_14[i,:],right=0)
            fary[:,ts+year*nb_dt] = np.exp(-(J_ary14*sza_rep[ts,:][:,None]).sum(axis=0)) 
        Js_week[ts] = (J_14[:,0]*sza_rep[ts,:]).sum()/sza_rep[ts,:].sum()
    fary[fary==0] = 1
    eps[np.isnan(eps)] = 0
    print(eps)
    Js_week[np.isnan(Js_week)] = 0
    return fary, eps, Js_week

##############################################################################
def recover_before_the_first_year(fary,eps,conc,d15N,D17O,f_cage=0.15):
    # make sure all fary and eps value are not nan value
    nbdt = eps.shape[0]
    fary[np.isnan(fary)] = 0
    eps[np.isnan(eps)] = 0
    num = int(fary.shape[0]/eps.shape[0])
    eps = np.tile(eps,(fary.shape[0],num-1))
    # for concentration:
    frem = (fary[:,nbdt:]+(1-fary[:,nbdt:])*f_cage).prod(axis=1)
    conc = conc/frem
    # for d15N:
    sigma_delta = (eps*(1-f_cage)*fary[:,nbdt:]*np.log(fary[:,nbdt:])/(fary[:,nbdt:]+f_cage-fary[:,nbdt:]*f_cage)).sum(axis=1)
    sigma_delta[np.isnan(sigma_delta)] = 0
    d15N -= sigma_delta
    prod_Delta = ((fary[:,nbdt:]+(1-fary[:,nbdt:])*f_cage)/(fary[:,nbdt:]+(1-fary[:,nbdt:])*f_cage*2/3)).prod(axis=1) # exchange of oxygen with H2O
    D17O = D17O*prod_Delta
    return conc,d15N,D17O      

def ozone_index(ozone):
    ozone0 = np.arange(0,1000,25)
    res = np.zeros((len(ozone)))
    for i in range(len(ozone)):
        idx = np.abs(ozone[i]-ozone0).argmin()
        res[i] = ozone0[idx]
    return res.astype(int)