# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 02:51:14 2022

@author: lenovo
"""

import numpy as np
import matplotlib.pyplot as plt
from .function import sza_repartition, calculating_f_ary, recover_before_the_first_year
from .Impurities import cal_gamma
##############################################################################
def Loop(mass,d15N,D17O,fary,eps,f_cage=0.15):
    nb_dt = eps.shape[0] # number of time step in one year
    total_year = int(d15N.shape[0]/nb_dt)
    SN_mass = np.zeros((fary.shape[0],nb_dt))
    SN_d15N = np.zeros_like(SN_mass)
    SN_D17O = np.zeros_like(SN_mass)
    # initial conditions
    SN_mass[:,-1] = mass
    SN_d15N[:,-1] = d15N
    SN_D17O[:,-1] = D17O
    # array save FP and FT
    FP_mass = np.zeros_like(SN_mass)
    FP_d15N = np.zeros_like(SN_mass)
    FD_mass = np.zeros((nb_dt,))
    FD_d15N = np.zeros_like(FD_mass)
    FD_D17O = np.zeros_like(FD_mass)
    for week in np.arange(nb_dt)[::-1]:
            # here we include all layers in the photic zone instead of a single year
            # recover snow depth profile before photolysis
            ly = (total_year-1)*nb_dt+week
            temp_mass = SN_mass[:ly+1,week]/(fary[:ly+1,week]+f_cage*(1-fary[:ly+1,week]))
            temp_d15N = SN_d15N[:ly+1,week]-fary[:ly+1,week]*(1-f_cage)*eps[week]*np.log(fary[:ly+1,week])/(fary[:ly+1,week]+
                        f_cage-f_cage*fary[:ly+1,week])
            temp_d15N[np.isnan(temp_d15N)] = 0           
            temp_D17O = SN_D17O[:ly+1,week]/(fary[:ly+1,week]+2/3*f_cage-2/3*f_cage*fary[:ly+1,week])*\
                        (fary[:ly+1,week]+f_cage-f_cage*fary[:ly+1,week])
            # set SN before deposition and photolysis in this time step and FD
            if week != 0:
                SN_mass[:ly,week-1] = temp_mass[1:]
                SN_d15N[:ly,week-1] = temp_d15N[1:]
                SN_D17O[:ly,week-1] = temp_D17O[1:]
            FD_mass[week] = temp_mass[0]
            FD_d15N[week] = temp_d15N[0]
            FD_D17O[week] = temp_D17O[0]
            # set FP in this time step
            FP_mass[:ly+1,week] = temp_mass[:ly+1]-SN_mass[:ly+1,week]
            FP_d15N[:ly+1,week] = temp_d15N-fary[:ly+1,week]*eps[week]*np.log(fary[:ly+1,week])/(1-fary[:ly+1,week])
    FP_d15N[np.isnan(FP_d15N)] = 0
    # calculate the weekly FP:
    FP_mass_week = FP_mass.sum(axis=0)
    FP_d15N_week = (FP_mass*FP_d15N).sum(axis=0)/FP_mass.sum(axis=0)
    FP_d15N_week[np.isnan(FP_d15N_week)] = 0
    df = {}
    df['SN_mass'] = SN_mass
    df['SN_d15N'] = SN_d15N
    df['SN_D17O'] = SN_D17O
    df['FD_mass'] = FD_mass
    df['FD_d15N'] = FD_d15N
    df['FD_D17O'] = FD_D17O
    df['FP_mass'] = FP_mass_week
    df['FP_d15N'] = FP_d15N_week
    return df

def mass_conservation(df,fexp):
    FPRI_mass = np.zeros_like(df['FP_mass'])
    FPRI_d15N = np.zeros_like(df['FP_mass'])
    FPRI_D17O = np.zeros_like(df['FP_mass'])
    for i in range(df['FP_mass'].shape[0]):
        FPRI_mass[i] = df['FD_mass'][i]-df['FP_mass'][i-1]*(1-fexp)
        FPRI_d15N[i] = (df['FD_mass'][i]*df['FD_d15N'][i]-
                        df['FP_mass'][i-1]*df['FP_d15N'][i-1]*(1-fexp))/FPRI_mass[i]
        FPRI_D17O[i] = (df['FD_mass'][i]*df['FD_D17O'][i]-
                        df['FP_mass'][i-1]*df['FP_D17O'][i-1]*(1-fexp))/FPRI_mass[i]        
    df['FPRI_mass'] = FPRI_mass
    df['FPRI_d15N'] = FPRI_d15N
    df['FPRI_D17O'] = FPRI_D17O
    df['FP_mass_annual'] = df['FP_mass'].sum()
    df['FP_d15N_annual'] = (df['FP_mass']*df['FP_d15N']).sum()/df['FP_mass_annual'] 
    df['FP_D17O_annual'] = (df['FP_mass']*df['FP_D17O']).sum()/df['FP_mass_annual'] 

##############################################################################
# main function used in one time calculation
def main(theta_min,sza_rep,accumulation,SN_conc,SN_d15N,SN_D17O,
         cBC,cdust,cHULIS,ssa,density,ozone,file_path,fexp,phi, 
         FP_D17O, fc=0.15, printf = True):
    '''
    Main function of the inverse model
    Parameters
    ----------
    theta_min : int
        DESCRIPTION: minimax solar zenith angle in degree
    sza_rep : np.array
        DESCRIPTION: sza repartition array in each ts in minutes
    accumulation : np.array
        DESCRIPTION: accumulated snow thickness in each ts, unit in m
    SN_conc : np.array
        DESCRIPTION: snow niatrate concentration in each snow layer. unit in ng g-1
    SN_d15N : np.array
        DESCRIPTION: snow niatrate d15N in each snow layer. unit in permil
    SN_D17O : np.array
        DESCRIPTION: snow niatrate D17O in each snow layer. unit in permil
    cBC : float
        DESCRIPTION: BC cocnentration in snowpack. unit in ng.g-1
    cdust : float
        DESCRIPTION: dust cocnentration in snowpack. unit in ng.g-1
    cHULIS : float
        DESCRIPTION: HULIS cocnentration in snowpack. unit in ng.g-1
    ssa : float
        DESCRIPTION: specific surface area of snow grain. unit in m2.kg-1
    density : float
        DESCRIPTION: snow density. unit in kg m-3
    ozone : np.array
        DESCRIPTION: total colum ozone in each time step. unit in DU
    file_path : string
        DESCRIPTION: flie path that save the output from tuv model
    fexp : float
        DESCRIPTION: export fraction or photolysis nitrate flux. dimensionless
    phi : float
        DESCRIPTION: quantum yield of nitrate photolysis
    FP_D17O : np.array
        DESCRIPTION: D17O of locally reformed nitrate   
    fc : float, optional
        DESCRIPTION. cage effect fraction during nitrate photolysis. The default is 0.15.

    Returns : dictionary save all nitrate flux in one year
    -------
    None.

    '''
    # step 1: set the total number of yaer
    wl = np.arange(280,350,1)*1e-9
    gamma = cal_gamma(wl,cBC,cdust,cHULIS,ssa,density)
    ze = 1/(0.6*gamma.mean()) # e_folding depth at 305 nm
    annual_accumulation = np.atleast_1d(accumulation).sum()
    nb_dt = accumulation.shape[0]
    total_year = 10+int(3*ze/annual_accumulation)
    print('ze: %.2f m, accu %.2f m'%(ze,annual_accumulation))
    print('total loop year: %s'%total_year)
    
    # step 2: make fary and eps array
    fary, eps, Jsurface = calculating_f_ary(total_year,nb_dt,ozone,cBC,cdust,cHULIS,accumulation,
                                            phi,ssa,density,sza_rep,theta_min,file_path)
    
    # step 3: make initial snowpack depth profile
    mass_temp = SN_conc*1e-9*accumulation*1*density*14/62 # kgN
    d15N_temp = SN_d15N.copy()
    D17O_temp = SN_D17O.copy()
    # used to calculate FA
    SN_mass = np.tile(mass_temp[::-1],total_year)
    # print(SN_mass.shape)
    SN_d15N = np.tile(SN_d15N[::-1],total_year)
    SN_D17O = np.tile(SN_D17O[::-1],total_year)

    # step 4: recover before the first year
    SN_mass,SN_d15N,SN_D17O = recover_before_the_first_year(fary,eps,SN_mass,SN_d15N,SN_D17O,fc)
    
    # step 5: recover snow profile in the first year
    df = Loop(SN_mass,SN_d15N,SN_D17O,fary,eps,f_cage=fc)
    # df['SN_mass'] = SN_mass
    df['FP_D17O'] = FP_D17O
    df['eps'] = eps
    df['Js'] = Jsurface
    df['fary'] = fary
    # stpe 6: solve mass balance equations for atmosphere box
    mass_conservation(df,fexp)
    
    # statistic
    # remaining fraction
    df['frem'] = mass_temp/df['FD_mass']
    # net remaining fraction
    df['frem_net'] = mass_temp.sum()/df['FPRI_mass'].sum()
    # plt.plot()
    df['SN_mass_temp'] = mass_temp
    df['FA_mass'] = mass_temp.sum()
    df['FA_d15N'] = (mass_temp*d15N_temp).sum()/mass_temp.sum()
    df['FA_D17O'] = (mass_temp*D17O_temp).sum()/mass_temp.sum()
    df['FPRI_mass_annual'] = df['FPRI_mass'].sum()
    df['FPRI_d15N_annual'] = (df['FPRI_mass']*df['FPRI_d15N']).sum()/df['FPRI_mass_annual']
    df['FPRI_D17O_annual'] = (df['FPRI_mass']*df['FPRI_D17O']).sum()/df['FPRI_mass_annual']
    # for FD:
    df['FD_mass_annual'] = df['FD_mass'].sum()
    df['FD_d15N_annual'] = (df['FD_mass']*df['FD_d15N']).sum()/df['FD_mass_annual']
    df['FD_D17O_annual'] = (df['FD_mass']*df['FD_D17O']).sum()/df['FD_mass_annual']
    # print summary
    if printf:
        print('################## SUMMARY ####################')
        print('mean FT_mass: %s kgN.m-2.yr-1'%df['FPRI_mass_annual'])
        print('mean FT_d15N: %.3f per mil'%(df['FPRI_d15N_annual']))
        print('mean FT_D17O: %.3f per mil'%(df['FPRI_D17O_annual']))
        print('mean alteration in FT d15N: %.3f'%(df['FA_d15N']-df['FPRI_d15N_annual']))
        print('mean alteration in FT D17O: %.3f'%(df['FA_D17O']-df['FPRI_D17O_annual']))
        print('annual net loss of FPRI: %.3f'%(1-df['frem_net']))
        print('mean alteration in FD d15N: %.3f'%(df['FA_d15N']-df['FD_d15N_annual']))
        print('mean alteration in FD D17O: %.3f'%(df['FA_D17O']-df['FD_D17O_annual']))
        print('annual loss of FD: %.3f'%(1-df['FA_mass']/df['FD_mass_annual']))
    return df





