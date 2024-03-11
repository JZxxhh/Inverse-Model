# -*- coding: utf-8 -*-
"""
Created on Sat May 11 15:49:35 2019

@author: pb130
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def transfer_conc(trans,Temp,Press):
    air_conc = (Press * 1E2 / (R_gas * Temp)) * 1E-6 * N_A
    return trans*air_conc*1e-9

def O3_NO_195_308K(T):
    'rate constant in cm3 molecule-1 s-1'

    # O3 + NO (Atkinson et al., 2004, preferred value)
    # T range = 195-308 K
    # Atkinson et al., 2004, ACP Evaluated kinetic and photochemical data for atmospheric chemistry: Volume I - gas phase reactions of Ox, HOx, NOx and SOx species
    # Reaction # 54
    # k = 1.4*E−12 * exp(-1310/T) cm3 molecule−1 s−1   
    return 1.4*10**-12 * np.exp(-1310/T)

def BrO_NO_220_430K(T):
    'rate constant in cm3 molecule-1 s-1'

    # BrO + NO (Atkinson et al., 2007, preferred value)
    # T range = 220-430 K
    # Atkinson et al., 2007, ACP Evaluated kinetic and photochemical data for atmospheric chemistry: Volume III – gas phase reactions of inorganic halogens
    # Reaction # 76
    # k = 8.7*E−12 * exp(260/T) cm3 molecule−1 s−1   
    return 8.7*10**-12 * np.exp(260/T)

def HO2_NO_200_400K(T):
    'rate constant in cm3 molecule-1 s-1'

    # HO2 + NO (Atkinson et al., 2004, preferred value)
    # T range = 200-400 K
    # Atkinson et al., 2004, ACP Evaluated kinetic and photochemical data for atmospheric chemistry: Volume I - gas phase reactions of Ox, HOx, NOx and SOx species
    # Reaction # 45
    # k = 3.6*E−12 * exp(270/T) cm3 molecule−1 s−1   
    return 3.6*10**-12 * np.exp(270/T)

def CH3O2_NO_200_430K(T):
    'rate constant in cm3 molecule-1 s-1'

    # CH3O2 + NO (Atkinson et al., 2006, preferred value)
    # T range = 200-430 K
    # Atkinson et al., 2006, ACP Evaluated kinetic and photochemical data for atmospheric chemistry: Volume II – gas phase reactions of organic species
    # Reaction # 120
    # k = 2.3*E−12 * exp(360/T) cm3 molecule−1 s−1   
    return 2.3*10**-12 * np.exp(360/T)

N_A =                     6.02214179E+23  # molecules.mol^(-1)        Avogadro's number
R_gas   =                     8.3144621       # J·K-1·mol-1               Gas constant
# calculate the D17O of nitrate formed in different weeks of 1 year
# import the monthly radical concentration
loc = ['Summit','Dome C','SP']
res = []
for name in loc:
    df = pd.read_excel('radical.xlsx',sheet_name=name)
    # df inclued OH, O3, HO2, BrO concentration, in unit of molecule/cm-3, ppb, molecule/cm-3, ppb respectively
    # transfer ppb to mole/cm-3
    OH = df['OH'].values
    O3 = df['O3'].values
    HO2 = df['HO2'].values
    BrO = df['BrO'].values
    Temp = df['Temp'].values
    Press = df['Pressure'].values
    O3 = transfer_conc(O3,Temp,Press)
    BrO = transfer_conc(BrO,Temp,Press)
    CH3O2 = 7/3*HO2 # following Erbland et al., 2015 method
    # calculate the alpha value
    AT_alpha = (O3_NO_195_308K(Temp) * O3 + BrO_NO_220_430K(Temp)* BrO) / (O3_NO_195_308K(Temp) 
                * O3 + BrO_NO_220_430K(Temp)* BrO + HO2_NO_200_400K(Temp) * HO2 + CH3O2_NO_200_430K(Temp) * CH3O2)
    N2_conc = 0.78 * Temp * 0.001 * 100000 * N_A / (R_gas * Temp * 1000000)
    k_OH = 3.3E-30 * (Temp/300)**-3.0 * N2_conc
    k_O3 = 1.4E-13 * np.exp(-2470/Temp)
    AT_frac_OH_oxidation = (k_OH * OH) / (k_OH * OH + k_O3 * O3)
    # calculate the D17O of nitrate
    AT_D17O_NO2_PSS = AT_alpha * (1.18 * 22.9 + 6.6)
    AT_D17O_O3_tro   = 1.23 * 22.9 + 9.02 # trotospheric ozone D17O
    AT_D17O_addO = AT_frac_OH_oxidation * 5 + (1 - AT_frac_OH_oxidation) * AT_D17O_O3_tro
    D17O_nitrate = 1/3*AT_D17O_addO + 2/3*AT_D17O_NO2_PSS
    plt.plot(np.arange(1,13),D17O_nitrate,label = name)
    res.append(D17O_nitrate)
    
plt.xlabel('Month')
plt.ylabel('D17O of nitrate')
plt.legend()
res = np.vstack((res[0],res[1],res[2]))
np.savetxt('D17O_of_recycled_nitrate.txt',res)




