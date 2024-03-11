# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 03:38:49 2022

@author: lenovo
"""

import numpy as np
import matplotlib.pyplot as plt

##############################################################################
# these impurity optical properties are adapted from the TARTES model (Libois et al. 2014)
# used for calculating the gamma value in Zatko et al. (2013) radiative transfer parameterization
class Soot(object):
    """class defining soot"""
    
    density = 1800.0
    
    @classmethod
    def refractive_index_imag(cls, wavelength):
        """return the imaginary part of the refracive index (= absorption) for soot."""
        # index from Chang (1990)
        wl_um=1e6*wavelength
        index_soot_real = 1.811+0.1263*np.log(wl_um)+0.027*np.log(wl_um)**2+0.0417*np.log(wl_um)**3
        index_soot_im = 0.5821+0.1213*np.log(wl_um)+0.2309*np.log(wl_um)**2-0.01*np.log(wl_um)**3
        m_soot = index_soot_real -1j * index_soot_im
        n = (m_soot**2 - 1) / (m_soot**2 + 2)  # absorption cross section of small particles (Bohren and Huffman, 1983)
        return n.imag

class HULIS(object):
    """class defining HULIS"""
    
    density = 1500.0

    @classmethod
    def refractive_index_imag(cls, wavelength):
        """return the imaginary part of the refracive index (= absorption) for HULIS."""

        # HULIS from Hoffer (2006)
        wl_um = 1e6 * wavelength
        m_hulis = 1.67 - 8e17 * 1j * (wl_um * 1e3)**(-7.0639) * 1e3 * cls.density * wl_um * 1e-6 / (4 * np.pi)
        n = (m_hulis**2 - 1) / (m_hulis**2 + 2)

        return n.imag

class CrocusDust(object):
    """class defining dust imaginary refractive index from Muller et al., 2011 (one of the higher bound of dust absorption found in the literrature) and 
    Skiles et al.,2014 (one of lower bound of dust absorption found in the literrature). Muller et al., 2011 is default
    François Tuzet, Marie Dumont, June 2018"""
    density = 2600.0
    
    @classmethod
    def refractive_index_imag(cls, wavelength):
        """return the absorption cross section of small particles (Bohren and Huffman, 1983) for a given type of dust

        :param wavelength: wavelength (in m)
        :param formulation: by default use "muller2011" but "skiles2014" is also available.
        """
        wavelength_interp_dust = [299e-3, 350e-3, 400e-3, 450e-3, 500e-3, 550e-3, 600e-3, 650e-3, 700e-3, 750e-3,
                              800e-3, 900e-3, 1000e-3, 1100e-3, 1200e-3, 1300e-3, 1400e-3, 1500e-3, 1600e-3, 1700e-3, 2501e-3]
        index_dust = {'muller2011': [0.038, 0.0312, 0.0193, 0.011, 0.0076, 0.0048, 0.003, 0.0025, 0.0021,
                                     0.002, 0.0018, 0.0017, 0.0016, 0.0016, 0.0016, 0.0015, 0.0015, 0.0015, 0.0014, 0.0014, 0.0014],
                      'skiles2014': [0.0019, 0.0018, 0.0016, 0.0013, 0.0011, 0.0009, 0.0008, 0.0007, 0.00067, 0.00064,
                                     0.00062, 0.00063, 0.00059, 0.00057, 0.00054, 0.00052, 0.00055, 0.00052, 0.0005, 0.00048, 0.00048]
                      }
        formulation = 'skiles2014' # choose the second one is more close to the results at Summit
        wl_um = 1e6 * wavelength
        index_dust_real = 1.53  # real part of the refractive index
        index_dust_im = np.exp(np.interp(np.log(wl_um),
                                         np.log(wavelength_interp_dust),
                                         np.log(index_dust[formulation])))
        m_dust = index_dust_real - 1j * index_dust_im
        n = (m_dust**2 - 1) / (m_dust**2 + 2)

        return n.imag


class Ice_grain(object):
    """class defining HULIS"""

    @classmethod
    def optical(cls, wavelength, ssa, rho_snow):
        """return the extinction and absorption coeffecient of ice grains"""
        rho_ice = 917
        mice = 1.3366 + 2e-11j
        re = 3/(ssa*rho_ice)
        # extinction coeffecient of snow from Wiscombe and Warren, 1980
        kext_ice = 3*rho_snow/(2*re*rho_ice)
        # signle albedo of ice
        c = 24*np.pi*mice.imag/(wavelength*rho_ice*ssa)
        wn = 0.0611 + 0.17*(mice.real-1.3)
        bn = 1.22 + 0.4*(mice.real-1.3)
        phin = 2/3*bn/(1-wn)
        cw_ice = 1/2*(1-wn)*(1-np.exp(-phin*c))
        kabs_ice = cw_ice*kext_ice   

        return kext_ice, kabs_ice

##############################################################################
# Zatko et al. (2013) ACP snow radiative transfer parameterization
# used to calculate the J(z) profile
# nitrate photolysis cross section
def cross_section_NO3():
    # cross section for 14NO3- and 15NO3- photolysis
    N_A    = 6.02214179E+23 # Avogadro Number
    factor = 1e-6
    # Calculate XS14 using the fit in Ayalneh et al 2014 on Chu et al 2003 278K data
    A = 192.5 # M-1
    C = 34052 # cm-1
    W = 3573  # cm-1
    S = 0.9
    # wavelength range    
    wl         = np.arange(280, 351, 1)
    E          = 1/wl * 1e7
    X          = (E - C) / W
    XSoverE    = (1000*np.log(10)/N_A) * factor * A * (1 - S*X) * np.exp(-X*X * (1 - S*X + 0.5*(S*X)*(S*X)))
    XS14_abs   = XSoverE * E
    
    wrf = 0.01  # width reduction factor
    DC  = -32.5 # cm-1
    
    A = A/(1-wrf)  # M-1
    C = C - DC # cm-1
    W = W*(1-wrf)  # cm-1
    S = 0.895
    
    E          = 1/wl * 1e7
    X          = (E - C) / W
    XSoverE    = (1000*np.log(10)/N_A) * factor * A * (1 - S*X) * np.exp(-X*X * (1 - S*X + 0.5*(S*X)*(S*X)))
    XS15_abs   = XSoverE * E
    return XS14_abs, XS15_abs

global XS14_abs, XS15_abs
XS14_abs, XS15_abs = cross_section_NO3()

def cal_gamma(wavelength,cBC,cdust,chulis,ssa,rho_snow):
    IOP_lis = [Soot,HULIS,CrocusDust]
    Con_lis = [cBC,chulis,cdust]
    kext, kabs = Ice_grain.optical(wavelength,ssa,rho_snow)
    for iop,conc in zip(IOP_lis,Con_lis):
        # print(iop,conc)
        n = iop.refractive_index_imag(wavelength)
        kabs += -6*np.pi*rho_snow/wavelength*conc/iop.density*n 
    return (kext*kabs)**0.5

def cal_corr(theta):
    th = np.array([10, 20, 30, 40, 50, 60, 70, 80, 85])
    co = np.array([1.061, 1.063, 1.063, 1.063, 1.058, 1.047, 1.023, 0.973, 6.993])   
    corr = np.interp(theta, th, co)
    return corr

def cal_Jz_depth_profile(cBC,cdust,chulis,AC_df,AC_dd,phi,ssa,rho_snow,theta_min,
                         theta_max=90):
    nb_sza = theta_max-theta_min+1
    # the radiative transfer parameterization from Zatok et al., 2013
    wavelength = np.arange(280,351)*1e-9 #nm
    gamma = cal_gamma(wavelength,cBC,cdust,chulis,ssa,rho_snow)
    ze = 1/gamma.mean()/0.6
    Iz = np.zeros((nb_sza,71,101)) # array to save the I(z) profile, row: wavelength, column: depth
    z = np.arange(101)*0.01
    for dtheta,theta in enumerate(range(theta_min,theta_max+1)): 
        Corruo = cal_corr(theta)
        uo = np.cos(np.deg2rad(theta))
        # calculate the Fdir and Fdiff 
        Fdiff = AC_df[:,dtheta][:,None]*3.831*np.exp(-0.6*(np.tile(z,(71,1))*gamma[:,None]))
        Fdir = AC_dd[:,dtheta][:,None]*3*(0.577+uo)*Corruo*np.exp(-0.6*(np.tile(z,(71,1))*gamma[:,None]))
        kcor = (0.577+uo)/0.577/uo*Corruo
        if kcor>1e5: # kcor threshold set to 1e5
            kcor = 412.8 # use the kcor in 89 degree
        Fdir[:,0] = AC_dd[:,dtheta]*kcor
        Fdir[:,1] = (Fdir[:,0]+Fdir[:,2])/2
        Iz[dtheta,:,:] = Fdir+Fdiff
    # calculate the Jz profile
    Jz14 = np.zeros((nb_sza,101))
    Jz15 = np.zeros((nb_sza,101))
    for theta in range(nb_sza):
        Jz14[theta,:] = (Iz[theta,:,:].T*XS14_abs).sum(axis=1)*phi
        Jz15[theta,:] = (Iz[theta,:,:].T*XS15_abs).sum(axis=1)*phi
    return ze, Jz14, Jz15