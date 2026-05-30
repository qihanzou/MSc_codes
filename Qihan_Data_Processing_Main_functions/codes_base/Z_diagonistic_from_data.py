# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 19:33:08 2023

@author: Qihan Zou
"""

import numpy as np 
Z_sun = 8.69

if __name__ != '__main__':
    ##2017
    N2_cal_data_2017   = np.loadtxt('C:/Users/qihan/Desktop/diagnositic_data/Curti17_N2.txt'  )
    O3N2_cal_data_2017 = np.loadtxt('C:/Users/qihan/Desktop/diagnositic_data/Curti17_O3N2.txt')
    O3O2_cal_data_2017 = np.loadtxt('C:/Users/qihan/Desktop/diagnositic_data/Curti17_O3O2.txt')
    R2_cal_data_2017   = np.loadtxt('C:/Users/qihan/Desktop/diagnositic_data/Curti17_R2.txt'  )
    R3_cal_data_2017   = np.loadtxt('C:/Users/qihan/Desktop/diagnositic_data/Curti17_R3.txt'  )
    R23_cal_data_2017  = np.loadtxt('C:/Users/qihan/Desktop/diagnositic_data/Curti17_R23.txt' )
    ##2020
    N2_cal_data_2020   = np.loadtxt('C:/Users/qihan/Desktop/diagnositic_data/Curti20_N2.txt'  )
    O3N2_cal_data_2020 = np.loadtxt('C:/Users/qihan/Desktop/diagnositic_data/Curti20_O3N2.txt')
    O3O2_cal_data_2020 = np.loadtxt('C:/Users/qihan/Desktop/diagnositic_data/Curti20_O3O2.txt')
    O3S2_cal_data_2020 = np.loadtxt('C:/Users/qihan/Desktop/diagnositic_data/Curti20_O3S2.txt')
    R2_cal_data_2020   = np.loadtxt('C:/Users/qihan/Desktop/diagnositic_data/Curti20_R2.txt'  )
    R3_cal_data_2020   = np.loadtxt('C:/Users/qihan/Desktop/diagnositic_data/Curti20_R3.txt'  )
    R23_cal_data_2020  = np.loadtxt('C:/Users/qihan/Desktop/diagnositic_data/Curti20_R23.txt' )
    RS32_cal_data_2020 = np.loadtxt('C:/Users/qihan/Desktop/diagnositic_data/Curti20_RS32.txt')
    S2_cal_data_2020   = np.loadtxt('C:/Users/qihan/Desktop/diagnositic_data/Curti20_S2.txt'  )


'''
line ratio:
    R2        [O II]λ3727,29/Hb
    R3        [O III]λ5007/Hb
    N2        [N II]λ6584/Ha
    S2        [S II]λ6717,31/Ha
    R23       ([O II]λ3727,29 + [O III]λ4959,5007)/Hb
    O3O2      [O III]λ5007/[O II]λ3727,29  (*** O32 in Curti+17 is same as O3O2 in 20)
    RS32      [O III]λ5007/Hb + [S II]λ6717,31/Ha
    O3S2      [O III]λ5007/Hb / [S II]λ6717,31/Ha
    O3N2      [O III]λ5007/Hb / [N II]λ6584/Ha

PHANGS-ALMA Data:
No.    Name         Ver    Type      Cards   Dimensions   Format
  0  OII3726         1 PrimaryHDU      23   (246, 666)   float32   
  1  OII3726_ERR     1 ImageHDU        24   (246, 666)   float32   
  2  OII3729         1 ImageHDU        24   (246, 666)   float32   
  3  OII3729_ERR     1 ImageHDU        24   (246, 666)   float32   
  4  HBETA           1 ImageHDU        24   (246, 666)   float32   
  5  HBETA_ERR       1 ImageHDU        24   (246, 666)   float32   
  6  OIII5007        1 ImageHDU        24   (246, 666)   float32   
  7  OIII5007_ERR    1 ImageHDU        24   (246, 666)   float32   
  8  HALPHA          1 ImageHDU        24   (246, 666)   float32   
  9  HALPHA_ERR      1 ImageHDU        24   (246, 666)   float32   
 10  NII6583         1 ImageHDU        24   (246, 666)   float32   
 11  NII6583_ERR     1 ImageHDU        24   (246, 666)   float32   
 12  SII6716         1 ImageHDU        24   (246, 666)   float32   
 13  SII6716_ERR     1 ImageHDU        24   (246, 666)   float32   
 14  SII6731         1 ImageHDU        24   (246, 666)   float32   
 15  SII6731_ERR     1 ImageHDU        24   (246, 666)   float32 
 
 PHANGS-MUSE data:
 Hb
 Ha
 NII
 SII
 OIII
 
 line_df: hdu list that contain all the different emission line data reduced from TYPHOON data cubes
 Z: array, the Metallicity using this diagnostic
 e_Z: array, Error in metallicity using this diagnostic, computed via linear error propagation.
'''

def compute_Z_R2_Curti17(line_df):
    # OII and Hb
    line_IDs  = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_OII3726     = line_df[line_IDs.index('OII3726')].data
    e_OII3726     = line_df[line_IDs.index('OII3726_ERR')].data
    f_OII3729     = line_df[line_IDs.index('OII3729')].data
    e_OII3729     = line_df[line_IDs.index('OII3729_ERR')].data
    f_Hb      = line_df[line_IDs.index('HBETA')].data
    e_Hb      = line_df[line_IDs.index('HBETA_ERR')].data
    # compute line ratios
    R2        = np.log10((f_OII3726+f_OII3729)/f_Hb)
    Z_R2_17 = np.interp(R2,R2_cal_data_2017[:,0],R2_cal_data_2017[:,1],left=np.nan, right=np.inf)
    x_R2_17 = Z_R2_17 - Z_sun
    return Z_R2_17

def compute_Z_R2_Curti20(line_df):
    # OII and Hb
    line_IDs  = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_OII3726     = line_df[line_IDs.index('OII3726')].data
    e_OII3726     = line_df[line_IDs.index('OII3726_ERR')].data
    f_OII3729     = line_df[line_IDs.index('OII3729')].data
    e_OII3729     = line_df[line_IDs.index('OII3729_ERR')].data
    f_Hb      = line_df[line_IDs.index('HBETA')].data
    e_Hb      = line_df[line_IDs.index('HBETA_ERR')].data
    # compute line ratios
    R2        = np.log10((f_OII3726+f_OII3729)/f_Hb)
    Z_R2_20 = np.interp(R2,R2_cal_data_2020[:,0],R2_cal_data_2020[:,1],left=np.nan, right=np.inf)
    x_R2_20 = Z_R2_20 - Z_sun
    return Z_R2_20

def compute_Z_R3_Curti17(line_df):
    # need R3=[O III]λ5007/Hb; OIII5007
    line_IDs  = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_OIII5007     = line_df[line_IDs.index('OIII5007')].data
    e_OIII5007     = line_df[line_IDs.index('OIII5007_ERR')].data
    f_Hb      = line_df[line_IDs.index('HBETA')].data
    e_Hb      = line_df[line_IDs.index('HBETA_ERR')].data
    # compute line ratios
    R3 = np.log10(f_OIII5007/f_Hb)
    Z_R3_17 = np.interp(R3,R3_cal_data_2017[:,0],R3_cal_data_2017[:,1],left=np.nan, right=np.inf)
    x_R3_17 = Z_R3_17 - Z_sun
    return Z_R3_17

def compute_Z_R3_Curti20(line_df):
    # need R3=[O III]λ5007/Hb; OIII5007
    line_IDs  = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_OIII5007     = line_df[line_IDs.index('OIII5007')].data
    e_OIII5007     = line_df[line_IDs.index('OIII5007_ERR')].data
    f_Hb      = line_df[line_IDs.index('HBETA')].data
    e_Hb      = line_df[line_IDs.index('HBETA_ERR')].data
    # compute line ratios
    R3 = np.log10(f_OIII5007/f_Hb)
    Z_R3_20 = np.interp(R3,R3_cal_data_2020[:,0],R3_cal_data_2020[:,1],left=np.nan, right=np.inf)
    x_R3_20 = Z_R3_20 - Z_sun
    return Z_R3_20

def compute_Z_N2_Curti17(line_df):
    #  N2=[N II]λ6584/Ha, NII6583
    line_IDs  = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_NII6583     = line_df[line_IDs.index('NII6583')].data
    e_NII6583     = line_df[line_IDs.index('NII6583_ERR')].data
    f_Ha      = line_df[line_IDs.index('HALPHA')].data
    e_Ha      = line_df[line_IDs.index('HALPHA_ERR')].data
    # compute line ratios
    N2 = np.log10(f_NII6583/f_Ha)
    Z_N2_17 = np.interp(N2,N2_cal_data_2017[:,0],N2_cal_data_2017[:,1],left=np.nan, right=np.inf)
    x_N2_17 = Z_N2_17 - Z_sun
    return Z_N2_17

def compute_Z_N2_Curti20(line_df):
    #  N2=[N II]λ6584/Ha, NII6583
    line_IDs  = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_NII6583     = line_df[line_IDs.index('NII6583')].data
    e_NII6583     = line_df[line_IDs.index('NII6583_ERR')].data
    f_Ha      = line_df[line_IDs.index('HALPHA')].data
    e_Ha      = line_df[line_IDs.index('HALPHA_ERR')].data
    # compute line ratios
    N2 = np.log10(f_NII6583/f_Ha)
    Z_N2_20 = np.interp(N2,N2_cal_data_2020[:,0],N2_cal_data_2020[:,1],left=np.nan, right=np.inf)
    x_N2_20 = Z_N2_20 - Z_sun
    return Z_N2_20

def compute_Z_O3N2_Curti17(line_df):
    # O3N2 = [O III]λ5007/Hb / [N II]λ6584/Ha; OIII5007; NII6583
    line_IDs  = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_OIII5007     = line_df[line_IDs.index('OIII5007')].data
    e_OIII5007     = line_df[line_IDs.index('OIII5007_ERR')].data
    f_NII6583     = line_df[line_IDs.index('NII6583')].data
    e_NII6583     = line_df[line_IDs.index('NII6583_ERR')].data
    f_Hb      = line_df[line_IDs.index('HBETA')].data
    e_Hb      = line_df[line_IDs.index('HBETA_ERR')].data
    f_Ha      = line_df[line_IDs.index('HALPHA')].data
    e_Ha      = line_df[line_IDs.index('HALPHA_ERR')].data
    # compute line ratios
    O3     = np.log10(f_OIII5007/f_Hb) 
    N2     = np.log10(f_NII6583/f_Ha)
    O3N2   = O3 - N2
    Z_O3N2_17 = np.interp(O3N2,O3N2_cal_data_2017[:,0],O3N2_cal_data_2017[:,1],left=np.nan, right=np.inf)
    x_O3N2_17 = Z_O3N2_17 - Z_sun
    return

def compute_Z_O3N2_Curti20(line_df):
    # O3N2 = [O III]λ5007/Hb / [N II]λ6584/Ha; OIII5007; NII6583
    line_IDs  = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_OIII5007     = line_df[line_IDs.index('OIII5007')].data
    e_OIII5007     = line_df[line_IDs.index('OIII5007_ERR')].data
    f_NII6583     = line_df[line_IDs.index('NII6583')].data
    e_NII6583     = line_df[line_IDs.index('NII6583_ERR')].data
    f_Hb      = line_df[line_IDs.index('HBETA')].data
    e_Hb      = line_df[line_IDs.index('HBETA_ERR')].data
    f_Ha      = line_df[line_IDs.index('HALPHA')].data
    e_Ha      = line_df[line_IDs.index('HALPHA_ERR')].data
    # compute line ratios
    O3     = np.log10(f_OIII5007/f_Hb) 
    N2     = np.log10(f_NII6583/f_Ha)
    O3N2   = O3 - N2
    Z_O3N2_20 = np.interp(O3N2,O3N2_cal_data_2020[:,0],O3N2_cal_data_2020[:,1],left=np.nan, right=np.inf)
    x_O3N2_20 = Z_O3N2_20 - Z_sun
    return Z_O3N2_20

def compute_Z_O3O2_Curti17(line_df):
    # [O III]λ5007/[O II]λ3727,29
    line_IDs  = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_OIII5007     = line_df[line_IDs.index('OIII5007')].data
    e_OIII5007     = line_df[line_IDs.index('OIII5007_ERR')].data
    f_OII3726     = line_df[line_IDs.index('OII3726')].data
    e_OII3726     = line_df[line_IDs.index('OII3726_ERR')].data
    f_OII3729     = line_df[line_IDs.index('OII3729')].data
    e_OII3729     = line_df[line_IDs.index('OII3729_ERR')].data
    # compute line ratios
    O3O2 = np.log10(f_OIII5007/(f_OII3726+f_OII3729))
    Z_O3O2_17 = np.interp(O3O2,O3O2_cal_data_2017[:,0],O3O2_cal_data_2017[:,1],left=np.nan, right=np.inf)
    x_O3O2_17 = Z_O3O2_17 - Z_sun
    return Z_O3O2_17

def compute_Z_O3O2_Curti20(line_df):
    # [O III]λ5007/[O II]λ3727,29
    line_IDs  = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_OIII5007     = line_df[line_IDs.index('OIII5007')].data
    e_OIII5007     = line_df[line_IDs.index('OIII5007_ERR')].data
    f_OII3726     = line_df[line_IDs.index('OII3726')].data
    e_OII3726     = line_df[line_IDs.index('OII3726_ERR')].data
    f_OII3729     = line_df[line_IDs.index('OII3729')].data
    e_OII3729     = line_df[line_IDs.index('OII3729_ERR')].data
    # compute line ratios
    O3O2 = np.log10(f_OIII5007/(f_OII3726+f_OII3729))
    Z_O3O2_20 = np.interp(O3O2,O3O2_cal_data_2020[:,0],O3O2_cal_data_2020[:,1],left=np.nan, right=np.inf)
    x_O3O2_20 = Z_O3O2_20 - Z_sun
    return Z_O3O2_20

def compute_Z_R23_Curti17(line_df):
    # R23=([O II]λ3727,29 + [O III]λ4959,5007)/Hb
    line_IDs  = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_OIII5007     = line_df[line_IDs.index('OIII5007')].data
    e_OIII5007     = line_df[line_IDs.index('OIII5007_ERR')].data
    f_OII3726     = line_df[line_IDs.index('OII3726')].data
    e_OII3726     = line_df[line_IDs.index('OII3726_ERR')].data
    f_OII3729     = line_df[line_IDs.index('OII3729')].data
    e_OII3729     = line_df[line_IDs.index('OII3729_ERR')].data
    f_Hb      = line_df[line_IDs.index('HBETA')].data
    e_Hb      = line_df[line_IDs.index('HBETA_ERR')].data
    # compute line ratios
    R23 = np.log10((f_OIII5007+f_OII3726+f_OII3729)/f_Hb)
    Z_R23_17 = np.interp(R23,R23_cal_data_2017[:,0],R23_cal_data_2017[:,1],left=np.nan, right=np.inf)
    x_R23_17 = Z_R23_17 - Z_sun
    return Z_R23_17

def compute_Z_R23_Curti20(line_df):
    # R23=([O II]λ3727,29 + [O III]λ4959,5007)/Hb
    line_IDs  = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_OIII5007     = line_df[line_IDs.index('OIII5007')].data
    e_OIII5007     = line_df[line_IDs.index('OIII5007_ERR')].data
    f_OII3726     = line_df[line_IDs.index('OII3726')].data
    e_OII3726     = line_df[line_IDs.index('OII3726_ERR')].data
    f_OII3729     = line_df[line_IDs.index('OII3729')].data
    e_OII3729     = line_df[line_IDs.index('OII3729_ERR')].data
    f_Hb      = line_df[line_IDs.index('HBETA')].data
    e_Hb      = line_df[line_IDs.index('HBETA_ERR')].data
    # compute line ratios
    R23 = np.log10((f_OIII5007+f_OII3726+f_OII3729)/f_Hb)
    Z_R23_20 = np.interp(R23,R23_cal_data_2020[:,0],R23_cal_data_2020[:,1],left=np.nan, right=np.inf)
    x_R23_20 = Z_R23_20 - Z_sun
    return Z_R23_20

def compute_Z_O3S2_Curti20(line_df):
    # O3S2=[O III]λ5007/Hb / [S II]λ6717,31/Ha; SII6716;  SII6731
    line_IDs  = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_OIII5007     = line_df[line_IDs.index('OIII5007')].data
    e_OIII5007     = line_df[line_IDs.index('OIII5007_ERR')].data
    f_SII6716      = line_df[line_IDs.index('SII6716')].data
    e_SII6716      = line_df[line_IDs.index('SII6716_ERR')].data
    f_SII6731      = line_df[line_IDs.index('SII6731')].data
    e_SII6731      = line_df[line_IDs.index('SII6731_ERR')].data
    f_Ha      = line_df[line_IDs.index('HALPHA')].data
    e_Ha      = line_df[line_IDs.index('HALPHA_ERR')].data
    f_Hb      = line_df[line_IDs.index('HBETA')].data
    e_Hb      = line_df[line_IDs.index('HBETA_ERR')].data
    # compute line ratios
    O3S2 = np.log10((f_OIII5007/f_Hb)/((f_SII6716+f_SII6731)/f_Ha))
    Z_O3S2_20 = np.interp(O3S2,O3S2_cal_data_2020[:,0],O3S2_cal_data_2020[:,1],left=np.nan, right=np.inf)
    x_O3S2_20 = Z_O3S2_20 - Z_sun
    return Z_O3S2_20

def compute_Z_RS32_Curti20(line_df):
    # RS32 = [O III]λ5007/Hb + [S II]λ6717,31/Ha
    line_IDs  = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_OIII5007     = line_df[line_IDs.index('OIII5007')].data
    e_OIII5007     = line_df[line_IDs.index('OIII5007_ERR')].data
    f_SII6716      = line_df[line_IDs.index('SII6716')].data
    e_SII6716      = line_df[line_IDs.index('SII6716_ERR')].data
    f_SII6731      = line_df[line_IDs.index('SII6731')].data
    e_SII6731      = line_df[line_IDs.index('SII6731_ERR')].data
    f_Ha      = line_df[line_IDs.index('HALPHA')].data
    e_Ha      = line_df[line_IDs.index('HALPHA_ERR')].data
    f_Hb      = line_df[line_IDs.index('HBETA')].data
    e_Hb      = line_df[line_IDs.index('HBETA_ERR')].data
    # compute line ratios
    RS32 = np.log10(f_OIII5007/f_Hb + (f_SII6716+f_SII6731)/f_Ha)
    Z_RS32_20 = np.interp(RS32,RS32_cal_data_2020[:,0],RS32_cal_data_2020[:,1],left=np.nan, right=np.inf)
    x_RS32_20 = Z_RS32_20 - Z_sun
    return Z_RS32_20

def compute_Z_S2_Curti20(line_df):
    # S2 = [S II]λ6717,31/Ha
    line_IDs  = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_SII6716      = line_df[line_IDs.index('SII6716')].data
    e_SII6716      = line_df[line_IDs.index('SII6716_ERR')].data
    f_SII6731      = line_df[line_IDs.index('SII6731')].data
    e_SII6731      = line_df[line_IDs.index('SII6731_ERR')].data
    f_Ha      = line_df[line_IDs.index('HALPHA')].data
    e_Ha      = line_df[line_IDs.index('HALPHA_ERR')].data
    # compute line ratios
    S2 = np.log10((f_SII6716+f_SII6731)/f_Ha)
    Z_S2_20 = np.interp(S2,S2_cal_data_2020[:,0],S2_cal_data_2020[:,1],left=np.nan, right=np.inf)
    x_S2_20 = Z_S2_20 - Z_sun
    return Z_S2_20


