'''
Z_diags.py

Running it once generates Z diagnostics
Once that has happened it can be imported/used by other functions.

Created by Benjamin Metha
'''

import numpy as np 
# from uncertainties import unumpy as unp #will write this another day maybe
Z_sun = 8.69

# Import memoisations
if __name__ != '__main__':
    O3N2_cal_data      = np.loadtxt('C:/Users/qihan/Desktop/q/Curti17_O3N2.txt')
    O3S2_cal_data_2017 = np.loadtxt('C:/Users/qihan/Desktop/q/Curti17_O3S2.txt')
    RS32_cal_data_2020 = np.loadtxt('C:/Users/qihan/Desktop/q/Curti20_RS32.txt')

def compute_Z_N2S2Ha_Dop16(line_df):
    '''
    Given a set of deredenned emission line maps+error, 
    compute metallicity maps+error, using the
    N2S2Ha diagnostic of Dopita+2016:
    
    
    Parameters
    ----------
    line_df: hdu list
        A big guy containing all the different emission line data reduced
        from TYPHOON data cubes
    https://ui.adsabs.harvard.edu/abs/2016Ap&SS.361...61D
    Returns
    -------
    Z: array
        Metallicity using this diagnostic
        
    e_Z: array
        Error in metallicity using this diagnostic, 
        computed via linear error propagation.
    '''
    # Unpack the wanted lines: f_NII, f_SII6717, f_SII6731, f_Ha,
    # and their errors.
    line_IDs = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_NII = line_df[line_IDs.index('NII6583')].data
    e_NII = line_df[line_IDs.index('NII6583_ERR')].data
    f_SII6716 = line_df[line_IDs.index('SII6716')].data
    e_SII6716 = line_df[line_IDs.index('SII6716_ERR')].data
    f_SII6731 = line_df[line_IDs.index('SII6731')].data
    e_SII6731 = line_df[line_IDs.index('SII6731_ERR')].data
    f_Ha = line_df[line_IDs.index('HALPHA')].data
    e_Ha = line_df[line_IDs.index('HALPHA_ERR')].data
    # compute relevant line ratios
    N2S2  = np.log10(f_NII/(f_SII6716+f_SII6731) )
    N2Ha  = np.log10(f_NII/f_Ha)
    N2S2Ha = N2S2 + 0.264*N2Ha
    Z_N2S2Ha_low = 8.77 + N2S2Ha
    Z_N2S2Ha_upper_correction = 0.45 * (N2S2Ha+0.3)**5
    Z_N2S2Ha = Z_N2S2Ha_low + Z_N2S2Ha_upper_correction*(Z_N2S2Ha_low > 9.05)
    # and errors
    dZ_dN2S2Ha = 1 + 2.25*(N2S2Ha+0.3)**4 *(Z_N2S2Ha_low > 9.05)
    dratio_dN2 = 1.264/(np.log(10)*f_NII)
    dratio_dS2 = 1.0/(np.log(10)*(f_SII6716+f_SII6731))
    dratio_dHa = 0.264/(np.log(10)*f_Ha)
    e_Z2 = dZ_dN2S2Ha**2 * ((dratio_dN2*e_NII)**2 + (dratio_dS2*e_SII6716)**2 \
                            + (dratio_dS2*e_SII6731)**2 + (dratio_dHa*e_Ha)**2)
    e_Z  = np.sqrt(e_Z2)
    return Z_N2S2Ha, e_Z
     
def compute_Z_O3N2_Curti17(line_df, kumari_correction = False):
    '''
    Given a set of deredenned emission line maps+error, 
    compute metallicity maps+error, using the
    N2S2Ha diagnostic of Curti+2017:
    https://ui.adsabs.harvard.edu/abs/2017MNRAS.465.1384C/abstract
    
    Optionally, add a DIG correction as devised by Kumari+19:
    https://ui.adsabs.harvard.edu/abs/2019MNRAS.485..367K
    
    Parameters
    ----------
    line_df: hdu list
        A big guy containing all the different emission line data reduced
        from TYPHOON data cubes
        
    kumari_correction: bool or str
        if False; no correction
        if 'N2': add the published correction for DIG as selected using an
            N2-BPT diagram
        if 'S2': add the published correction for DIG as selected using an
            S2-BPT diagram
    
    Returns
    -------
    Z: array
        Metallicity using this diagnostic
        
    e_Z: array
        Error in metallicity using this diagnostic, 
        computed via linear error propagation.
    '''
    # Unpack the wanted lines: f_NII, f_OIII, f_Ha, f_Hb,
    # and their errors.
    line_IDs = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_NII  = line_df[line_IDs.index('NII6583')].data
    e_NII  = line_df[line_IDs.index('NII6583_ERR')].data
    f_OIII = line_df[line_IDs.index('OIII5007')].data
    e_OIII = line_df[line_IDs.index('OIII5007_ERR')].data
    f_Ha   = line_df[line_IDs.index('HALPHA')].data
    e_Ha   = line_df[line_IDs.index('HALPHA_ERR')].data
    f_Hb   = line_df[line_IDs.index('HBETA')].data
    e_Hb   = line_df[line_IDs.index('HBETA_ERR')].data
    # define line ratios
    O3     = np.log10(f_OIII/f_Hb) 
    N2     = np.log10(f_NII/f_Ha)
    O3N2   = O3 - N2 
    # compute Z_O3N2 by interpolating the inverse function
    Z_o3n2 = np.interp(O3N2,O3N2_cal_data[:,0],O3N2_cal_data[:,1],left=np.nan, right=np.inf)
    x_o3n2 = Z_o3n2 - Z_sun
    # compute error function
    dZ_dratio= 1.0/(4.765+4.536*x_o3n2)
    if kumari_correction==False:
        e_Z2 = (dZ_dratio/np.log(10))**2 * ((e_NII/f_NII)**2 + (e_OIII/f_OIII)**2 + (e_Hb/f_Hb)**2+ (e_Ha/f_Ha)**2)
    if kumari_correction=='S2':
        # Then add a DIG correction based on the O3 line;
        # and our square error is a little different.
        Z_o3n2 = Z_o3n2 + 0.156*O3
        e_Z2 = (dZ_dratio/np.log(10))**2 * ((e_NII/f_NII)**2 + (1.156*e_OIII/f_OIII)**2 + (1.156*e_Hb/f_Hb)**2+ (e_Ha/f_Ha)**2)
    if kumari_correction=='N2':
        # same idea
        Z_o3n2 = Z_o3n2 + 0.033+0.127*O3
        e_Z2 = (dZ_dratio/np.log(10))**2 * ((e_NII/f_NII)**2 + (1.127*e_OIII/f_OIII)**2 + (1.127*e_Hb/f_Hb)**2+ (e_Ha/f_Ha)**2)
    e_Z = np.sqrt(e_Z2)
    return Z_o3n2, e_Z



def compute_Z_O3S2_Curti17(line_df, kumari_correction=False):
    '''
    Given a set of deredenned emission line maps+error, 
    compute metallicity maps+error, using the
    N2S2Ha diagnostic of Curti+2017:
    https://ui.adsabs.harvard.edu/abs/2017MNRAS.465.1384C/abstract
    
    Optionally, add a DIG correction as devised by Kumari+19:
    https://ui.adsabs.harvard.edu/abs/2019MNRAS.485..367K
    
    Parameters
    ----------
    line_df: hdu list
        A big guy containing all the different emission line data reduced
        from TYPHOON data cubes
        
    kumari_correction: bool or str
        if False; no correction
        if 'N2': add the published correction for DIG as selected using an
            N2-BPT diagram
        if 'S2': add the published correction for DIG as selected using an
            S2-BPT diagram
    
    Returns
    -------
    Z: array
        Metallicity using this diagnostic
        
    e_Z: array
        Error in metallicity using this diagnostic, 
        computed via linear error propagation.
    '''
    # Unpack the wanted lines: f_SII6717, f_SII6731, f_OIII, f_Ha, f_Hb,
    # and their errors.
    line_IDs = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_SII6716 = line_df[line_IDs.index('SII6716')].data
    e_SII6716 = line_df[line_IDs.index('SII6716_ERR')].data
    f_SII6731 = line_df[line_IDs.index('SII6731')].data
    e_SII6731 = line_df[line_IDs.index('SII6731_ERR')].data
    f_OIII = line_df[line_IDs.index('OIII5007')].data
    e_OIII = line_df[line_IDs.index('OIII5007_ERR')].data
    f_Ha   = line_df[line_IDs.index('HALPHA')].data
    e_Ha   = line_df[line_IDs.index('HALPHA_ERR')].data
    f_Hb   = line_df[line_IDs.index('HBETA')].data
    e_Hb   = line_df[line_IDs.index('HBETA_ERR')].data
    # define handy dandy ratios
    O3     = np.log10(f_OIII/f_Hb) 
    combo_flux = f_OIII/f_Hb + (f_SII6716+f_SII6731)/f_Ha
    combo_S2 = f_SII6716+f_SII6731
    O3S2   = np.log10(combo_flux)
    Z_o3s2 = np.interp(O3S2,O3S2_cal_data_2017[:,0],O3S2_cal_data_2017[:,1],left=np.nan, right=np.inf)
    x_o3s2 = Z_o3s2 - Z_sun
    # compute error function
    dZ_dratio= 1.0/(-2.223 - 2.146*x_o3s2 + 1.599* x_o3s2**2 )
    if kumari_correction==False:
        e_Z2 = (dZ_dratio/(np.log(10)*combo_flux))**2 * ( (e_OIII/f_Hb)**2 + (f_OIII*e_Hb/(f_Hb**2))**2+(combo_S2*e_Ha/(f_Ha**2))**2 +(e_SII6716/f_Hb)**2+ (e_SII6731/f_Hb)**2)
    if kumari_correction=='S2':
        # Then add a DIG correction based on the O3 line;
        # and the associated error.
        Z_o3s2 = Z_o3s2 + 0.075 + 0.309*O3 + 0.208*(O3**2)
        # first account for Ha,S2 contributions:
        e_Z2 = (dZ_dratio/(np.log(10)*combo_flux))**2 * ( (combo_S2*e_Ha/(f_Ha**2))**2 +(e_SII6716/f_Hb)**2+ (e_SII6731/f_Hb)**2)
        dZ_dO3 = 0.309+0.416*O3
        e_Z2_from_OIII = (dZ_dratio/(np.log(10)*combo_flux*f_Hb) +dZ_dO3/(np.log(10)*f_OIII)  )**2 * e_OIII**2
        e_Z2_from_Hb = (dZ_dratio*f_OIII/(np.log(10)*combo_flux*f_Hb*f_Hb) +dZ_dO3/(np.log(10)*f_Hb)  )**2 * e_Hb**2
        e_Z2 = e_Z2 + e_Z2_from_OIII + e_Z2_from_Hb
    if kumari_correction=='N2':
        # Then add a DIG correction based on the O3 line;
        # and the associated error
        Z_o3s2 = Z_o3s2 + 0.113 + 0.229*O3
        e_Z2 = (dZ_dratio/(np.log(10)*combo_flux))**2 * ( (combo_S2*e_Ha/(f_Ha**2))**2 +(e_SII6716/f_Hb)**2+ (e_SII6731/f_Hb)**2)
        dZ_dO3 = 0.229
        e_Z2_from_OIII = (dZ_dratio/(np.log(10)*combo_flux*f_Hb) +dZ_dO3/(np.log(10)*f_OIII)  )**2 * e_OIII**2
        e_Z2_from_Hb = (dZ_dratio*f_OIII/(np.log(10)*combo_flux*f_Hb*f_Hb) +dZ_dO3/(np.log(10)*f_Hb)  )**2 * e_Hb**2
        e_Z2 = e_Z2 + e_Z2_from_OIII + e_Z2_from_Hb
    e_Z = np.sqrt(e_Z2)
    return Z_o3s2, e_Z

def compute_Z_RS32_Curti20(line_df):
    '''
    Given a set of deredenned emission line maps+error, 
    compute metallicity maps+error, using the
    N2S2Ha diagnostic of Curti+20:
    https://ui.adsabs.harvard.edu/abs/2020MNRAS.491..944C/abstract
    
    Parameters
    ----------
    line_df: hdu list
        A big guy containing all the different emission line data reduced
        from TYPHOON data cubes
    
    Returns
    -------
    Z: array
        Metallicity using this diagnostic
        
    e_Z: array
        Error in metallicity using this diagnostic, 
        computed via linear error propagation.
    '''
    line_IDs = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_SII6716 = line_df[line_IDs.index('SII6716')].data
    e_SII6716 = line_df[line_IDs.index('SII6716_ERR')].data
    f_SII6731 = line_df[line_IDs.index('SII6731')].data
    e_SII6731 = line_df[line_IDs.index('SII6731_ERR')].data
    f_OIII = line_df[line_IDs.index('OIII5007')].data
    e_OIII = line_df[line_IDs.index('OIII5007_ERR')].data
    f_Ha   = line_df[line_IDs.index('HALPHA')].data
    e_Ha   = line_df[line_IDs.index('HALPHA_ERR')].data
    f_Hb   = line_df[line_IDs.index('HBETA')].data
    e_Hb   = line_df[line_IDs.index('HBETA_ERR')].data
    # define handy dandy ratios
    O3     = np.log10(f_OIII/f_Hb) 
    combo_flux = f_OIII/f_Hb + (f_SII6716+f_SII6731)/f_Ha
    combo_S2 = f_SII6716+f_SII6731
    O3S2   = np.log10(combo_flux)
    # Interpolate metallicities
    Z_o3s2 = np.interp(O3S2,RS32_cal_data_2020[:,0],RS32_cal_data_2020[:,1],left=np.nan, right=np.inf)
    x_o3s2 = Z_o3s2 - Z_sun
    # compute error function
    # 1/Derivative (d/dx) of:  -0.054 -2.546*x -1.970*(x**2) + 0.082* (x**3) + 0.222*(x**4)
    dZ_dratio= 1.0/(-2.546 - 3.940*x_o3s2 + 0.246*(x_o3s2**2) + 0.888*(x_o3s2**3) )
    # Compute error using linear error propagation
    e_Z2 = (dZ_dratio/(np.log(10)*combo_flux))**2 * ( (e_OIII/f_Hb)**2 + (f_OIII*e_Hb/(f_Hb**2))**2+(combo_S2*e_Ha/(f_Ha**2))**2 +(e_SII6716/f_Hb)**2+ (e_SII6731/f_Hb)**2)
    e_Z = np.sqrt(e_Z2)
    return Z_o3s2, e_Z


def compute_Z_N2O2_Dop13(line_df):
    '''
    Given a set of deredenned emission line maps+error, 
    compute metallicity maps+error, using the
    N2O2 diagnostic of Dopita+2013:
    https://ui.adsabs.harvard.edu/abs/2013ApJS..208...10D/abstract
    
    Parameters
    ----------
    line_df: hdu list
        A big guy containing all the different emission line data reduced
        from TYPHOON data cubes

    Returns
    -------
    Z: array
        Metallicity using this diagnostic
        
    e_Z: array
        Error in metallicity using this diagnostic, 
        computed via linear error propagation.
    '''
    # Unpack the wanted lines: f_NII, f_OII3727, f_OII3729,
    # and their errors.
    line_IDs = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_OII3727 = line_df[line_IDs.index('OII3726')].data
    e_OII3727 = line_df[line_IDs.index('OII3726_ERR')].data
    f_OII3729 = line_df[line_IDs.index('OII3729')].data
    e_OII3729 = line_df[line_IDs.index('OII3729_ERR')].data
    f_NII  = line_df[line_IDs.index('NII6583')].data
    e_NII  = line_df[line_IDs.index('NII6583_ERR')].data
    # add the hecking, uh, parts.
    R = np.log10(f_NII/(f_OII3727+f_OII3729))
    Q = 1.5402+1.26602*R+0.167977*(R**2)
    Z = np.log10(Q) + 8.93
    # find error using linear error propagation
    dZ_dR  = (1.26602 + 0.335954*R)/(np.log(10)*Q)
    dR_dN2 = 1/(np.log(10)*f_NII)
    dR_dO2 = 1/(np.log(10)*(f_OII3727+f_OII3729))
    e_Z2 = (dZ_dR*dR_dN2*e_NII)**2 + (dZ_dR*dR_dO2)**2 *(e_OII3727**2 + e_OII3729**2)
    e_Z = np.sqrt(e_Z2)
    return Z, e_Z

def compute_O3(line_df):
    '''
    Just give me log_10(OIII/Ha) and error; no diagnostics.
    '''
    line_IDs = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_OIII   = line_df[line_IDs.index('OIII5007')].data
    e_OIII   = line_df[line_IDs.index('OIII5007_ERR')].data
    f_Hb     = line_df[line_IDs.index('HBETA')].data
    e_Hb     = line_df[line_IDs.index('HBETA_ERR')].data
    O3       = np.log10(f_OIII/f_Hb) 
    e_O3_sq  = np.log10(np.e)**2 * ( (e_Hb/f_Hb)**2 + (e_OIII/f_OIII)**2 )
    e_O3     = np.sqrt(e_O3_sq)
    return O3, e_O3

if __name__=='__main__':
    # Run this to generate/save memoisations.
    x = np.arange(7.6, 8.9, 0.01) - Z_sun 
    # for accuracies of 0.01 in Zworld
    # Range over which metallicities are valid taken from Curti+2020.
    O3N2 = 0.281 -4.765*x -2.268*(x**2)
    O3N2_cal_data = np.vstack((O3N2[::-1], x[::-1]+Z_sun)).T
    np.savetxt('C:/Users/qihan/Desktop/q/Curti17_O3N2.txt',O3N2_cal_data)
    # Aaaaand for O3S2
    O3S2_2017 = -0.046 -2.223*x -1.073*(x**2) + 0.533* (x**3) 
    O3S2_cal_data_2017 = np.vstack((O3S2_2017[::-1], x[::-1]+Z_sun)).T
    np.savetxt('C:/Users/qihan/Desktop/q/Curti17_O3S2.txt',O3S2_cal_data_2017)
    O3S2_2020 = -0.054 -2.546*x -1.970*(x**2) + 0.082* (x**3) + 0.222*(x**4)
    O3S2_cal_data_2020 = np.vstack((O3S2_2020[::-1], x[::-1]+Z_sun)).T
    np.savetxt('C:/Users/qihan/Desktop/q/Curti20_O3S2.txt',O3S2_cal_data_2020)