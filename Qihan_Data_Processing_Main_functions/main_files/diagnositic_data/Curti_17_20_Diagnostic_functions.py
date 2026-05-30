# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 17:56:06 2023

@author: Qihan Zou
"""

import numpy as np 
Z_sun = 8.69

if __name__=='__main__':
    '''
    If the name is something_2020 then is from Curti+20:
    https://ui.adsabs.harvard.edu/abs/2020MNRAS.491..944C/abstract
    The following calculation follow the formula 
    log(R)=sum c_n * x^n where x = z - 8.69
    R is a given diagnositic
    x is the oxygen aboundance normialized to the solar value 12+log(O/H)_sun = 8.69
    
    If the name is something_2017 then is from Curti+17:
    https://ui.adsabs.harvard.edu/abs/2017MNRAS.465.1384C/abstract
    The following calculation follow the formula 
    log(R)=sum c_n * x^n where x = z - 8.69
    R is a given diagnositic
    x is the oxygen aboundance normialized to the solar value 12+log(O/H)_sun = 8.69
    
    a = '1234'
    a[::-1]
    '4321'
    
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

    
    
    '''
    x = np.arange(7.6, 8.9, 0.01) - Z_sun 
    '''
    This part is for Curti+20:
    '''
    # R2_2020
    R2_2020    = 0.435 -1.362*x -5.655*(x**2) -4.851*(x**3) -0.478*(x**4) +0.736*(x**5)
    R2_2020_cal_data = np.vstack((R2_2020[::-1], x[::-1]+Z_sun)).T
    np.savetxt('C:/Users/qihan/Desktop/diagnositic_data/Curti20_R2.txt',R2_2020_cal_data)
    # R3_2020
    R3_2020    = -0.277 -3.549*x -4.593*(x**2) -0.981*(x**3)
    R3_2020_cal_data = np.vstack((R3_2020[::-1], x[::-1]+Z_sun)).T
    np.savetxt('C:/Users/qihan/Desktop/diagnositic_data/Curti20_R3.txt',R3_2020_cal_data)
    # O3O2_2020
    O3O2_2020  = -0.691 -2.944*x -1.308*(x**2)
    O3O2_2020_cal_data = np.vstack((O3O2_2020[::-1], x[::-1]+Z_sun)).T
    np.savetxt('C:/Users/qihan/Desktop/diagnositic_data/Curti20_O3O2.txt',O3O2_2020_cal_data)
    # R23_2020
    R23_2020   = 0.527 -1.569*x - 1.652*(x**2) - 0.421*(x**3)
    R23_2020_cal_data = np.vstack((R23_2020[::-1], x[::-1]+Z_sun)).T
    np.savetxt('C:/Users/qihan/Desktop/diagnositic_data/Curti20_R23.txt',R23_2020_cal_data)
    # N2_2020
    N2_2020    = -0.489 + 1.513*x - 2.554**(x**2) - 5.293*(x**3) -2.867**(x**4)
    N2_2020_cal_data = np.vstack((N2_2020[::-1], x[::-1]+Z_sun)).T
    np.savetxt('C:/Users/qihan/Desktop/diagnositic_data/Curti20_N2.txt',N2_2020_cal_data)
    # O3N2_2020
    O3N2_2020  = 0.281 -4.765*x -2.268*(x**2)
    O3N2_2020_cal_data = np.vstack((O3N2_2020[::-1], x[::-1]+Z_sun)).T
    np.savetxt('C:/Users/qihan/Desktop/diagnositic_data/Curti20_O3N2.txt',O3N2_2020_cal_data)
    # S2_2020
    S2_2020    = -0.442 -0.360*x -6.271*(x**2)-8.339*(x**3)-3.559*(x**4)
    S2_2020_cal_data = np.vstack((S2_2020[::-1], x[::-1]+Z_sun)).T
    np.savetxt('C:/Users/qihan/Desktop/diagnositic_data/Curti20_S2.txt',S2_2020_cal_data)
    # RS32_2020
    RS32_2020  = -0.054 -2.546*x -1.970*(x**2)+0.082*(x**3)+0.222**(x**4)
    RS32_2020_cal_data = np.vstack((RS32_2020[::-1], x[::-1]+Z_sun)).T
    np.savetxt('C:/Users/qihan/Desktop/diagnositic_data/Curti20_RS32.txt',RS32_2020_cal_data)
    # O3S2_2020
    O3S2_2020  = 0.191 -4.292*x -2.538*(x**2) +0.053*(x**3) +0.332*(x**4)
    O3S2_2020_cal_data = np.vstack((O3S2_2020[::-1], x[::-1]+Z_sun)).T
    np.savetxt('C:/Users/qihan/Desktop/diagnositic_data/Curti20_O3S2.txt',O3S2_2020_cal_data)
    '''
    Next part is for Curti+17:
    '''
    # R2_2017
    R2_2017 = 0.418-0.961*x-3.505*(x**2)-1.949*(x**3)
    R2_2017_cal_data = np.vstack((R2_2017[::-1], x[::-1]+Z_sun)).T
    np.savetxt('C:/Users/qihan/Desktop/diagnositic_data/Curti17_R2.txt',R2_2017_cal_data)
    # R3_2017
    R3_2017 = -0.277-3.549*x-3.593*(x**2)-0.981*(x**3)
    R3_2017_cal_data = np.vstack((R3_2017[::-1], x[::-1]+Z_sun)).T
    np.savetxt('C:/Users/qihan/Desktop/diagnositic_data/Curti17_R3.txt',R3_2017_cal_data)
    # O3O2_2017
    O3O2_2017 = -0.691-2.944*x-1.308*(x**2)
    O3O2_2017_cal_data = np.vstack((O3O2_2017[::-1], x[::-1]+Z_sun)).T
    np.savetxt('C:/Users/qihan/Desktop/diagnositic_data/Curti17_O3O2.txt',O3O2_2017_cal_data)
    # R23_2017
    R23_2017 = 0.527-1.569*x-1.652*(x**2)-0.421*(x**3)
    R23_2017_cal_data = np.vstack((R23_2017[::-1], x[::-1]+Z_sun)).T
    np.savetxt('C:/Users/qihan/Desktop/diagnositic_data/Curti17_R23.txt',R23_2017_cal_data)
    # N2_2017
    N2_2017 = -0.489 +1.513*x - 2.554*(x**2)-5.293*(x**3)-2.867*(x**4)
    N2_2017_cal_data = np.vstack((N2_2017[::-1], x[::-1]+Z_sun)).T
    np.savetxt('C:/Users/qihan/Desktop/diagnositic_data/Curti17_N2.txt',N2_2017_cal_data)
    # O3N2_2017
    O3N2_2017 = 0.281 -4.765*x -2.268*(x**2)
    O3N2_2017_cal_data = np.vstack((O3N2_2017[::-1], x[::-1]+Z_sun)).T
    np.savetxt('C:/Users/qihan/Desktop/diagnositic_data/Curti17_O3N2.txt',O3N2_2017_cal_data)
    
    
    
    
    
    
    