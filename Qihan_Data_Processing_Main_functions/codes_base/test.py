# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 09:37:58 2024

@author: qihan
"""

from MUSE_data_processing import *
MUSE_plot_galaxy(gal_name = 'N4321', file_path = r"C:\Users\qihan\Desktop\q\MUSE_pkl\N4321_ver1_copt.pkl", 
                     DIG_CUT="N2_BPT", Z_train = "Z_N2S2Ha", 
                     galaxy_info_path = "C:/Users/qihan/Desktop/q/galaxydata.xlsx")




MUSE_plot_galaxy(gal_name = 'N4321', file_path = r"C:\Users\qihan\Desktop\q\MUSE_pkl\N4321_ver2_copt.pkl", 
                     DIG_CUT="N2_BPT", Z_train = "Z_N2S2Ha", 
                     galaxy_info_path = "C:/Users/qihan/Desktop/q/galaxydata.xlsx")



gal_df = fits.open(r"C:\Users\qihan\Desktop\q\MUSE\NGC4321_MAPS_copt_1.16asec.fits")
gal_df.info()


'''
Filename: C:\Users\qihan\Desktop\q\MUSE\NGC4321_MAPS_copt_1.16asec.fits
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU     104   ()      
  1  BIN_ID        1 ImageHDU        24   (1231, 999)   float64   
  2  V_STARS       1 ImageHDU        26   (1231, 999)   float64   
  3  FORM_ERR_V_STARS    1 ImageHDU        26   (1231, 999)   float64   
  4  SIGMA_STARS    1 ImageHDU        26   (1231, 999)   float64   
  5  FORM_ERR_SIGMA_STARS    1 ImageHDU        26   (1231, 999)   float64   
  6  HB4861_FLUX    1 ImageHDU        26   (1231, 999)   float64   
  7  HB4861_FLUX_ERR    1 ImageHDU        26   (1231, 999)   float64   
  8  HB4861_VEL    1 ImageHDU        26   (1231, 999)   float64   
  9  HB4861_VEL_ERR    1 ImageHDU        26   (1231, 999)   float64   
 10  HB4861_SIGMA    1 ImageHDU        26   (1231, 999)   float64   
 11  HB4861_SIGMA_ERR    1 ImageHDU        26   (1231, 999)   float64   
 12  OIII4958_FLUX    1 ImageHDU        26   (1231, 999)   float64   
 13  OIII4958_FLUX_ERR    1 ImageHDU        26   (1231, 999)   float64   
 14  OIII4958_VEL    1 ImageHDU        26   (1231, 999)   float64   
 15  OIII4958_VEL_ERR    1 ImageHDU        26   (1231, 999)   float64   
 16  OIII4958_SIGMA    1 ImageHDU        26   (1231, 999)   float64   
 17  OIII4958_SIGMA_ERR    1 ImageHDU        26   (1231, 999)   float64   
 18  OIII5006_FLUX    1 ImageHDU        26   (1231, 999)   float64   
 19  OIII5006_FLUX_ERR    1 ImageHDU        26   (1231, 999)   float64   
 20  OIII5006_VEL    1 ImageHDU        26   (1231, 999)   float64   
 21  OIII5006_VEL_ERR    1 ImageHDU        26   (1231, 999)   float64   
 22  OIII5006_SIGMA    1 ImageHDU        26   (1231, 999)   float64   
 23  OIII5006_SIGMA_ERR    1 ImageHDU        26   (1231, 999)   float64   
 24  NII6548_FLUX    1 ImageHDU        26   (1231, 999)   float64   
 25  NII6548_FLUX_ERR    1 ImageHDU        26   (1231, 999)   float64   
 26  NII6548_VEL    1 ImageHDU        26   (1231, 999)   float64   
 27  NII6548_VEL_ERR    1 ImageHDU        26   (1231, 999)   float64   
 28  NII6548_SIGMA    1 ImageHDU        26   (1231, 999)   float64   
 29  NII6548_SIGMA_ERR    1 ImageHDU        26   (1231, 999)   float64   
 30  HA6562_FLUX    1 ImageHDU        26   (1231, 999)   float64   
 31  HA6562_FLUX_ERR    1 ImageHDU        26   (1231, 999)   float64   
 32  HA6562_VEL    1 ImageHDU        26   (1231, 999)   float64   
 33  HA6562_VEL_ERR    1 ImageHDU        26   (1231, 999)   float64   
 34  HA6562_SIGMA    1 ImageHDU        26   (1231, 999)   float64   
 35  HA6562_SIGMA_ERR    1 ImageHDU        26   (1231, 999)   float64   
 36  NII6583_FLUX    1 ImageHDU        26   (1231, 999)   float64   
 37  NII6583_FLUX_ERR    1 ImageHDU        26   (1231, 999)   float64   
 38  NII6583_VEL    1 ImageHDU        26   (1231, 999)   float64   
 39  NII6583_VEL_ERR    1 ImageHDU        26   (1231, 999)   float64   
 40  NII6583_SIGMA    1 ImageHDU        26   (1231, 999)   float64   
 41  NII6583_SIGMA_ERR    1 ImageHDU        26   (1231, 999)   float64   
 42  SII6716_FLUX    1 ImageHDU        26   (1231, 999)   float64   
 43  SII6716_FLUX_ERR    1 ImageHDU        26   (1231, 999)   float64   
 44  SII6716_VEL    1 ImageHDU        26   (1231, 999)   float64   
 45  SII6716_VEL_ERR    1 ImageHDU        26   (1231, 999)   float64   
 46  SII6716_SIGMA    1 ImageHDU        26   (1231, 999)   float64   
 47  SII6716_SIGMA_ERR    1 ImageHDU        26   (1231, 999)   float64   
 48  SII6730_FLUX    1 ImageHDU        26   (1231, 999)   float64   
 49  SII6730_FLUX_ERR    1 ImageHDU        26   (1231, 999)   float64   
 50  SII6730_VEL    1 ImageHDU        26   (1231, 999)   float64   
 51  SII6730_VEL_ERR    1 ImageHDU        26   (1231, 999)   float64   
 52  SII6730_SIGMA    1 ImageHDU        26   (1231, 999)   float64   
 53  SII6730_SIGMA_ERR    1 ImageHDU        26   (1231, 999)   float64   

'''