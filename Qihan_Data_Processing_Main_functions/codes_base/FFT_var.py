# -*- coding: utf-8 -*-
"""
Created on Sat Jun 29 20:47:43 2024

@author: qihan
"""

import math
import numpy as np
from scipy.fft import fft2, ifft2





def FFT_variograms(x1,x2):
    n = x1.shape[0]
    p = x1.shape[1]
    
    n2 = x2.shape[0]
    #p2 = x2.shape[1]
    
    nrows = 2*n-1
    ncols = 2*p-1
    nr2 = math.ceil(nrows/8)*8
    nc2 = math.ceil(ncols/8)*8
    
    # Replace missing values with 0
    x1_id = ~np.isnan(x1)
    x1[~x1_id] = 0
    

    # Compute Fourier transforms
    fx1 = np.fft.fft2(x1, nr2, nc2)
    fx1_x1 = np.fft.fft2(x1 * x1, nr2, nc2)
    fx1_id = np.fft.fft2(x1_id, nr2, nc2)
    nh11 = np.round(np.real(ifft2(np.conj(fx1_id) * fx1_id)))

    # Handle x2 if defined
    if n2 > 0:
       x2_id = ~np.isnan(x2)
       x2[~x2_id] = 0
       fx2 = fft2(x2, nr2, nc2)
       fx2_x2 = fft2(x2 * x2, nr2, nc2)
       fx2_id = fft2(x2_id, nr2, nc2)
       nh22 = np.round(np.real(ifft2(np.conj(fx2_id) * fx2_id)))
       fx12_id = fft2(x1_id * x2_id, nr2, nc2)
       nh12 = np.round(np.real(ifft2(np.conj(fx12_id) * fx12_id)))

   # Compute structural functions
    if n2 > 0:
       gh22 = np.real(ifft2(np.conj(fx2_id) * fx2_x2 + np.conj(fx2_x2) * fx2_id - 2 * np.conj(fx2) * fx2)) / np.maximum(nh22, 1) / 2
       t1 = fft2(x1 * x2_id, nr2, nc2)
       t2 = fft2(x2 * x1_id, nr2, nc2)
       t12 = fft2(x1 * x2, nr2, nc2)
       gh12 = np.real(ifft2(np.conj(fx12_id) * t12 + np.conj(t12) * (fx12_id) - np.conj(t1) * t2 - t1 * np.conj(t2))) / np.maximum(nh12, 1) / 2
       
    gh11 = np.real(ifft2(np.conj(fx1_id) * fx1_x1 + np.conj(fx1_x1) * fx1_id - 2 * np.conj(fx1) * fx1)) / np.maximum(nh11, 1) / 2

    # Reduce matrices to required size
    nh11 = np.concatenate((nh11[0:n, 0:p], nh11[0:n, nc2-p+1:nc2], nh11[nr2-n+1:nr2, 0:p], nh11[nr2-n+1:nr2, nc2-p+1:nc2]), axis=1)
    gh11 = np.concatenate((gh11[0:n, 0:p], gh11[0:n, nc2-p+1:nc2], gh11[nr2-n+1:nr2, 0:p], gh11[nr2-n+1:nr2, nc2-p+1:nc2]), axis=1)
    if n2 > 0:
       nh22 = np.concatenate((nh22[0:n, 0:p], nh22[0:n, nc2-p+1:nc2], nh22[nr2-n+1:nr2, 0:p], nh22[nr2-n+1:nr2, nc2-p+1:nc2]), axis=1)
       gh22 = np.concatenate((gh22[0:n, 0:p], gh22[0:n, nc2-p+1:nc2], gh22[nr2-n+1:nr2, 0:p], gh22[nr2-n+1:nr2, nc2-p+1:nc2]), axis=1)
       nh12 = np.concatenate((nh12[0:n, 0:p], nh12[0:n, nc2-p+1:nc2], nh12[nr2-n+1:nr2, 0:p], nh12[nr2-n+1:nr2, nc2-p+1:nc2]), axis=1)
       gh12 = np.concatenate((gh12[0:n, 0:p], gh12[0:n, nc2-p+1:nc2], gh12[nr2-n+1:nr2, 0:p], gh12[nr2-n+1:nr2, nc2-p+1:nc2]), axis=1)

    # Shift matrices to center the 0 lag
    nh11 = np.fft.fftshift(nh11)
    gh11 = np.fft.fftshift(gh11)
    if n2 > 0:
       nh22 = np.fft.fftshift(nh22)
       gh22 = np.fft.fftshift(gh22)
       nh12 = np.fft.fftshift(nh12)
       gh12 = np.fft.fftshift(gh12)
    return gh11,nh11,gh12,nh12,gh22,nh22


#gh11,nh11,gh12,nh12,gh22,nh22 = FFT_variograms(m1,m2)
#m1 = matrix([[3,6,5], [7,2,2], [4,NaN,0]])
#m2 = matrix([[10,NaN,5], [NaN,8,7], [5,9,11]])


def bartrand_fast_semivariogram(datain):
	'''
	An implementation of the fast semivariogram algorithm by Jonah Bartrand,
	from the Colorado School of Mines.
	
	Based on Marcotte's algorithm, but this one might be correct.
	
	Source:
	https://wtools.readthedocs.io/en/latest/_modules/wtools/geostats/raster.html#raster_to_struct_grid
	
	Parameters
	----------
	datain: np.array
		A random field.
		
	rtol: float
		Relative tolerance. Default is 1e-10.
		
	Returns
	-------
	ND_svg: np.array
		The N-dimensional semivariogram computed using this algorithm.
		
	n_pairs: np.array
		
	'''
	# Padding
	data_dims = datain.shape
	nDim = len(data_dims)
	out_dims = [2*d-1 for d in data_dims]
	
	# Build indicator matrix
	data_loc_ind = ~np.isnan(datain) * ~np.isinf(datain)
	# In data matrix, replace missing values by 0:
	datain[~data_loc_ind] = 0
	
	# Construct the F.T'd matrices that we will need:
	fD	= np.fft.fftn(datain, s=out_dims)
	#fDD = np.fft.fftn(datain*datain, s=out_dims)				#These are unused -- but are used in Marcotte's algorithm?
	fI	= np.fft.fftn(data_loc_ind, s=out_dims)
	#fID = np.fft.fftn(datain*data_loc_ind, s=out_dims)
	
	# Compute number of pairs at all lags
	n_pairs = np.real(np.fft.ifftn(np.abs(fI)**2)).astype(int)
	
	# Compute covariance with something that looks like Marcotte's algorithm,
	# but isn't. 
	cov = np.real(	np.fft.ifftn(np.abs(fD)**2) /
					np.fft.ifftn(np.abs(fI)**2) -
					np.fft.ifftn(np.conj(fD)*fI) *
					np.fft.ifftn(np.conj(fI)*fD) /
					(np.fft.ifftn(np.abs(fI)**2))**2
				)
	# Using this, compute the semivariogram
	ND_svg = np.max(cov)-cov
	
	# Reduce matrix to required size 
	# and shift so that the 0 lag appears at the center of each matrix
	unpad_ind	= [[int(d/2),int(3*d/2)] for d in data_dims]
	unpad_list	= [np.arange(*l) for l in unpad_ind]
	unpad_coord = np.meshgrid(*unpad_list, indexing='ij')

	ND_svg	= np.fft.fftshift(ND_svg) [tuple(unpad_coord)]
	n_pairs = np.fft.fftshift(n_pairs)[tuple(unpad_coord)]
	
	return ND_svg, n_pairs

def raster_to_struct_grid(datain, imeas='covar', rtol=1e-10):
    """Create an auto-variogram or auto-covariance map from 1D or 2D rasters.
    This computes auto-variogram or auto-covariance maps from
    1D or 2D rasters. This function computes variograms/covariances in the
    frequency domain via the Fast Fourier Transform (``np.fftn``).

    Note:
        For viewing the results, please use the ``plot_struct_grid`` method
        from the ``plots`` module.

    Note:
        Missing values, flagged as ``np.nan``, are allowed.

    Args:
        datain (np.ndarray): input arrray with raster in GeoEas format
        imeas (str): key indicating which structural measure to compute:
            ``'var'`` for semi-variogram or ``'covar'`` for covariogram.
        gridspecs (list(GridSpec)): array with grid specifications using
            ``GridSpec`` objects
        rtol (float): the tolerance. Default is 1e-10

    Return:
        tuple(np.ndarray, np.ndarray):
            output array with variogram or covariogram map, depending
            on variogram choice, with size: in 1D: ( 2*nxOutHalf+1 ) or in 2D:
            ( 2*nxOutHalf+1 x 2*nxOutHalf+1 ).

            output array with number of pairs available in each lag,
            of same size as outStruct

    References:
        Originally implemented in MATLAB by:
            Phaedon Kyriakidis,
            Department of Geography,
            University of California Santa Barbara,
            May 2005

        Reimplemented into Python by:
            Jonah Bartrand,
            Department of Geophysics,
            Colorado School of Mines,
            October 2018

        Algorith based on:
            Marcotte, D. (1996): Fast Variogram Computation with FFT,
            Computers & Geosciences, 22(10), 1175-1186.
    """
    # Check imeas
    itypes = ['covar', 'var']
    if isinstance(imeas, int) and imeas < 2 and imeas > -1:
        imeas = itypes[imeas]
    if imeas not in itypes:
        raise RuntimeError("imeas argument must be one of 'covar' for covariogram or 'var' for semi-variance. Not {}".format(imeas))

    data_dims = datain.shape
    nDim = len(data_dims)

    ## Get appropriate dimensions
    # find the closest multiple of 8 to obtain a good compromise between
    # speed (a power of 2) and memory required
    out_dims = [2*d-1 for d in data_dims]#[int(np.ceil((2*d-1)/8)*8) for d in data_dims]

    ## Form an indicator  matrix:
    # 0's for all data values, 1's for missing values
    missing_data_ind = np.isnan(datain);
    data_loc_ind = np.logical_not(missing_data_ind)
    # In data matrix, replace missing values by 0;
    datain[missing_data_ind] = 0  # missing replaced by 0

    ## FFT of datain
    fD = np.fft.fftn(datain, s=out_dims)

    ## FFT of datain*datain
    fDD = np.fft.fftn(datain*datain, s=out_dims)

    ## FFT of the indicator matrix
    fI = np.fft.fftn(data_loc_ind, s=out_dims)

    ## FFT of datain*indicator
    fID = np.fft.fftn(datain*data_loc_ind, s=out_dims)

    ## Compute number of pairs at all lags
    outNpairs = np.real(np.fft.ifftn(np.abs(fI)**2)).astype(int)
    #Edit remove single formating for matlab v6
    #outNpairs = single(outNpairs);

    cov = np.real(  np.fft.ifftn(np.abs(fD)**2) /
                    np.fft.ifftn(np.abs(fI)**2) -
                    np.fft.ifftn(np.conj(fD)*fI) *
                    np.fft.ifftn(np.conj(fI)*fD) /
                    (np.fft.ifftn(np.abs(fI)**2))**2
                )

    if imeas == 'var':
        outStruct = np.max(cov)-cov
    else:
        outStruct = cov

    ## Reduce matrix to required size and shift,
    # so that the 0 lag appears at the center of each matrix

    unpad_ind = [[int(d/2),int(3*d/2)] for d in data_dims]
    unpad_list = [np.arange(*l) for l in unpad_ind]
    unpad_coord = np.meshgrid(*unpad_list, indexing='ij')

    outStruct=np.fft.fftshift(outStruct)[tuple(unpad_coord)]
    outNpairs=np.fft.fftshift(outNpairs)[tuple(unpad_coord)]

    indzeros = outNpairs<(np.max(outNpairs)*rtol)
    outStruct[indzeros] = np.nan

    return outStruct, outNpairs



A = np.zeros((10, 10))

# Assign specific values to elements of the matrix
A[0, 0] = 10
A[0, 4] = 14
A[0, 9] = 8
A[1, 2] = 12
A[3, 0] = 14
A[5, 2] = 18
A[6, 0] = 17
A[9, 0] = 14

print(A)

