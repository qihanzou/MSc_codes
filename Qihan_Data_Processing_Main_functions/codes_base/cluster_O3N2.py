# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 23:08:27 2024

@author: qihan
"""

from TYPHOON_wrangling import *
import pandas as pd
import numpy as np
import logging
from astropy.io import fits
from SFR_calculation import *
from surface_gas_density_calculation import *
from stellar_mass_calculation import *
from Four_relations_calculation_sun2023 import *
import pickle
from AST2 import *
from Internal_Fun import *
from BFuns import *
from sklearn.datasets import make_blobs
from sklearn.preprocessing import StandardScaler
import numpy as np
from sklearn import metrics
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from dec_ra_to_xy import *
import scipy.spatial as spatial
import itertools
# curve-fit() function imported from scipy
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt


meta = meta_getter('N5236')

'''
For O3N2
'''
data_galaxy = pd.read_pickle('C:/Users/qihan/Desktop/q/N5236_25_1_2024.pkl')
index_above_0_Z = data_galaxy["Z_O3N2"] > 0 # ignore NA Z_N2S2Ha Z_RS32 Z_O3N2
data_galaxy = data_galaxy[index_above_0_Z]
#index_below_1_BPT = data_galaxy["N2_BPT"] < 1 # BPT cut for Hii regions
#data_galaxy = data_galaxy[index_below_1_BPT] # result: 7193 rows of data
index_CHII = data_galaxy["S2_DIG"] > 0.9 # BPT cut for Hii regions
data_galaxy = data_galaxy[index_CHII] # result: 7193 rows of data
RA_galaxy_list = data_galaxy['RA']
DEC_galaxy_list = data_galaxy['DEC']
proj_dist_galaxy_list = data_galaxy['proj_dist']
coor_full_galaxy = RA_DEC_to_xy(RA_galaxy_list, DEC_galaxy_list, meta)
X = np.transpose(coor_full_galaxy)

plt.scatter(X[:, 0], X[:, 1], c = 'blue', marker='.', s=1)
ax = plt.gca()
ax.set_xlim([-15, 15])
ax.set_ylim([-12, 12])
plt.title("Hii regions")
plt.grid()
plt.show()




r = 1
critical_number_of_points = 150
empty_set = list()
index_set = list()


for ii in range(X.shape[0]):
    target_point =  X[ii].reshape((1,2)) 
    X_tree = spatial.cKDTree(X)
    positions = X_tree.query_ball_point(target_point, r)
    positions = np.concatenate(positions).astype(None)
    number_of_points = len(positions) 
    if number_of_points >= critical_number_of_points:
       empty_set.append(target_point)
       index_set.append(ii)
     
a = array(empty_set)
result = a[:, 0, :]

plt.scatter(X[:, 0], X[:, 1], c = 'blue', marker='.', s=1)
ax = plt.gca()
ax.set_xlim([-15, 15])
ax.set_ylim([-12, 12])
plt.scatter(result[:,0], result[:, 1], c = 'red', marker='.', s=1)
plt.title("Hii regions and spiral arms")
plt.grid()
plt.show()        


db = DBSCAN(eps=1, min_samples=10).fit(result)
labels = db.labels_

# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
n_noise_ = list(labels).count(-1)

print("Estimated number of clusters: %d" % n_clusters_)
print("Estimated number of noise points: %d" % n_noise_)


unique_labels = set(labels)
core_samples_mask = np.zeros_like(labels, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True

colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = [0, 0, 0, 1]

    class_member_mask = labels == k

    xy = result[class_member_mask & core_samples_mask]
    plt.plot(
        xy[:, 0],
        xy[:, 1],
        "o",
        markerfacecolor=tuple(col),
        markeredgecolor="k",
        markersize=5,
    )

    xy = result[class_member_mask & ~core_samples_mask]
    plt.plot(
        xy[:, 0],
        xy[:, 1],
        "*",
        markerfacecolor=tuple(col),
        markeredgecolor="k",
        markersize=5,
    )

plt.title(f"Estimated number of clusters: {n_clusters_}")
plt.show()

i, = np.where(labels == 0)
j, = np.where(labels == 1)
cluster1 = result[i]
cluster2 = result[j]

plt.scatter(X[:, 0], X[:, 1], c = 'blue', marker='.', s=1)
ax = plt.gca()
ax.set_xlim([-15, 15])
ax.set_ylim([-12, 12])
plt.scatter(cluster1[:,0], cluster1[:, 1], c = 'orange', marker='.', s=1)
plt.scatter(cluster2[:,0], cluster2[:, 1], c = 'yellow', marker='.', s=1)
plt.title("Hii regions and spiral arms in different colors")
plt.grid()
plt.show()

index_set1 = list()
for jj in range(len(labels)):
    if labels[jj] == 0:
       index_set1.append(1)
    if labels[jj] == 1:
       index_set1.append(2)
       
final_indicator = list(itertools.repeat(0, X.shape[0]))
for kk in range(len(index_set)):
    index = index_set[kk]
    value = index_set1[kk]
    final_indicator[index] = value
    
final_indicator = array(final_indicator)
#np.savetxt("N5236_O3N2_indicator_012.csv", final_indicator, delimiter=",", header="I_arms")
