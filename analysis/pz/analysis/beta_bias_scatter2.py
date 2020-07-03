import os
import sys
import numpy as np

import scipy.stats

import h5py

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import FLARE.plt as fplt
import FLARE.filters


input = h5py.File('../../../data2/beta.hdf5', 'r')
pz = h5py.File('../../../data2/beta_pz.hdf5', 'r')

SN = pz['obs/Euclid.NISP.H/flux'][:]/pz['obs/Euclid.NISP.H/error'][:]
beta = input['beta'][:]

s = (SN>5)

z_input = input['z'][s]
z_EAZY = pz['EAZY/z_m1'][s]
dz = z_EAZY - z_input
beta = beta[s]


br = [-3, 2]
beta_binw = 0.5
beta_bins = np.arange(br[0]+beta_binw/2., br[1], beta_binw) # really bin centres

zr = [5, 12]
z_binw = 0.5
z_bins = np.arange(zr[0]+z_binw/2., zr[1], z_binw) # really bin centres

ar = len(beta_bins)/len(z_bins)

print(ar)



# # not the same bins

# fig = plt.figure(figsize = (3,3))
#
# ax = fig.add_axes((0.15, 0.15, 0.8, 0.8))
# met = np.zeros((len(beta_bins), len(z_bins)))
#
# N, xedges, yedges = np.histogram2d(z_input, beta, bins = [z_bins, beta_bins])
#
# ax.imshow(N, cmap = cmap, vmin = 0, vmax = 100, aspect = 'auto', extent = [*zr, *br])
#
# fig.savefig(f'figures/beta_N.pdf')
# fig.clf()




# --- bias

vmin = -5
vmax = 5
cmap = mpl.cm.RdBu_r
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

left = 0.15
bottom = 0.15
width = 0.8
height = 0.65
cbar_space = 0.15 # space needed by cbar, not

fig = plt.figure(figsize = (3, 3*(ar/width)))

ax = fig.add_axes((left, bottom, width, height))
met = np.zeros((len(beta_bins), len(z_bins)))

for j,beta_bin in enumerate(beta_bins):
    sel = np.fabs(beta-beta_bin)<beta_binw/2
    dz = z_EAZY[sel] - z_input[sel]
    met[j,:] = np.array([np.median(dz[(np.fabs(z_input[sel]-z)<z_binw/2)]) for z in z_bins])

ax.imshow(met, vmin=vmin, vmax=vmax, cmap = cmap, aspect = 'auto', extent = [*zr, *br])


cbar_ax = fig.add_axes((left, bottom+height, width, 0.025))
cbar = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, orientation='horizontal')
cbar.set_label(r'$\rm P_{50}(\Delta z)$', size=8)
cbar.ax.tick_params(labelsize=5)
cbar_ax.xaxis.set_ticks_position('top')
cbar_ax.xaxis.set_label_position('top')

ax.set_xlabel(r'$\rm z$')
ax.set_ylabel(r'$\rm \beta$')

fig.savefig(f'figures/beta_bias.pdf')
fig.clf()




# --- scatter

vmin = 0
vmax = 10
cmap = mpl.cm.Reds
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

left = 0.15
bottom = 0.15
width = 0.8
height = 0.65
cbar_space = 0.15 # space needed by cbar, not

fig = plt.figure(figsize = (3, 3*(ar/width)))

ax = fig.add_axes((left, bottom, width, height))
met = np.zeros((len(beta_bins), len(z_bins)))

for j,beta_bin in enumerate(beta_bins):
    sel = np.fabs(beta-beta_bin)<beta_binw/2
    dz = z_EAZY[sel] - z_input[sel]
    met[j,:] = np.array([np.percentile(dz[(np.fabs(z_input[sel]-(z+0.5))<0.5)], 84.) - np.percentile(dz[(np.fabs(z_input[sel]-(z+0.5))<0.5)], 16.) for z in z_bins])


ax.imshow(met, vmin=vmin, vmax=vmax, cmap = cmap, aspect = 'auto', extent = [*zr, *br])


cbar_ax = fig.add_axes((left, bottom+height, width, 0.025))
cbar = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, orientation='horizontal')
cbar.set_label(r'$\rm P_{84}(\Delta z)-P_{16}(\Delta z)$', size=8)
cbar.ax.tick_params(labelsize=5)
cbar_ax.xaxis.set_ticks_position('top')
cbar_ax.xaxis.set_label_position('top')

ax.set_xlabel(r'$\rm z$')
ax.set_ylabel(r'$\rm \beta$')

fig.savefig(f'figures/beta_scatter.pdf')
fig.clf()
