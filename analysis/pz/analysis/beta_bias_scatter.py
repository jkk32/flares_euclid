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


a = np.arange(-3.,2.,0.5)
beta_bins = np.array([a[:-1], a[1:]]).T

print(beta_bins)

zr = [5,12]
z_binw = 0.5
z_bins = np.arange(zr[0]+z_binw/2., zr[1], z_binw)





fig = plt.figure(figsize = (6, 1.5))

panel_width = 0.3


ax = {}
ax['bias'] = fig.add_axes((0.1, 0.1, panel_width, 0.8))
ax['scatter'] = fig.add_axes((0.1 + panel_width + (1-0.2-panel_width*2), 0.1, panel_width, 0.8))

scale_ax = {}
scale_ax['bias'] = fig.add_axes((0.08, 0.1, 0.015, 0.8))
scale_ax['scatter'] = fig.add_axes((0.905, 0.1, 0.015, 0.8))


for j,beta_bin in enumerate(beta_bins):
    fig.text(0.5, 0.8*(len(beta_bins)-j)/(len(beta_bins)), fplt.ml(fr'log_{{10}}(\beta)\in[{",".join(map(str, beta_bin))})'), fontsize = 8, ha = 'center', va = 'center', color = 'k') # add filter labels


for fig_type in ['bias','scatter']:

    ax[fig_type].set_xticks([])
    ax[fig_type].set_yticks([])

    if fig_type == 'bias':
        vmin = -5
        vmax = 5
        cmap = 'RdBu_r'
    if fig_type == 'scatter':
        vmin = 0
        vmax = 10
        cmap = 'Reds'

    scale_ax[fig_type].imshow(np.expand_dims(np.linspace(vmax,vmin,100), axis=0).T, cmap = cmap, vmin = vmin, vmax = vmax, aspect = 'auto', extent=[0,1,vmin,vmax])

    if fig_type == 'bias': scale_ax[fig_type].set_ylabel(fplt.ml('bias\ P_{50}'))
    if fig_type == 'scatter': scale_ax[fig_type].set_ylabel(fplt.ml('scatter\ P_{84}-P_{16}'))
    scale_ax[fig_type].set_xticks([])
    if fig_type == 'scatter':
        scale_ax[fig_type].yaxis.tick_right()
        scale_ax[fig_type].yaxis.set_label_position("right")





    met = np.zeros((len(beta_bins), len(z_bins)))

    for j,beta_bin in enumerate(beta_bins):

        sel = (beta>beta_bin[0]) & (beta<beta_bin[1])

        dz = z_EAZY[sel] - z_input[sel]
        if fig_type == 'bias': met[j,:] = np.array([np.median(dz[(np.fabs(z_input[sel]-(z+0.5))<0.5)]) for z in z_bins])
        if fig_type == 'scatter': met[j,:] = np.array([np.percentile(dz[(np.fabs(z_input[sel]-(z+0.5))<0.5)], 84.) - np.percentile(dz[(np.fabs(z_input[sel]-(z+0.5))<0.5)], 16.) for z in z_bins])

    ax[fig_type].imshow(met, cmap = cmap, vmin = vmin, vmax = vmax, aspect = 'auto')

    for i,z in enumerate(np.arange(zr[0], zr[1]+1, 1.0)):
        ax[fig_type].text(i/z_binw-0.5, -0.6, f'${z:.0f}$', fontsize = 5, ha = 'center', va = 'bottom')
        ax[fig_type].text(i/z_binw-0.5, len(beta_bins)+0.1, f'${z:.0f}$', fontsize = 5, ha = 'center', va = 'bottom')


fig.savefig(f'figures/beta.pdf')
fig.clf()
