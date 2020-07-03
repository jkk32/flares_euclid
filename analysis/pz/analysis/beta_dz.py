

import os
import sys
import numpy as np

import scipy.stats

import h5py

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import FLARE.plt as fplt

input = h5py.File('../../../data2/beta.hdf5', 'r')
pz = h5py.File('../../../data2/beta_pz.hdf5', 'r')



SN = pz['obs/Euclid.NISP.H/flux'][:]/pz['obs/Euclid.NISP.H/error'][:]
beta = input['beta'][:]

s = (SN>5)&(beta<0)&(SN<20)

print(len(SN[s]))


fig, ax = fplt.simple_fig()

z_input = input['z'][s]
z_EAZY = pz['EAZY/z_m1'][s]
dz = z_EAZY - z_input

# --- plot lines of different tau

zr = [5,12]

z_binw = 0.05
bins = np.arange(zr[0]+z_binw/2., zr[1], z_binw)


N = [len(dz[(np.fabs(z_input-b)<z_binw/2.)]) for b in bins]
print(N)

med = [np.median(dz[(np.fabs(z_input-b)<z_binw/2.)]) for b in bins]
P16 = [np.percentile(dz[(np.fabs(z_input-b)<z_binw/2.)], 16.) for b in bins]
P84 = [np.percentile(dz[(np.fabs(z_input-b)<z_binw/2.)], 84.) for b in bins]

c='k'

ax.plot(bins, med, c=c, label=r'$\rm median\ (P_{50})$')
ax.fill_between(bins, P16, P84, alpha = 0.2, color=c, label=r'$\rm [P_{16},P_{86}]$')


# ax.legend(fontsize = 7)
ax.plot(zr,[0,0], alpha = 0.2, c='k', lw=1) # 1:1 line
ax.set_xlim(zr)
ax.set_ylim([-5,5])
ax.set_xlabel(fplt.fancy('z'))
ax.set_ylabel(fplt.fancy('\Delta z = z_{EAZY} - z'))
ax.legend(fontsize=7)
fig.savefig(f'figures/beta_dz.pdf')
fig.clf()
