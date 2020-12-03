import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.cm as cm
mpl.use('Agg')
import matplotlib.pyplot as plt

import flares

import FLARE.photom as photom
import FLARE.plt as fplt


# --- define colour scale

cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin=9., vmax=11.)


# --- used to define lower x-limit

plot_flux_limit = 2. # log10(nJy)
print(f'log10(flux_limit/nJy): {plot_flux_limit}')




# --- read in FLARS data

fl = flares.flares('../../data/flares.hdf5', sim_type='FLARES')
halo = fl.halos

# -- define data of interest
# -- CREATE A LOOP OVER PARAMETERS OF INTEREST


df = pd.read_csv('../../../flares/weight_files/weights_grid.txt')
weights = np.array(df['weights'])


fig, axes = plt.subplots(len(fl.tags), 3, figsize=(3*2, len(fl.tags)*2), sharey = True, sharex = True)

plt.subplots_adjust(left=0.1, bottom=0.05, right=0.9, top=0.95, wspace=0.0, hspace=0.0)


filters = ['FUV', 'NUV', 'V']

Mstar = fl.load_dataset('Mstar_30', arr_type=f'Galaxy') # M_sol

H = fl.load_dataset('H', arr_type=f'Galaxy/BPASS_2.2.1/Chabrier300/Flux/DustModelI/Euclid/NISP/')


phot_types = ['Pure_Stellar','Intrinsic', 'DustModelI']
phot_labels = ['Stellar', 'Stellar + Nebular', 'Stellar + Nebular + Dust']

for j, phot_type in enumerate(phot_types):

    L = {}
    for filter in filters:

        #inst = '/'.join(filter.split('/')[:2])

        f = filter #.split('/')[-1]
        L[filter] = fl.load_dataset(f, arr_type=f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/{phot_type}/') # ergs/s/Hz


    for i, tag in enumerate(fl.tags): # loop over snapshots (redshifts)

        z = fl.zeds[i]

        ws = np.array([])
        mstar = np.array([])
        h = np.array([])

        l = {}
        for filter in filters:
            l[filter]=np.array([])

        for ii in range(len(halo)):

            ws = np.append(ws, np.ones(np.shape(L['V'][halo[ii]][tag]))*weights[ii])
            mstar = np.append(mstar, Mstar[halo[ii]][tag] + 10.)
            h = np.append(h, H[halo[ii]][tag])

            for filter in filters:
                l[filter] = np.append(l[filter], L[filter][halo[ii]][tag])



        c1 = -2.5*np.log10(l['FUV']/l['NUV'])
        c2 = -2.5*np.log10(l['NUV']/l['V'])

        s = h>10**plot_flux_limit # dust attenuated luminosity
        # ws = ws[s]

        print(len(s[s==True]), np.mean(c1[s]), np.mean(c2[s]))

        ax = axes[i, j]

        ax.hexbin(c1[s], c2[s], gridsize=50, extent=[-0.5, 0.8, -0.5, 2.4], norm=mpl.cm.colors.LogNorm(vmin=1, vmax=400),
                  cmap=mpl.cm.Greys, linewidths=0.1)

        #ax.scatter(c1[s], c2[s], c = norm(mstar[s]), s=2, alpha=0.25)

        if j==0: axes[i, 0].text(0.05, 0.9, rf'$\rm z={z:.0f}$', ha = 'left', va = 'baseline', transform=axes[i, 0].transAxes)

    axes[0, j].text(0.5, 1.05, phot_labels[j], ha = 'center', va = 'baseline', transform=axes[0, j].transAxes)



fig.text(0.01, 0.5, 'FUV-NUV', ha = 'left', va = 'center', rotation = 'vertical')
fig.text(0.5,0.01, 'NUV-V', ha = 'center', va = 'bottom')


# fig.add_subplot(111, frameon=False)
# plt.xlabel('H-ch1')
# plt.ylabel('ch1-ch2')
# plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
#


fig.savefig(f'figures/paper/colours_rest_FUV_NUV_V.pdf')
fig.clf()
