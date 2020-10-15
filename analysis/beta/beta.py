import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.cm as cm
mpl.use('Agg')
import matplotlib.pyplot as plt

import flares

import FLARE.photom as photom
import FLARE.plt as fplt

import FLARE

cosmo = FLARE.default_cosmo()

def beta(l1, l2, lum1, lum2):
    # Wilkins+15
    #return (2.5 * np.log10(l1 / l2)) ** (-1) * (m1 - m2) - 2
    return (np.log10(l1 / l2)) ** (-1) * np.log10(lum1/lum2) - 2

# --- define colour scale

cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin=9., vmax=11.)


# --- used to define lower x-limit

plot_flux_limit = 1. # log10(nJy)
print(f'log10(flux_limit/nJy): {plot_flux_limit}')




# --- read in FLARS data

fl = flares.flares('../../data/flares.hdf5', sim_type='FLARES')
halo = fl.halos

# -- define data of interest
# -- CREATE A LOOP OVER PARAMETERS OF INTEREST


df = pd.read_csv('../../../flares/weight_files/weights_grid.txt')
weights = np.array(df['weights'])


###
###
###
# CHANGE THE REST TO REFLECT THE BETA MEASUREMENTS

# https://www.aanda.org/articles/aa/pdf/2016/12/aa27782-15.pdf

#l_nuv = 2315.7 #\AA
#l_fuv = 1538.6 #\AA

l_nuv = 2500. #\AA
l_fuv = 1500. #\AA


fig, axes = plt.subplots(len(fl.tags), 3, figsize=(3*2, len(fl.tags)*2), sharey = True, sharex = True)

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.0, hspace=0.0)


filters = ['FUV', 'NUV']

Mstar = fl.load_dataset('Mstar_30', arr_type=f'Galaxy') # M_sol

H = fl.load_dataset('H', arr_type=f'Galaxy/BPASS_2.2.1/Chabrier300/Flux/DustModelI/Euclid/NISP/')


phot_types = ['Pure_Stellar','Intrinsic', 'DustModelI']
phot_labels = ['Stellar', 'Stellar + Nebular', 'Stellar + Nebular + Dust']

betas = {}
zetas = []

beta_stellar = {}
beta_nebular = {}
beta_summary = {}

cmap_p = mpl.cm.plasma


beta_all = {}

for j, phot_type in enumerate(phot_types):

    L = {}
    for filter in filters:

        #inst = '/'.join(filter.split('/')[:2])

        f = filter #.split('/')[-1]
        L[filter] = fl.load_dataset(f, arr_type=f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/{phot_type}/') # ergs/s/Hz


    beta_all[phot_type] = {}


    for i, tag in enumerate(fl.tags): # loop over snapshots (redshifts)

        z = fl.zeds[i]

        ws = np.array([])
        mstar = np.array([])
        h = np.array([])

        l = {}
        for filter in filters:
            l[filter]=np.array([])

        for ii in range(len(halo)):

            ws = np.append(ws, np.ones(np.shape(L['NUV'][halo[ii]][tag]))*weights[ii])
            mstar = np.append(mstar, Mstar[halo[ii]][tag] + 10.)
            h = np.append(h, H[halo[ii]][tag])

            for filter in filters:
                l[filter] = np.append(l[filter], L[filter][halo[ii]][tag])



        c1 = np.log10(h)
        c2 = beta(l_fuv, l_nuv, l['FUV'], l['NUV'])


        beta_all[phot_type][str(z)] = {'beta': c2, 'log10H': c1}


        s = h>10**plot_flux_limit # dust attenuated luminosity
        # ws = ws[s]

        if phot_type == 'DustModelI':
            betas[str(z)] = {'beta': c2[s], 'log10H': c1[s]}
            zetas.append(z)

        print(np.mean(photom.flux_to_m(photom.lum_to_flux(l['NUV'], cosmo=cosmo, z=z))), np.mean(photom.flux_to_m(photom.lum_to_flux(l['FUV'], cosmo=cosmo, z=z))))
        print(len(s[s==True]), np.mean(c1[s]), np.mean(c2[s]))

        ax = axes[i, j]

        cmap = mpl.cm.colors.Colormap('Greys', N=256)

        vmin = np.log10(1.)
        vmax = np.log10(1000.)

        ax.hexbin(c1[s], c2[s], gridsize=50, extent=[1., 3.7, -2.9, -0.6], norm=mpl.cm.colors.LogNorm(vmin=1, vmax=400), cmap=mpl.cm.Greys, linewidths=0.1)
        #ax.scatter(c1[s], c2[s], c = norm(mstar[s]), s=2, alpha=0.25)

        bin_edges = np.linspace(.9, 3.7, 30)
        dh_bin = bin_edges[1] - bin_edges[0]
        beta_spine = []
        beta_84 = []
        beta_16 = []
        bin_centres = []
        for bin in bin_edges[:-1]:
            try:
                sd = abs(c1[s] - (bin + dh_bin/2)) < dh_bin
                beta_spine.append(np.percentile(c2[s][sd], 50))
                beta_84.append(np.percentile(c2[s][sd], 84))
                beta_16.append(np.percentile(c2[s][sd], 16))
                bin_centres.append(bin + dh_bin/2)
            except:
                continue

        if phot_type == 'DustModelI':
            beta_summary[str(z)] = {'beta': beta_spine, 'log10H': bin_centres}

        if phot_type == 'Pure_Stellar':
            beta_stellar[str(z)] = {'beta': beta_spine, 'log10H': bin_centres}

        if phot_type == 'Intrinsic':
            beta_nebular[str(z)] = {'beta': beta_spine, 'log10H': bin_centres}

        ax.axhline(-2, color='k', ls='-', alpha=0.3)

        ax.axhline(np.median(c2[s]), color='r', ls='--', alpha=0.3)

        ax.fill_between(bin_centres, beta_16, beta_84, color=cmap_p(1.5 * (float(z) - 5) / 10), alpha=0.3)
        ax.plot(bin_centres, beta_spine, c=cmap_p(1.5 * (float(z) - 5) / 10), alpha=0.3)

        if j==0: axes[i, 0].text(0.05, 0.9, rf'$\rm z={z:.0f}$', ha = 'left', va = 'baseline', transform=axes[i, 0].transAxes)

    axes[0, j].text(0.5, 1.05, phot_labels[j], ha = 'center', va = 'baseline', transform=axes[0, j].transAxes)

plt.xlim(1)

fig.text(0.01, 0.5, r'$\rm \beta$', ha = 'left', va = 'center', rotation = 'vertical')
#fig.text(0.5,0.01, r'$\rm \log_{10}(L_{NUV} \; / \; erg \; s^{-1})$', ha = 'center', va = 'bottom')
fig.text(0.5,0.07, r'$\rm \log_{10}(f_{H} \; / \; nJy)$', ha = 'center', va = 'bottom')

# fig.add_subplot(111, frameon=False)
# plt.xlabel('H-ch1')
# plt.ylabel('ch1-ch2')
# plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
#


fig.savefig(f'figures/beta_rest_FUV_NUV__H__deeper.pdf', bbox_inches="tight")
fig.clf()


binw = 0.25
log10H_bin_edges = np.arange(2., 4.25, binw)

cmap = mpl.cm.plasma


'''
fig = plt.figure(10)
for z in betas.keys():
    medians = np.array([])
    b_84 = np.array([])
    b_16 = np.array([])
    for bin in log10H_bin_edges:
        s = abs(betas[z]['log10H']-(bin+binw/2)) < binw/2

        try:
            print(len(betas[z]['beta']))
            medians = np.append(medians, np.percentile(betas[z]['beta'][s], 50))
            b_84 = np.append(b_84, np.percentile(betas[z]['beta'][s], 84))
            b_16 = np.append(b_16, np.percentile(betas[z]['beta'][s], 16))

        except:
            medians = np.append(medians, np.nan)
            b_84 = np.append(b_84, np.nan)
            b_16 = np.append(b_16, np.nan)

    plt.plot(log10H_bin_edges+binw/2, medians, c=cmap((float(z)-5)/10), alpha=0.6)
    #plt.fill_between(log10H_bin_edges+binw/2, b_16, b_84, color=cmap((float(z)-5)/10), alpha=0.6)


plt.ylabel(r'$\rm \beta$')
plt.xlabel(r'$\rm \log_{10}(f_{H} \; / \; nJy)$')

fig.savefig(f'figures/beta_summary.pdf')
fig.clf()
'''

fig = plt.figure(figsize=(4,6))

for z in beta_summary.keys():
    if z == '5.0' or z == '7.0' or z == '9.0':
        plt.plot(beta_summary[z]['log10H'], beta_summary[z]['beta'], '-', c=cmap(1.5*(float(z) - 5) / 10), alpha=0.6, label=f'z = {int(float(z))}')
        plt.plot(beta_summary[z]['log10H'], beta_nebular[z]['beta'], '--', c=cmap(1.5 * (float(z) - 5) / 10), alpha=0.6)
        plt.plot(beta_summary[z]['log10H'], beta_stellar[z]['beta'], linestyle='dashdot', c=cmap(1.5 * (float(z) - 5) / 10), alpha=0.6)


plt.plot([0.], [0.], c='k', linestyle='-', alpha=0.4, label='stellar + nebular + dust')
plt.plot([0.], [0.], c='k', linestyle='--', alpha=0.4, label='stellar + nebular')
plt.plot([0.], [0.], c='k', linestyle='dashdot', alpha=0.4, label='stellar')


plt.ylim(-2.85, -1.85)
plt.xlim(1.0, 2.8)
plt.ylabel(r'$\rm \beta$')
plt.xlabel(r'$\rm \log_{10}(f_{H} \; / \; nJy)$')

plt.legend(loc='upper left')

fig.savefig(f'figures/beta_summary_V3_deeper.pdf', bbox_inches="tight")
fig.clf()


import _pickle as pickle

pickle.dump(beta_all, open('beta_all.p', 'wb'))

