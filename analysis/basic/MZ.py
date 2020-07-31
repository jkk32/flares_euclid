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

cmap = mpl.cm.plasma
norm = mpl.colors.Normalize(vmin=5., vmax=10.)




# --- read in FLARS data

fl = flares.flares('../../data/flares.hdf5', sim_type='FLARES')
halo = fl.halos

iz = -4

tag = fl.tags[iz]
z = fl.zeds[iz]


# -- define data of interest
# -- CREATE A LOOP OVER PARAMETERS OF INTEREST

M = fl.load_dataset('Mstar_30', arr_type='Galaxy')




df = pd.read_csv('/cosma/home/dp004/dc-wilk2/data/flare/modules/flares/weight_files/weights_grid.txt')
weights = np.array(df['weights'])









fig = plt.figure(figsize=(3,3))

left  = 0.2
bottom = 0.2
width = 0.75
height = 0.75

ax = fig.add_axes((left, bottom, width, height))



for P_type, ls in zip(['S','G'],[':','-']):

    ws, x, y = np.array([]), np.array([]), np.array([])
    for ii in range(len(halo)):
        ws = np.append(ws, np.ones(np.shape(M[halo[ii]][tag]))*weights[ii])
        x = np.append(x, np.log10(M[halo[ii]][tag]))


        y_d = fl.get_particles(p_str=f'{P_type}_Z', halo=halo[ii], tag=tag)
        m = fl.get_particles(p_str=f'{P_type}_Mass', halo=halo[ii], tag=tag)

        y_a = np.array([np.sum(y_d[k][f'{P_type}_Z']*m[k][f'{P_type}_Mass'])/np.sum(m[k][f'{P_type}_Mass']) for k in y_d.keys()])

        y = np.append(y, np.log10(y_a))

    x += 10.


    bins = np.arange(7., 12., 0.2)
    bincen = (bins[:-1]+bins[1:])/2.
    out = flares.binned_weighted_quantile(x,y,ws,bins,[0.84,0.50,0.16])

    N, bin_edges = np.histogram(x, bins=bins)

    Ns = N>5

    ax.plot(bincen[Ns], out[:,1][Ns], c='k', ls = ls)
    ax.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color='k', alpha = 0.2)




ax.set_xlim([7., 12.])

ax.set_xlabel(r'$\rm log_{10}(M^{\star}/M_{\odot})$')
ax.set_ylabel(r'$\rm log_{10}(<Z>)$')

fig.savefig(f'figures/MZ_{z}.pdf')
fig.clf()
