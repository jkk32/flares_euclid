# load and print differences in beta

import _pickle as pickle
import numpy as np

beta_all = pickle.load(open('beta_all.p', 'rb'))

h_lim = 2. # log10(h) limit

diffs = {}
['Pure_Stellar','Intrinsic', 'DustModelI']

for z in beta_all['Intrinsic'].keys():
    s = beta_all['DustModelI'][z]['log10H'] > h_lim
    diffs[z] = {'dust-nebular': beta_all['DustModelI'][z]['beta'][s] - beta_all['Intrinsic'][z]['beta'][s], 'nebular-stellar': beta_all['Intrinsic'][z]['beta'][s] - beta_all['Pure_Stellar'][z]['beta'][s]}


for z in diffs.keys():
    print(f'z = {z}; dust-nebular = {np.mean(diffs[z]["dust-nebular"]):.2f}, nebular-stellar = {np.mean(diffs[z]["nebular-stellar"]):.2f}')