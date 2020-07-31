

import numpy as np
import flares


fl = flares.flares('../../data/flares.hdf5', sim_type='FLARES')
halo = fl.halos

tag = fl.tags[0]
isim = 0


p_str = 'G_Z'

d = fl.get_particles(p_str=p_str, halo=halo[isim], tag=tag)
m = fl.get_particles(p_str='G_Mass', halo=halo[isim], tag=tag)

print(d.keys())
print(d[0].keys())



keys = np.array(list(d.keys()))

print(keys)

a = np.array([np.mean(d[k][p_str]*m[k]['G_Mass'])/np.sum(m[k]['G_Mass']) for k in d.keys()])

print(a.shape)
