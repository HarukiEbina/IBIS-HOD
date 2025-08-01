import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from astropy.io import fits
from classy import Class
from scipy.interpolate import interp1d


import os, sys,json
sys.path.append('..')
from make_mask import SurveyMask
from misc import *
from calc_wt import *


# compute the 2D angular correlation function of of clustering sample, with 5 bins

os.makedirs('plots', exist_ok=True)
os.makedirs('wR_v4', exist_ok=True)


cosmo = get_cosmo_P18()
h = cosmo.h()

m_bands = ['M411','M438','M464','M490','M517']
band_z = {
    'M411':[2.26,2.56],
    'M438':[2.47,2.77],
    'M464':[2.68,2.98],
    'M490':[2.89,3.2],
    'M517':[3.1,3.41]
}

def chi(zz):
    if zz>0:return (cosmo.comoving_distance(zz)*h)
    else: return (0)

radius = 1.4
ra0,dec0 = 35.75, -4.75

loc_dat = '/pscratch/sd/h/hebina/IBIS/selection_clustering/catalogs/clustering_cat_v4.fits' # target file
loc_ran = '../ibis-xmm.ran.fits' # random file
dat = Table.read(loc_dat)
dat = dat[angular_distance(dat['RA'],dat['DEC'],ra0,dec0)<radius]
ran = Table.read(loc_ran)


### 5 bin analysis
dict1 = {}
fig, ax = plt.subplots(1,1,figsize=(6,4))
rbins = np.geomspace(0.4,45,9)[:-1]
# rbins = np.geomspace(0.4,40,9)
dict1['r'] = np.sqrt(rbins[1:]*rbins[:-1]).tolist()

for i,band in enumerate(m_bands):
    zmin, zmax = band_z[band][0],band_z[band][1]
    zmid = (zmin+zmax)/2
    chimin, chimax = chi(zmin),chi(zmax)
    bins = rbins / (cosmo.comoving_distance(zmid) * h) * 180/np.pi

    # tdat = dat[(dat['RRJM_Z']>zmin)&(dat['RRJM_Z']<zmax)]
    tdat = dat[dat['SEL_%s'%band]]
    # tran = ran[(ran['CHI']>chimin)&(ran['CHI']<chimax)]
    tran = ran
    print(tdat['RA'].shape,tran['RA'].shape)
    twt, wwt = calc_wt(tdat,tran[::10],bins=bins)
    # dict1['theta'] = twt.tolist()
    dict1['z%d'%(i+1)] = {
        'zmin': zmin,
        'zmax': zmax,
        'wt': wwt.tolist(),
    }
    # ax.plot(twt,wwt,marker='.',label=m_bands[i])
    ax.plot(dict1['r'],np.array(dict1['r'])*wwt,marker='.',label=m_bands[i])
    
ax.legend()
ax.grid()
ax.set_xscale('log')
ax.set_xlabel(r'$R$')
ax.set_ylabel(r'$R w_\theta(R)$')
fig.tight_layout()    
fig.savefig('plots/wR_5bin_v4.png')
fig.show()
out_file = open("wR_v4/lae_5bins.json", "w")
json.dump(dict1, out_file, indent = 6)
out_file.close()

