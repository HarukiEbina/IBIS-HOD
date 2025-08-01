import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from astropy.io import fits
from classy import Class
from scipy.interpolate import interp1d

from scipy.interpolate import make_smoothing_spline
import os, sys
sys.path.append('..')
from make_mask import SurveyMask
from misc import *

'''
Compute the interloper fraction of dataset
'''

os.makedirs('fint_v4', exist_ok=True)

m_bands = ['M411','M438','M464','M490','M517']
band_z = {
    'M411':[2.26,2.56],
    'M438':[2.47,2.77],
    'M464':[2.68,2.98],
    'M490':[2.89,3.2],
    'M517':[3.1,3.41]
}

cosmo = get_cosmo_P18()
h = cosmo.h()
def chi(zz):
    if zz>0:return (cosmo.comoving_distance(zz)*h)
    else: return (0)

fn = "/pscratch/sd/h/hebina/IBIS/selection_clustering/catalogs/clustering_cat_v4.fits"
d = Table.read(fn)

# cut to 1.4deg mask
mask = SurveyMask('../ibis_tertiary44_msk.fits')
ww = mask(d['RA'],d['DEC'])
d = d[ww]

# define min and max of redshift ranges
zmin, zmax = band_z['M411'][0], band_z['M517'][1]
chimin, chimax = chi(zmin),chi(zmax)

# restrict to objects with redshift and 2hrs observation
d = d[d['EFFTIME']>7200]
n_tot = d['RA'].size
d = d[d['RRJM_DELTACHI2']>30]
n_rr = d['RA'].size
# restrict to objects with desired redshift range (2.26,3.41)
d = d[(d['RRJM_Z']>zmin)&(d['RRJM_Z']<zmax)]
n_good = d['RA'].size

print('pessimistic %.3f'%(1-n_good/n_tot))
print('optimistic %.3f'%(1-n_good/n_rr))


opt_list = []
pes_list = []
for band in m_bands:
    fn = "/pscratch/sd/h/hebina/IBIS/selection_clustering/catalogs/clustering_cat_v4.fits"
    d = Table.read(fn)
    
    # cut to 1.4deg
    mask = SurveyMask('../ibis_tertiary44_msk.fits')
    ww = mask(d['RA'],d['DEC'])
    d = d[ww]
    
    # define min and max of redshift ranges
    if band in band_z.keys():
        zmin, zmax = band_z[band][0], band_z[band][1]
        d = d[d['SEL_%s'%band]]
    else:
        zmin, zmax = band_z['M411'][0], band_z['M517'][1]
    chimin, chimax = chi(zmin),chi(zmax)

    n_targ = d['RA'].size
    # restrict to objects with redshift and 2hrs observation
    d = d[d['EFFTIME']>7200]
    n_tot = d['RA'].size
    d = d[d['RRJM_DELTACHI2']>30]
    n_rr = d['RA'].size
    # restrict to objects with desired redshift range (2.26,3.41)
    d = d[(d['RRJM_Z']>zmin)&(d['RRJM_Z']<zmax)]
    n_good = d['RA'].size

    print(band)
    print('pessimistic %.3f'%(1-n_good/n_tot))
    print('optimistic %.3f'%(1-n_good/n_rr))

    opt_list.append(1-n_good/n_rr)
    pes_list.append(1-n_good/n_tot)

np.savetxt('fint_v4/optimistic.txt',opt_list)
np.savetxt('fint_v4/pessimistic.txt',pes_list)
    
opt_list = []
pes_list = []
for band in m_bands:
    fn = "/pscratch/sd/h/hebina/IBIS/selection_clustering/catalogs/clustering_cat_v4_bright.fits"
    d = Table.read(fn)
    
    # cut to 1.4deg
    mask = SurveyMask('../ibis_tertiary44_msk.fits')
    ww = mask(d['RA'],d['DEC'])
    d = d[ww]
    
    # define min and max of redshift ranges
    if band in band_z.keys():
        zmin, zmax = band_z[band][0], band_z[band][1]
        d = d[d['SEL_%s'%band]]
    else:
        zmin, zmax = band_z['M411'][0], band_z['M517'][1]
    chimin, chimax = chi(zmin),chi(zmax)
    
    n_targ = d['RA'].size

    # restrict to objects with redshift and 2hrs observation
    d = d[d['EFFTIME']>7200]
    n_tot = d['RA'].size
    d = d[d['RRJM_DELTACHI2']>30]
    n_rr = d['RA'].size
    # restrict to objects with desired redshift range (2.26,3.41)
    d = d[(d['RRJM_Z']>zmin)&(d['RRJM_Z']<zmax)]
    n_good = d['RA'].size

    print(band)
    print('pessimistic %.3f'%(1-n_good/n_tot))
    print('optimistic %.3f'%(1-n_good/n_rr))

    opt_list.append(1-n_good/n_rr)
    pes_list.append(1-n_good/n_tot)


np.savetxt('fint_v4/optimistic_bright.txt',opt_list)
np.savetxt('fint_v4/pessimistic_bright.txt',pes_list)


