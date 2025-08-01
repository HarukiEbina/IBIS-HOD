import numpy as np

import json
from   astropy.table        import Table, Column
from   scipy.interpolate    import InterpolatedUnivariateSpline as Spline
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from scipy.interpolate import make_smoothing_spline

import sys, os
# basedir = "/pscratch/sd/h/hebina/AnalyzeLAE/mocks/"
# for subdir in [""]:
#     sys.path.append(basedir+subdir)
basedir = "/pscratch/sd/h/hebina/AbacusLBG/"
for subdir in ["",'ibis_tertiary44/LAE_auto_v2']:
    sys.path.append(basedir+subdir)

from rotate_to import rotate_to
from   fake_lbg import MockLBG
from   calc_xi  import calc_xi
from   calc_wt  import calc_wt
from   make_mask import SurveyMask
from misc import  *

# from mc_forecast import mc_forecast
# from forecast_utils import *
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
# fn = "/global/cfs/cdirs/desi/users/raichoor/laelbg/ibis/analysis/daily-tmp/ibis-xmm-ar-djs-he-rr.fits"
d = Table.read(fn)
radius = 1.4
ra0,dec0 = 35.75, -4.75
area = 1.4**2 * np.pi
coord_cut = angular_distance(d['RA'],d['DEC'],ra0,dec0)<radius
d = d[coord_cut]

fn = "/pscratch/sd/h/hebina/IBIS/selection_clustering/catalogs/redshift_cat_v4.fits"
# fn = "/pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto/ibis-xmm-ar-djs-he-rr.fits"
d_spec = Table.read(fn)
radius = 1.4
ra0,dec0 = 35.75, -4.75
area = 1.4**2 * np.pi
coord_cut = angular_distance(d_spec['RA'],d_spec['DEC'],ra0,dec0)<radius
d_spec = d_spec[coord_cut]

zmin, zmax = band_z['M411'][0],band_z['M517'][1]
chimin, chimax = chi(zmin), chi(zmax)

bin_width = 10
chivec = np.linspace(chimin,chimax,100000)

f2_dict = {}
for i, band in enumerate(m_bands):
    print(band)
    dtmp_spec = d_spec[d_spec['SEL_%s'%band]]
    dtmp = d[d['SEL_%s'%band]]
    # ndensity = dtmp['RA'].size/mask.area()
    ndensity = dtmp['RA'].size/6.11
    
    ### dNdchi
    zmin_tmp, zmax_tmp = band_z[band][0], band_z[band][1]
    chimin_tmp, chimax_tmp = chi(zmin_tmp), chi(zmax_tmp)

    bins = np.arange(chimin,chimax, bin_width) # spaced out by chi=10
    hchi = np.histogram(dtmp_spec['CHI'],bins = bins) # normalize
    x = (hchi[1][1:]+hchi[1][:-1])/2
    y = hchi[0]
    spl = make_smoothing_spline(x, y, lam=1e3)
    probabilities = spl(chivec)
    probabilities[probabilities<0] = 0

    # recent change
    probabilities[chivec<chimin_tmp] *= 0
    probabilities[chivec>chimax_tmp] *= 0
    

    integral = np.trapz(probabilities,chivec)
    probabilities /= integral
    print(np.trapz(probabilities**2,chivec))
    f2_dict[band] = np.trapz(probabilities**2,chivec)

import json
json.dump(f2_dict,open('f2.json','w'))