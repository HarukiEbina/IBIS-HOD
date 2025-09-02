#!/usr/bin/env python3
#
# For each setting, save the HOD particles as a .fits file to read in the fiber assignment code
#
import numpy as np

import json
from   astropy.table        import Table, Column
from   scipy.interpolate    import InterpolatedUnivariateSpline as Spline
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

import sys, os
sys.path.append('..')

from rotate_to import rotate_to
from   fake_lbg import MockLBG
from   calc_xi  import calc_xi
from   calc_wt  import calc_wt
from   make_mask import SurveyMask

# from mc_forecast import mc_forecast
from forecast_utils import *


'''
This is the save ran function 
Only for 2D distribution for now. Extending this to 3D is easy
'''

if __name__=="__main__":
    config_file = sys.argv[1]
    with open(config_file) as f:
        config = json.load(f)
    print(config_file)
    chi0 = config['chi0']
    chirange = config['chirange']
    hod_params = config['hod_params']
    ndensity_band = config['ndensity_band'] # target density
    ndensity_band = np.array(ndensity_band) 
    # m_bands = ['M411','M438','M464','M490','M517']
    m_bands = config['bands']
    chivec = np.array(config['chi_vector'])
    name = config['name']
    bin_width = config['bin_width']
    nbox = config['nbox']
    maskfn = config['maskfn']
    use_mask = config['use_mask']
    fint_list = np.array(config['interloper_band'])
    gal_type = config['gal_type']

    ntarg_ratio = config['ntarg_ratio'] # use this for now
    # nspec_ratio = config['nspec_ratio']
    custom_ntarg = config['custom_ntarg'] # use this for now
    # custom_nspec = config['custom_nspec']

    cra = config['cra']
    cdc = config['cdc']
    diam = config['diam']

    zbox = config['zbox']
    sim_name = config['sim_name']

    if ntarg_ratio == -1 and custom_ntarg == -1: ntarg_ratio = 1.
    elif ntarg_ratio != -1: assert custom_ntarg == -1
    else: ntarg_ratio = custom_ntarg / np.sum(ndensity_band*(1-fint_list))

    hod_string = get_hod_string(hod_params)
    
    if use_mask: 
        assert os.path.exists(maskfn)
        mask = SurveyMask(maskfn)
        area = mask.area()
    else: 
        area = (diam * 180/np.pi /2)**2 * np.pi
    
    def in_field(tab):
        if use_mask:
            nra,ndc    = rotate_to(tab['RA'],tab['DEC'],cra,cdc)
            ww         = mask(nra,ndc)
            return ww, nra[ww], ndc[ww]
        else:
            # rad2 = ( (tab['RA']-Lx)**2 + (tab['DEC'])**2 )*(np.pi/180.)**2
            rad2 = ( (tab['RA'])**2 + (tab['DEC'])**2 )*(np.pi/180.)**2
            nra,ndc    = rotate_to(tab['RA'],tab['DEC'],cra,cdc)
            ww = rad2<diam**2/4
            return ww, nra[ww], ndc[ww]
        
    outdir = '{:s}/z{:.1f}/{:s}/{:s}/inputcats'.format(sim_name,zbox,name,hod_string)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    diam_ang = diam * 180/np.pi 
    rng    = np.random.default_rng(2) # use a different seed than save_hod just in case
    chimin, chimax = chi0 -chirange/2, chi0 + chirange/2
    # nran = 8000000
    # nran = 4000000
    nran = 1000000
    for ibox in range(nbox):
        ran  = {}
        ran['RA' ] = rng.uniform(low=-diam_ang/2.,high=diam_ang/2. ,size=nran) # + Lx
        ran['DEC'] = rng.uniform(low=-diam_ang/2.,high=diam_ang/2. ,size=nran)
        ran['CHI'] = rng.uniform(low=chimin,high=chimax,size=nran)

        ww, nra, ndc = in_field(ran)
        ran['RA'] = nra
        ran['DEC'] = ndc
        ran['CHI'] = ran['CHI'][ww]
        print("Random size ",len(ran['RA']),flush=True)
    
        ran['SEL'] = np.zeros(ran['RA'].size,dtype=bool)
        rand = rng.uniform(low=0,high=1,size=ran['RA'].size)
    
        vmax = 0
        for i, band in enumerate(m_bands):
            dNdc = np.array(config['dNdchi_band'])[i]
            dNdc_func = interp1d(chivec,dNdc,bounds_error=False,fill_value=0)
            vmax = max(np.max(dNdc_func(ran['CHI'])),vmax)
    
        factor = 1/vmax # adjust maximum probability to be 1
        for i, band in enumerate(m_bands):
            dNdc = np.array(config['dNdchi_band'])[i]
            dNdc_func = interp1d(chivec,dNdc,bounds_error=False,fill_value=0)
            ran['SEL_%s'%band] = rand < dNdc_func(ran['CHI'] ) * factor 
            ran['SEL'] |= ran['SEL_%s'%band]
        ww = ran['SEL']
        for key in ran.keys():
            ran[key] = ran[key][ww]

        
        ran_table = Table()
        ran_table['RA'] = ran['RA' ]
        ran_table['DEC'] = ran['DEC']
        ran_table['CHI'] = ran['CHI']
        for i, band in enumerate(m_bands):
            ran_table['SEL_%s'%band] = ran['SEL_%s'%band]
        ran_table.write(os.path.join(outdir,'ran%d.fits'%ibox),overwrite=True)

        del ran, ran_table




