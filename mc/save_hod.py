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
import asdf

# from mc_forecast import mc_forecast
from forecast_utils import *


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
    # else: ntarg_ratio = custom_ntarg / np.sum(ndensity_band)
    else: ntarg_ratio = custom_ntarg / np.sum(ndensity_band*(1-fint_list))
    # else: 
        # if custom_ntarg!=-1: ntarg_ratio = custom_ntarg / np.sum(ndensity_band*(1-fint_list))
        # else: ntarg_ratio = custom_ntarg / np.sum(ndensity_band*(1-fint_list))

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
    
    if 'Acent' in hod_params.keys() or 'Bcent' in hod_params.keys() or 'Asat' in hod_params.keys() or 'Bsat' in hod_params.keys(): want_AB = True
    else: want_AB = False
    if 's' in hod_params.keys() or 's_v' in hod_params.keys() or 's_p' in hod_params.keys() or 's_r' in hod_params.keys(): want_rank = True
    else: want_rank = False

    if want_AB: yamlfn = 'lae_AB.yaml'
    elif want_rank: yamlfn = 'lae_ranks.yaml'
    else: yamlfn = 'lae_base.yaml'
    rng    = np.random.default_rng(1)
    lbgs   = MockLBG(yamlfn,None,chi0,chirange,zbox=zbox,sim_name=sim_name,gal_type=gal_type)
    
    lbgs.set_hod(hod_params)
    print(lbgs.sim_params)
    # lbgs.generate()
    lbgs.generate(load_path=f'/pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/hod/boxes/z{zbox:.1f}/{gal_type}/{hod_string}.asdf')
    # lbgs.assign_lum(0.25) # this means nothing
    lbgs.select(diam,[0.,0.,0.])

    print('galaxy catalog size',lbgs.xang.size)
    # Lside  = lbgs.d['Lside']
    # Lx     = Lside/chi0 * 180./np.pi
    # Ly     = Lside/chi0 * 180./np.pi
    # ichi   = 1.0/  chi0 * 180./np.pi
    diam_ang = diam * 180/np.pi 
    sq_area = (diam_ang)**2
    # Now we want to determine the sampling fraction to
    # get the right angular number density.
    # ntarget, nbar = ndensity, []
    nbar = []
    print('ntarget',ndensity_band)
    for i in range(25):
        # Generate the galaxies.
        offset = rng.uniform(low=-0.5,high=0.5,size=3)
        lbgs.select(diam,offset)
        nobj   = lbgs.d['nkeep']
        nbar.append(nobj / (diam_ang**2))
    # fsamp_band = ndensity_band / np.median(nbar)
    # chimin,chimax = np.min(lbgs.zpos+chi0),np.max(lbgs.zpos+chi0)
    # chimin,chimax = np.min(lbgs.zchi),np.max(lbgs.zchi)
    chimin, chimax = chi0 -chirange, chi0 + chirange
    print('chis',chimin,chimax)
    # Generate a uniform random catalog.
    
    # nran = 8000000
    # nran = 4000000
    nran = 1000000
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
    # for i, band in enumerate(m_bands):
    #     dNdc = np.array(config['dNdchi_band'])[i]
    #     dNdc_func = interp1d(chivec,dNdc,bounds_error=False,fill_value=0)
    #     # ran['SEL_%s'%band] = rand < ndensity_band[i]/np.max(ndensity_band) * dNdc_func(ran['CHI']) / np.max(dNdc_func(ran['CHI'])) # old (wrong?) normalization
    #     ran['SEL_%s'%band] = rand < dNdc_func(ran['CHI']) / ( np.median(nbar) * bin_width / (2*chirange) ) * 100 # same normalization as dat; ran should be 100x more populous as dat
    #     ran['SEL'] |= ran['SEL_%s'%band]

    vmax = 0
    for i, band in enumerate(m_bands):
        dNdc = np.array(config['dNdchi_band'])[i]
        dNdc_func = interp1d(chivec,dNdc,bounds_error=False,fill_value=0)
        vmax = max(np.max(dNdc_func(ran['CHI']) / ( np.median(nbar) * bin_width / (2*chirange) )),vmax)

    factor = 1/vmax # adjust maximum probability to be 1
    for i, band in enumerate(m_bands):
        dNdc = np.array(config['dNdchi_band'])[i]
        dNdc_func = interp1d(chivec,dNdc,bounds_error=False,fill_value=0)
        ran['SEL_%s'%band] = rand < dNdc_func(ran['CHI']) / ( np.median(nbar) * bin_width / (2*chirange) ) * factor 
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
    ran_table.write(os.path.join(outdir,'ran.fits'),overwrite=True)

    del ran, ran_table

    rval,wts ,ngals = None,[],[]
    for i in range(nbox):
        print(i)
        # Generate the galaxies.
        offset = rng.uniform(low=-0.5,high=0.5,size=3)
        lbgs.select(diam,offset)
        dat        = {}
        dat['RA' ] = lbgs.xang # + Lx
        dat['DEC'] = lbgs.yang
        dat['CHI'] = lbgs.zchi

        # use for check - otherwise take out to speed up
        # dat['x'] = lbgs.xpos
        # dat['y'] = lbgs.ypos
        # dat['z'] = lbgs.zpos


        dat['SEL'] = np.zeros(dat['RA'].size,dtype=bool)
        '''
        for j, band in enumerate(m_bands):
            rand = rng.uniform(low=0,high=1,size=dat['RA'].size)
            dNdc = np.array(config['dNdchi_band'])[j]
            dNdc_func = interp1d(chivec,dNdc,bounds_error=False,fill_value=0)
            # dat['SEL_%s'%band] = rand < dNdc_func(dat['CHI']) / ( np.median(nbar) * bin_width / (2*chirange) ) * ntarg_ratio # allow ndensity to modulate later
            dat['SEL_%s'%band] = rand < dNdc_func(dat['CHI']) / ( np.average(nbar) * bin_width / (2*chirange) ) * ntarg_ratio # allow ndensity to modulate later
            dat['SEL'] |= dat['SEL_%s'%band]
        ww = dat['SEL']
        '''
        
        vmax = 0
        for j, band in enumerate(m_bands):
            dNdc = np.array(config['dNdchi_band'])[j]
            dNdc_func = interp1d(chivec,dNdc,bounds_error=False,fill_value=0)
            vmax = max(np.max(dNdc_func(dat['CHI']) / ( np.median(nbar) * bin_width / (2*chirange) )),vmax)
    
        factor = 1/vmax # adjust maximum probability to be 1

        # downsample to dNdchi shape
        for j, band in enumerate(m_bands):
            rand = rng.uniform(low=0,high=1,size=dat['RA'].size)
            dNdc = np.array(config['dNdchi_band'])[j]
            dNdc_func = interp1d(chivec,dNdc,bounds_error=False,fill_value=0)
            dat['SEL_%s'%band] = rand < dNdc_func(dat['CHI']) / ( np.average(nbar) * bin_width / (2*chirange) ) * factor # allow ndensity to modulate later
            dat['SEL'] |= dat['SEL_%s'%band]

        print('compare data density and mock density',np.sum(ndensity_band*(1-fint_list))* ntarg_ratio,(dat['SEL'].sum()/sq_area)) # area right now is a square
        # assert np.sum(ndensity_band*(1-fint_list))<(dat['SEL'].sum()/sq_area)
        if np.sum(ndensity_band*(1-fint_list)) * ntarg_ratio>(dat['SEL'].sum()/sq_area): print('warning')
        # downsample to number density
        rand = rng.uniform(low=0,high=1,size=dat['RA'].size)
        ww = rand < np.sum(ndensity_band*(1-fint_list))/(dat['SEL'].sum()/sq_area) * ntarg_ratio
        ww &= dat['SEL']
        
        for key in dat.keys():
            dat[key] = dat[key][ww]
        print((dat['SEL'].sum()/sq_area))

        if i==0:        
            fig, ax = plt.subplots(1,1,figsize=(5,3.5))
            bins = np.arange(round(chimin,-1),round(chimax,-1), 10)
            for j, band in enumerate(m_bands):
                plt.hist(dat['CHI'][dat['SEL_%s'%band]],bins=bins,color='C%d'%j,alpha=.5)
            plt.xlim(chimin,chimax)
            plt.savefig('plots/tmp.png')
            plt.show()

        print('mock density after downsampling',dat['RA'].size/sq_area)
        
        ww, nra, ndc = in_field(dat)
        dat['RA' ] = nra
        dat['DEC'] = ndc
        for key in dat.keys():
            if key=='RA' or key=='DEC': continue
            dat[key] = dat[key][ww]
        # dat['CHI'] = dat['CHI'][ww]
        print('mock density after in_field',dat['RA'].size/area)

        dat_table = Table()
        for key in dat.keys():
            dat_table[key] = dat[key]
        dat_table.write(os.path.join(outdir,'dat%d.fits'%i),overwrite=True)

        del dat, dat_table

    del lbgs