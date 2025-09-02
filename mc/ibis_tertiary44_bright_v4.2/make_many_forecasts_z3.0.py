
basedir = "/pscratch/sd/h/hebina/AbacusLBG/"
import sys, os
from astropy.table import Table, vstack
from astropy.io import fits

for subdir in ["",'ibis_tertiary44/LAE_auto_v2']:
    sys.path.append(basedir+subdir)
from misc import *
from matplotlib import pyplot as plt
sys.path.append('..')
from mc_forecast import mc_forecast
from forecast_utils import *
import time

cosmo = get_cosmo_P18()
h = cosmo.h()

fn = "/pscratch/sd/h/hebina/IBIS/selection_clustering/catalogs/clustering_cat_v4_bright.fits"
# fn = "/global/cfs/cdirs/desi/users/raichoor/laelbg/ibis/analysis/daily-tmp/ibis-xmm-ar-djs-he-rr.fits"
d = Table.read(fn)
radius = 1.4
ra0,dec0 = 35.75, -4.75
area = 1.4**2 * np.pi
coord_cut = angular_distance(d['RA'],d['DEC'],ra0,dec0)<radius
d = d[coord_cut]

fn = "/pscratch/sd/h/hebina/IBIS/selection_clustering/catalogs/redshift_cat_v4_bright.fits"
# fn = "/pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto/ibis-xmm-ar-djs-he-rr.fits"
d_spec = Table.read(fn)
radius = 1.4
ra0,dec0 = 35.75, -4.75
area = 1.4**2 * np.pi
coord_cut = angular_distance(d_spec['RA'],d_spec['DEC'],ra0,dec0)<radius
d_spec = d_spec[coord_cut]

use_mask = True
m_bands = ['M411','M438','M464','M490','M517']
band_z = {
    'M411':[2.26,2.56],
    'M438':[2.47,2.77],
    'M464':[2.68,2.98],
    'M490':[2.89,3.2],
    'M517':[3.1,3.41]
}

# if len(sys.argv) > 1:
    # lgMc_list = [float(sys.argv[1])]
# else:
    # lgMc_list = [11.00,11.25,11.50,11.75,12.0,12.25,12.5]


# lgMc_list = [11.00,11.25,11.50,11.75,12.0,12.25,12.5]
lgMc_list = [11.25,11.50,11.75,12.0,12.25,12.5]
lgMc_list = [11.00] # done with sbatch
# lgMc_list = [11.25] # done with login
alph_list = [.33,.5,.66]#,1.]
sigm_list = [.33,.5,.66]#,1.]
kapp_list = [1]
plat_list = [5,10]

# lgMc_list = [11.75]
# alph_list = [0.0,.33,.5,.66,1.]
# sigm_list = [.66]
# kapp_list = [0,1]
# plat_list = [5,10]


# enumerate through HOD parameters and generate HODs
for lgMcut in lgMc_list:
  for plateau in plat_list:
    lgM1 = lgMcut + np.log10(plateau)
    for alph in alph_list:
      for sigm in sigm_list:
        for kappa in kapp_list:
          # hod_params = {'logM_cut':11.75,'logM1':12.45,'sigma':0.66,'kappa':1.00,'alpha':1.0}
          hod_params = {'logM_cut':lgMcut,'logM1':lgM1,'sigma':sigm,'kappa':kappa,'alpha':alph}           
          os.system('date')
          hod_string = get_hod_string(hod_params)
          print(hod_string)
            
          # # first two bins
          # zbox = 2.5
          # simfn = f'../../hod/boxes/z{zbox:.1f}/LRG/{hod_string}.asdf'
          # print(f'checking sim box at {simfn}')
          # while True: 
          #     # check every minute to see if hod for this forecast is done
          #     if os.path.exists(simfn): break
          #     time.sleep(10)
          # print(f'box exists. proceeding')
          # forecast = mc_forecast.initialize('ibis_tertiary44_bright_v4.2',cosmo,d,d_spec , '../../ibis_tertiary44_msk.fits',35.75,-4.75, 2.8*np.pi/180,
          #                        band_z['M411'][0], band_z['M517'][1], hod_params, use_mask, overwrite=False,ntarg_ratio=1.,
          #                        fsamp=d_spec['RA'].size/d['RA'].size,basedir='..',nbox=256,m_bands=['M411','M438'],zbox=2.5)
          # # forecast.make_hod(overwrite=False,jobtype='debug')
          # # forecast.make_hod(overwrite=False,jobtype='regular')
          # forecast.make_hod(overwrite=False,sbatch=False)
          # print('{:s} z1 submit job'.format(get_hod_string(hod_params)))
          # os.system('date')
            
          # next three bins
          zbox = 3.0
          simfn = f'../../hod/boxes/z{zbox:.1f}/LRG/{hod_string}.asdf'
          print(f'checking sim box at {simfn}')
          while True: 
              # check every minute to see if hod for this forecast is done
              if os.path.exists(simfn): break
              time.sleep(10)
          print(f'box exists. proceeding')
          forecast = mc_forecast.initialize('ibis_tertiary44_bright_v4.2',cosmo,d,d_spec , '../../ibis_tertiary44_msk.fits',35.75,-4.75, 2.8*np.pi/180,
                                 band_z['M411'][0], band_z['M517'][1], hod_params, use_mask, overwrite=False,ntarg_ratio=1.,
                                 fsamp=d_spec['RA'].size/d['RA'].size,basedir='..',nbox=256,m_bands=['M464','M490','M517'],zbox=3.0)
          # forecast.make_hod(overwrite=False,jobtype='debug')
          # forecast.make_hod(overwrite=False,jobtype='regular')
          forecast.make_hod(overwrite=False,sbatch=False)
          print('{:s} z2 submit job'.format(get_hod_string(hod_params)))
          os.system('date')

for lgMcut in lgMc_list:
  for plateau in plat_list:
    lgM1 = lgMcut + np.log10(plateau)
    for alph in alph_list:
      for sigm in sigm_list:
        for kappa in kapp_list:
          hod_params = {'logM_cut':lgMcut,'logM1':lgM1,'sigma':sigm,'kappa':kappa,'alpha':alph}        
            
          # # first two bins
          # forecast = mc_forecast.initialize('ibis_tertiary44_bright_v4.2',cosmo,d,d_spec , '../../ibis_tertiary44_msk.fits',35.75,-4.75, 2.8*np.pi/180,
          #                        band_z['M411'][0], band_z['M517'][1], hod_params, use_mask, overwrite=False,ntarg_ratio=1.,
          #                        fsamp=d_spec['RA'].size/d['RA'].size,basedir='..',nbox=256,m_bands=['M411','M438'],zbox=2.5)
          # while True: 
          #     # check every minute to see if hod for this forecast is done
          #     if forecast.check_hod_exists(): break
          #     time.sleep(10)
          # # make catalog including interlopers when the HODs are done
          # forecast.make_cat(overwrite=False)

          # # compute wt
          # # forecast.compute_wt()
          # forecast.compute_wR(overwrite=False)
          # # # # compute xi
          # # forecast.compute_xi()
          # # for j, band in enumerate(m_bands):
          # #     forecast.compute_xi(zbins=[j])
            
          # print('{:s} z1 done'.format(get_hod_string(hod_params)))

          # next three bins
          forecast = mc_forecast.initialize('ibis_tertiary44_bright_v4.2',cosmo,d,d_spec , '../../ibis_tertiary44_msk.fits',35.75,-4.75, 2.8*np.pi/180,
                                 band_z['M411'][0], band_z['M517'][1], hod_params, use_mask, overwrite=False,ntarg_ratio=1.,
                                 fsamp=d_spec['RA'].size/d['RA'].size,basedir='..',nbox=256,m_bands=['M464','M490','M517'],zbox=3.0)

          while True: 
              # check every minute to see if hod for this forecast is done
              if forecast.check_hod_exists(): break
              time.sleep(10)
          # make catalog including interlopers when the HODs are done
          forecast.make_cat(overwrite=False)

          # compute wt
          # forecast.compute_wt()
          forecast.compute_wR(overwrite=False)
          # # # compute xi
          # forecast.compute_xi()
          # for j, band in enumerate(m_bands):
          #     forecast.compute_xi(zbins=[j])
            
          print('{:s} z2 done'.format(get_hod_string(hod_params)))

