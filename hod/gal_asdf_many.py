#!/usr/bin/env python3
#
import numpy as np
import yaml
import json

import asdf
from abacusnbody.hod.abacus_hod import AbacusHOD
from abacusnbody.metadata       import get_meta
from abacusnbody.analysis.power_spectrum       import get_field,get_interlaced_field_fft,_normalize,get_W_compensated
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
from scipy.interpolate import interp1d

from matplotlib import pyplot as plt

from scipy.fft import fftn, ifftn, fftfreq

import os, pathlib, sys
# These are the keys we will copy from the simulation meta-data into
# the output JSON file.
simkeys = ['n_s', 'omega_b', 'omega_cdm', 'omega_ncdm', 'N_ncdm', 'N_ur',\
           'H0', 'w0', 'wa', 'w',\
           'Omega_DE', 'Omega_K', 'Omega_M', 'Omega_Smooth', \
           'Redshift', 'ScaleFactor', \
           'OmegaNow_DE', 'OmegaNow_K', 'OmegaNow_m', \
           'f_growth', 'fsmooth', 'Growth', 'Growth_on_a_n', \
           'SimComment', 'SimName', 'SimSet', \
           'BoxSize', 'NP', 'BoxSizeHMpc', 'HubbleTimeHGyr', \
           'ParticleMassHMsun']

    #
def get_hod_string(hod_params,use_extension=False):
    # hod_keys = ['logM_cut','logM1','sigma','kappa','alpha'] # base hod
    # if use_extension: hod_keys += ['Acent','Bcent','Asat','Bsat','alpha_c','alpha_s','s','s_v','s_p','s_r']

    hod_keys = hod_params.keys()
    out_str = ''
    for i,key in enumerate(hod_keys): 
        if i==0: out_str += '{:s}_{:.2f}'.format(key,hod_params[key])
        else: out_str += '_{:s}_{:.2f}'.format(key,hod_params[key])
    return out_str

#
# Load the config file and parse in relevant parameters
path2config= './lae_ext.yaml'
config     = yaml.safe_load(open(path2config))

zbox = float(sys.argv[1])
config['sim_params']['z_mock'] = zbox

sim_params = config['sim_params']
HOD_params = config['HOD_params']
clustering_params = config['clustering_params']
#
print(sim_params)
# Get the metaparameters for the simulation.
meta = get_meta(sim_params['sim_name'],redshift=sim_params['z_mock'])
#
# additional parameter choices
want_zcv      = False
want_rsd      = HOD_params['want_rsd'] & False
write_to_disk = HOD_params['write_to_disk']

# # parameters for fft field
# nmesh = 637
# compensated = False; interlaced = False
# paste = 'TSC'


lgMc_list = [11.00,11.25,11.50,11.75,12.0,12.25,12.5]
alph_list = [.33,.5,.66]#,1.]
sigm_list = [.33,.5,.66]#,1.]
kapp_list = [1]
plat_list = [5,10]

outsuf = 's' if want_rsd else 'r'
outsim = sim_params['sim_name']
# outsim = outsim[outsim.find("_c0"):] + '_z%.1f_'%zbox
if want_zcv: outsuf += '_zcv'

outfn_base = 'boxes/z%.1f/LRG/'%zbox
if not os.path.exists(outfn_base):
    os.makedirs(outfn_base)

# os.makedirs(outfn, exist_ok=True)          # classic stdlib way
# pathlib.Path(outfn).parent.mkdir(parents=True, exist_ok=True)

hodkeys = ['logM_cut','logM1','sigma','kappa','alpha']
kk = 0

for lgMcut in lgMc_list:
  for alph in alph_list:
    for sigm in sigm_list:
      for kappa in kapp_list:
        for plateau in plat_list:
            lgM1 = lgMcut + np.log10(plateau)
            vec = ([lgMcut,lgM1, sigm,kappa,alph])
    
            for key,val in zip(hodkeys,vec):
              HOD_params['LRG_params'][key] = val            
            hod      = [HOD_params['LRG_params'][k] for k in hodkeys]
            print(HOD_params)
            tmp = {}
            for key in hodkeys: tmp[key] = HOD_params['LRG_params'][key]
            hod_string = get_hod_string(tmp)
    
            outfn = outfn_base + f'{hod_string}.asdf'
            if os.path.exists(outfn): continue
    
            newBall  = AbacusHOD(sim_params,HOD_params,clustering_params)
            mock_dict= newBall.run_hod(newBall.tracers,want_rsd,\
                                       write_to_disk,Nthread=16)
    
            ### directly output the mock_dict
            af = asdf.AsdfFile(mock_dict)
            af.write_to(outfn) 
            af.close()
    
            del newBall, mock_dict
    
