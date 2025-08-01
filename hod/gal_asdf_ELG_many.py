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
#
# Load the config file and parse in relevant parameters
path2config= './lae_ext.yaml'
config     = yaml.safe_load(open(path2config))

zbox = float(sys.argv[1])
config['sim_params']['z_mock'] = zbox

config['HOD_params']['tracer_flags']['LRG'] = False
config['HOD_params']['tracer_flags']['ELG'] = True

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

lgMc_list = [11.00,11.25,11.50,11.75,12.0]
# lgMc_list = list(reversed(lgMc_list))
# lgMc_list += [11.3,11.4,11.6,11.7,11.8,11.9]
# alph_list = [.33,.5,.66]#,1.]
# sigm_list = [.33,.5,.66]#,1.]
alph_list = [.33,.5,.66]#,1.]
sigm_list = [.33,.66]#,1.]
kapp_list = [1]
plat_list = [10,30,100]

Q_list = [20,100]
gamma_list = [1,5]
pm_list = [.33,.66]


Ac_list = [0] # default 0
Bc_list = [0] # default 0
As_list = [0] # default 0
Bs_list = [0] # default 0

alphac_list = [0] # default 0; >=0
alphas_list = [1.] # default 1

s_list = [1.] # default 1
sv_list = [1.] # default 1
sp_list = [1.] # default 1
sr_list = [1.] # default 1


outsuf = 's' if want_rsd else 'r'
outsim = sim_params['sim_name']
# outsim = outsim[outsim.find("_c0"):] + '_z%.1f_'%zbox
if want_zcv: outsuf += '_zcv'

outfn_base = 'boxes/z%.1f/ELG/'%zbox
if not os.path.exists(outfn_base):
    os.makedirs(outfn_base)


maxobj  = round(1e7) # Should be an integer.

hodkeys = ['logM_cut','logM1','sigma','kappa','alpha',\
           # default values
           # 0 0 0 0
           'Q','gamma','p_max'] 
kk = 0

for lgMcut in lgMc_list:
  for plateau in plat_list:
    lgM1 = lgMcut + np.log10(plateau)
    for alph in alph_list:
      for sigm in sigm_list:
        for kappa in kapp_list:
          for Q in Q_list:
            for gamma in gamma_list:
                for pmax in pm_list:
                    vec = ([lgMcut,lgM1, sigm,kappa,alph,\
                        Q,gamma,pmax])
                      
                    for key,val in zip(hodkeys,vec):
                      HOD_params['ELG_params'][key] = val            
                    hod      = [HOD_params['ELG_params'][k] for k in hodkeys]
                    print(HOD_params)
                    tmp = {}
                    for key in hodkeys: tmp[key] = HOD_params['ELG_params'][key]
                    hod_string = get_hod_string(tmp)
                    
                    outfn = outfn_base + f'{hod_string}.asdf'
                    if os.path.exists(outfn): continue

                    newBall  = AbacusHOD(sim_params,HOD_params,clustering_params)
                    mock_dict= newBall.run_hod(newBall.tracers,want_rsd,\
                                               write_to_disk,Nthread=16)

                    # downsample for ELGs due to large catalog size
                    nobj  = mock_dict['ELG']['mass'].size
                    if nobj>maxobj:
                        print("Have nobj=",nobj," downsampling to ",maxobj,flush=True)
                        rng  = np.random.default_rng()
                        inds = rng.choice(nobj,size=maxobj,replace=False)
                        for k in ['x','y','z','vx','vy','vz','mass','id']:
                            mock_dict['ELG'][k] = mock_dict['ELG'][k][inds]
                    
                    af = asdf.AsdfFile(mock_dict)
                    af.write_to(outfn) 
                    af.close()
                
                    del newBall, mock_dict

