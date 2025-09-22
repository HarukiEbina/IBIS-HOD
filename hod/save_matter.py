#!/usr/bin/env python3
#
import numpy as np
import yaml
import json
import gc
import random

from abacusnbody.data.bitpacked import unpack_rvint

import asdf
from abacusnbody.hod.abacus_hod import AbacusHOD
from abacusnbody.metadata       import get_meta
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
from scipy.interpolate import interp1d
from astropy.table import Table, vstack

from matplotlib import pyplot as plt


import os, sys
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



def get_metadata(simdir,simname,redshift,proxy):
    """This routine gets the metadata for the simulation.  For simulations
    in the Abacus database this duplicates the get_meta method, but that
    database is not updated frequently and for some simulations some keys
    are missing.
    Proxy is the name of a simulation with the same cosmology that has all
    of the keys, used as a work-around for missing data."""
    # First fill in the "missing" information from the proxy simulation.
    dd = {}
    pr = get_meta(proxy,redshift)
    for k in ['n_s','omega_b','omega_cdm','omega_ncdm',\
              'N_ncdm','N_ur','SimSet']:
        dd[k] = pr[k]
    # Now get the remaining information from the ASDF files directly.
    fn = simdir+simname+\
         '/halos/z{:.3f}/halo_info/halo_info_000.asdf'.format(float(redshift))
    with asdf.open(fn) as af:
        meta = af['header']
        for k in simkeys:
            if k in meta: dd[k] = meta[k]
    return(dd)
    #
    
def load_matter(sim_dir,Lbox,z,n_chunks,N_parts,f_down,type_AB='A'):
    """Loads the matter particles from the sub-sample files."""
    # estimate total number of particles and preselect indices
    N_all = 0
    N_offset = np.zeros(n_chunks, dtype=int)
    N_file = np.zeros(n_chunks, dtype=int)
    for i_chunk in range(n_chunks):
        print(i_chunk, n_chunks)
        # halo and field particles
        fn_halo = sim_dir+f'/halos/z{z:.3f}/halo_rv_{type_AB}/halo_rv_{type_AB}_{i_chunk:03d}.asdf'
        fn_field = sim_dir+f'/halos/z{z:.3f}/field_rv_{type_AB}/field_rv_{type_AB}_{i_chunk:03d}.asdf'
        N_this = asdf.open(fn_halo)['data']['rvint'].shape[0]+asdf.open(fn_field)['data']['rvint'].shape[0]
        N_offset[i_chunk] = N_all
        N_file[i_chunk] = N_this
        N_all += N_this
    gc.collect()
    #print("offsets", N_offset, flush=True)
    #print("per file", N_file, flush=True)
    print("Total particles in halo and field files: ",N_all,flush=True)

    # global indices to keep
    inds_keep = random.sample(range(N_all), N_all//f_down)
    N_keep = len(inds_keep)
    print("N_keep", N_keep, flush=True)
    pos_down = np.zeros((N_keep, 3), dtype=np.float32)
    vel_down = np.zeros((N_keep, 3), dtype=np.float32)

    # load the matter particles
    count = 0
    for i_chunk in range(n_chunks):
        print("Reading ",i_chunk," of ",n_chunks,flush=True)
        # indices to keep in this chunk
        inds_keep_this = inds_keep - N_offset[i_chunk]
        inds_keep_this = inds_keep_this[(inds_keep_this >= 0) & (inds_keep_this < N_file[i_chunk])]

        # halo and field particles
        fn_halo = sim_dir+f'/halos/z{z:.3f}/halo_rv_{type_AB}/halo_rv_{type_AB}_{i_chunk:03d}.asdf'
        fn_field = sim_dir+f'/halos/z{z:.3f}/field_rv_{type_AB}/field_rv_{type_AB}_{i_chunk:03d}.asdf'

        halo_data = (asdf.open(fn_halo)['data'])['rvint']
        pos_halo, vel_halo = unpack_rvint(halo_data, Lbox, float_dtype=np.float32, velout=None)
        #print("pos_halo = ", pos_halo[:5],flush=True)

        field_data = (asdf.open(fn_field)['data'])['rvint']
        pos_field, vel_field = unpack_rvint(field_data, Lbox, float_dtype=np.float32, velout=None)
        #print("pos_field = ", pos_field[:5],flush=True)

        # stack halo and field particles
        pos_both = np.vstack((pos_halo, pos_field))
        vel_both = np.vstack((vel_halo, vel_field))

        pos_down[count:count+len(inds_keep_this)] = pos_both[inds_keep_this]
        vel_down[count:count+len(inds_keep_this)] = vel_both[inds_keep_this]
        count += len(inds_keep_this)
        del halo_data, pos_halo, vel_halo, field_data, pos_field, vel_field
        gc.collect()
    print("These two must be the same: ",count,pos_down.shape[0],flush=True)
    pos_down = pos_down[:count]
    vel_down = vel_down[:count]
    return pos_down, vel_down


'''
Script to save the matter particles for a simulation
'''

#
# Load the config file and parse in relevant parameters
path2config= './lae_ext.yaml'
config     = yaml.safe_load(open(path2config))

zbox = float(sys.argv[1]) # 2.5
config['sim_params']['z_mock'] = zbox

sim_params = config['sim_params']
HOD_params = config['HOD_params']
clustering_params = config['clustering_params']
#
# Get the metaparameters for the simulation.
meta = get_meta(sim_params['sim_name'],redshift=sim_params['z_mock'])
n_chunks = 34
#
sim_name = sim_params['sim_name']

# sampfact = 10
sampfact = 10
mpos,vcel = load_matter(sim_params['sim_dir']+sim_params['sim_name'],\
              meta['BoxSize'],sim_params['z_mock'],n_chunks,\
              meta['NP'],sampfact,type_AB='A')
# matter_dict['matter'] = {}
# matter_dict['matter']['x'] = mpos[:,0]
# matter_dict['matter']['y'] = mpos[:,1]
# matter_dict['matter']['z'] = mpos[:,2]
# ndm = matter_dict['matter']['x'].size
matter_dict = {}
matter_dict['x'] = mpos[:,0]
matter_dict['y'] = mpos[:,1]
matter_dict['z'] = mpos[:,2]

dat_table = Table()
for key in matter_dict.keys():
    dat_table[key] = matter_dict[key]
dat_table.write('matter_{:s}_z{:.1f}.fits'.format(sim_name,zbox),overwrite=True)
del matter_dict, dat_table

