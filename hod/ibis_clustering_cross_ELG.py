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


import os
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

Script to generate the gg, gm spectra with ZCV - still in progress
Missing capability for zcv on cross-spectrum, but have zcv for auto-spectra

mm spectra is generated in another script (since it only needs to be generated once)

'''

#
# Load the config file and parse in relevant parameters
path2config= './lae_ext.yaml'
config     = yaml.safe_load(open(path2config))

zbox = 3.0
config['sim_params']['z_mock'] = zbox

config['HOD_params']['tracer_flags']['LRG'] = False
config['HOD_params']['tracer_flags']['ELG'] = True

sim_params = config['sim_params']
HOD_params = config['HOD_params']
clustering_params = config['clustering_params']
print(sim_params)

# Get the metaparameters for the simulation.
meta = get_meta(sim_params['sim_name'],redshift=sim_params['z_mock'])
n_chunks = 34
#
# additional parameter choices
want_zcv      = False
want_rsd      = HOD_params['want_rsd']  & False
write_to_disk = HOD_params['write_to_disk']
#
# Load the rp pi binning from the config file. Note pimax and pi_bin_size are ints.
bin_params = clustering_params['bin_params']
rpbins     = np.logspace(bin_params['logmin'],\
                         bin_params['logmax'],bin_params['nbins']+1)
pimax,dpi  = clustering_params['pimax'],clustering_params['pi_bin_size']
# work out the rp and pi bin centers (assume log binning as above)
Rcen = np.sqrt(rpbins[1:]*rpbins[:-1])
Zcen = np.arange(0.0,float(pimax),dpi) + 0.5*dpi
#
# Now loop over HODs writing out HOD, wp(R), xi0, etc. for each.
#
std_hod_list = []

if zbox==2.5:
    # z=2.5 
    # total
    print('error')
    assert True == False
    
elif zbox==3.:
    # z=3.0
    # total
    std_hod_list.append( {'logM_cut': 11.0,
                          'logM1': 13.0,
                          'sigma': 0.33,
                          'kappa': 1,
                          'alpha': 0.66,
                          'Q': 100,
                          'gamma': 5,
                          'p_max': 0.66}) 
    # bright
    # std_hod_list.append( {'logM_cut': 11.0,
    #                       'logM1': 12.477121254719663,
    #                       'sigma': 0.33,
    #                       'kappa': 1,
    #                       'alpha': 0.33,
    #                       'Q': 20,
    #                       'gamma': 5,
    #                       'p_max': 0.66})
    
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
#
outsuf = 's' if want_rsd else 'r'
outsim = sim_params['sim_name']
outsim = outsim[outsim.find("_c0"):] + '_z%.1f_'%zbox
if want_zcv: outsuf += '_zcv'

outfn = "lae_clustering_ELG"+outsim+outsuf+"_gm.json"

if os.path.exists(outfn):
    with open(outfn, 'r') as file:
        d = json.load(file)
    
    dats = d['mocks']
    vecs = d['mock_vecs']
else:
    dats, vecs = [], []

maxobj  = round(1e7) # Should be an integer.
# maxobj  = 8000000 # Should be an integer.
# maxobj  = 4000000 # Should be an integer.
hodkeys = ['logM_cut','logM1','sigma','kappa','alpha',\
           # default values
           # 0 0 0 0
           'Q','gamma','p_max'] 
kk = 0

# sampfact = 10
# mpos,vcel = load_matter(sim_params['sim_dir']+sim_params['sim_name'],\
#               meta['BoxSize'],sim_params['z_mock'],n_chunks,\
#               meta['NP'],sampfact,type_AB='A')
# matter_dict['matter'] = {}
# matter_dict['matter']['x'] = mpos[:,0]
# matter_dict['matter']['y'] = mpos[:,1]
# matter_dict['matter']['z'] = mpos[:,2]
sim_name = sim_params['sim_name']
d = Table.read('matter_{:s}_z{:.1f}.fits'.format(sim_name,zbox))
matter_dict = {}
matter_dict['matter'] = {}
for key in d.keys():
    matter_dict['matter'][key] = d[key]
del d

ndm = matter_dict['matter']['x'].size

for std_dict in std_hod_list:
    lgMcut = std_dict['logM_cut']
    lgM1 = std_dict['logM1']
    sigm = std_dict['sigma']
    kappa = std_dict['kappa']
    alph = std_dict['alpha']
    Q = std_dict['Q']
    gamma = std_dict['gamma']
    pmax = std_dict['p_max']
    vec = ([lgMcut,lgM1, sigm,kappa,alph,\
        Q,gamma,pmax])

    ### check for duplicates in vecs
    if vec in vecs: continue
      
    for key,val in zip(hodkeys,vec):
      HOD_params['ELG_params'][key] = val            
    hod      = [HOD_params['ELG_params'][k] for k in hodkeys]
    print(HOD_params)
    newBall  = AbacusHOD(sim_params,HOD_params,clustering_params)
    mock_dict= newBall.run_hod(newBall.tracers,want_rsd,\
                               write_to_disk,Nthread=16)

    nobj  = mock_dict['ELG']['mass'].size
    ncen  = mock_dict['ELG']['Ncent']
    fsat  = 1-float(ncen)/float(nobj)
    print(nobj,ncen,fsat)

    #
    if nobj>maxobj:
        print("Have nobj=",nobj," downsampling to ",maxobj,flush=True)
        rng  = np.random.default_rng()
        inds = rng.choice(nobj,size=maxobj,replace=False)
        for k in ['x','y','z','vx','vy','vz','mass','id']:
            mock_dict['ELG'][k] = mock_dict['ELG'][k][inds]
    xiell = newBall.compute_multipole(mock_dict,rpbins,pimax,\
                      sbins=rpbins,nbins_mu=11,orders=[0,2,4])['ELG_ELG']
    wpR   = xiell[0*len(Rcen):1*len(Rcen)]
    xi0   = xiell[1*len(Rcen):2*len(Rcen)]
    xi2   = xiell[2*len(Rcen):3*len(Rcen)]
    xi4   = xiell[3*len(Rcen):4*len(Rcen)]
    bb    = 0.0
    #
    mock_dict['matter'] = matter_dict['matter']

    if want_zcv:
        # Compute variance reduced spectra.
        zcv_dict = newBall.apply_zcv(mock_dict,config)
        kk = zcv_dict['k_binc']
        pkl= zcv_dict['Pk_tr_tr_ell_zcv'] # variance-reduced multipoles
        pk0= pkl[0]
        pk2= pkl[1]
        bb = zcv_dict['bias'][1] + 1.0 # Convert to Eulerian bias.
    else:
        pk3d = newBall.compute_power(mock_dict,nbins_k=50,nbins_mu=11,\
                 k_hMpc_max=0.5,logk=False,poles=[0,2],num_cells=512,\
                 paste='TSC',compensated=True,interlaced=True)
        # pk3d = newBall.compute_power(mock_dict,nbins_k=20,nbins_mu=3,\
        #          k_hMpc_max=0.5,logk=False,poles=[0,2],num_cells=512,\
        #          paste='TSC',compensated=True,interlaced=True)
        print(pk3d.keys())
        kk   = pk3d['k_binc']
        pkl  = pk3d['ELG_matter_ell']
        pk0  = pkl[:,0]
        pk2  = pkl[:,1]
        bb   = 0.0 # Placeholder.
    #
    dats.append({'hod':hod,'nobj':nobj,'fsat':fsat,'bias':bb,\
                 'wp':wpR.tolist(),\
                 'xi0':xi0.tolist(),'xi2':xi2.tolist(),'xi4':xi4.tolist(),\
                 'pk0':pk0.tolist(),'pk2':pk2.tolist()})
    vecs.append(vec)
    del newBall, mock_dict


# Now write out our answers.
outdict = {}
for k in simkeys: outdict[k]=meta[k]
outdict['in_red' ] = want_rsd
outdict['R'      ] = Rcen.tolist()
outdict['k'      ] = kk.tolist()
outdict['hodkeys'] = hodkeys
outdict['mocks'  ] = dats
outdict['mock_vecs'  ] = vecs
with open(outfn,"w") as fout:
    json.dump(outdict,fout,indent=2)
