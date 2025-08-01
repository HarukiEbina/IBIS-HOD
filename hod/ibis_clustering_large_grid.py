#!/usr/bin/env python3
#
import numpy as np
import yaml
import json

import asdf
from abacusnbody.hod.abacus_hod import AbacusHOD
from abacusnbody.metadata       import get_meta
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
from scipy.interpolate import interp1d

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
print(zbox,sim_params)
# Get the metaparameters for the simulation.
meta = get_meta(sim_params['sim_name'],redshift=sim_params['z_mock'])
#
# additional parameter choices
want_zcv      = False
want_rsd      = HOD_params['want_rsd'] & False
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
lgMc_list = [11.00,11.25,11.50,11.75,12.0,12.25,12.5]
alph_list = [.33,.5,.66]#,1.]
sigm_list = [.33,.5,.66]#,1.]
kapp_list = [1]
plat_list = [5,10]

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

outfn = "large_grid"+outsim+outsuf+".json"
print(outfn)
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
           'Acent','Bcent','Asat','Bsat',\
           # 0 1
           'alpha_c','alpha_s',\
           # 1 1 1 1
           's','s_v','s_p','s_r']
kk = 0

for lgMcut in lgMc_list:
  for alph in alph_list:
    for sigm in sigm_list:
      for kappa in kapp_list:
        for plateau in plat_list:
          lgM1 = lgMcut + np.log10(plateau)
          for Acent in Ac_list:
            for Bcent in Bc_list:
              for Asat in As_list:
                for Bsat in Bs_list:
                  for alphac in alphac_list:
                    for alphas in alphas_list:
                      for s in s_list:
                        for s_v in sv_list:
                          for s_p in sp_list:
                            for s_r in sr_list:
                                vec = ([lgMcut,lgM1, sigm,kappa,alph,\
                                        Acent,Bcent,Asat,Bsat,\
                                        alphac,alphas,\
                                        s,s_v,s_p,s_r])
        
                                ### check for duplicates in vecs
                                if vec in vecs: continue
                              
                                for key,val in zip(hodkeys,vec):
                                  HOD_params['LRG_params'][key] = val            
                                hod      = [HOD_params['LRG_params'][k] for k in hodkeys]
                                print(HOD_params)
                                newBall  = AbacusHOD(sim_params,HOD_params,clustering_params)
                                mock_dict= newBall.run_hod(newBall.tracers,want_rsd,\
                                                           write_to_disk,Nthread=16)

                                nobj  = mock_dict['LRG']['mass'].size
                                ncen  = mock_dict['LRG']['Ncent']
                                fsat  = 1-float(ncen)/float(nobj)
                                print(nobj,ncen,fsat)

                                # The first ``'Ncent'`` galaxies in the catalog are always centrals and the rest are satellites.

                                # cen_vel = mock_dict['LRG']['vz'][:ncen]
                                # sat_vel = mock_dict['LRG']['vz'][ncen:]
                                # common_ids, idx1, idx2 = np.intersect1d(mock_dict['LRG']['id'][:ncen], mock_dict['LRG']['id'][ncen:], return_indices=True)
                                # delta_vz = sat_vel[idx2]-cen_vel[idx1]

                                # cen_vel = mock_dict['LRG']['vz'][:ncen]
                                # counter = 0
                                # sat_vel = mock_dict['LRG']['vz'][ncen:]
                                # sat_id = mock_dict['LRG']['id'][ncen:]
                                # delta_vz = None
                                # while True:
                                #     print(sat_vel.size)
                                #     _, indices = np.unique(sat_id, return_inverse=True)
                                #     common_ids, idx1, idx2 = np.intersect1d(mock_dict['LRG']['id'][:ncen], sat_id[indices], return_indices=True)
                                #     if delta_vz is None: delta_vz = sat_vel[indices][idx2]-cen_vel[idx1]
                                #     else: delta_vz = np.vstack(delta_vz,sat_vel[indices][idx2]-cen_vel[idx1])
                                #     print(sat_vel.size, indices.size,idx1.size)
                                #     if indices.size == sat_vel.size: break
                                #     sat_vel = np.delete(sat_vel,indices)
                                #     sat_id = np.delete(sat_vel,indices)
                                
                                # ids1 = mock_dict['LRG']['id'][:ncen]
                                # ids2 = mock_dict['LRG']['id'][ncen:]
                                # # Use broadcasting to compare all ids
                                # matches = ids1[:, None] == ids2  # shape (len(list1), len(list2))                                
                                # # Get indices where IDs match
                                # i1, i2 = np.where(matches)
                                # delta_vz = sat_vel[i2]-cen_vel[i1]
                                # print('out of {:d} satellite galaxies, {:d} found matches'.format(ids2.size,i2.size))
                                # hist, bin_edges = np.histogram(delta_vz, bins=np.linspace(-1000,1000,31))

                                # print('out of {:d} satellite galaxies, {:d} found matches'.format(mock_dict['LRG']['vz'][ncen:].size,delta_vz.size))
                                
                                # delta_vz = np.abs(delta_vz)
                                # hist, bin_edges = np.histogram(delta_vz, bins=np.linspace(0,1350,28))
                                #
                                if nobj>maxobj:
                                    print("Have nobj=",nobj," downsampling to ",maxobj,flush=True)
                                    rng  = np.random.default_rng()
                                    inds = rng.choice(nobj,size=maxobj,replace=False)
                                    for k in ['x','y','z','vx','vy','vz','mass','id']:
                                        mock_dict['LRG'][k] = mock_dict['LRG'][k][inds]
                                xiell = newBall.compute_multipole(mock_dict,rpbins,pimax,\
                                                  sbins=rpbins,nbins_mu=11,orders=[0,2,4])['LRG_LRG']
                                wpR   = xiell[0*len(Rcen):1*len(Rcen)]
                                xi0   = xiell[1*len(Rcen):2*len(Rcen)]
                                xi2   = xiell[2*len(Rcen):3*len(Rcen)]
                                xi4   = xiell[3*len(Rcen):4*len(Rcen)]
                                bb    = 0.0
                                #
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
                                    kk   = pk3d['k_binc']
                                    pkl  = pk3d['LRG_LRG_ell']
                                    pk0  = pkl[:,0]
                                    pk2  = pkl[:,1]
                                    bb   = 0.0 # Placeholder.
                                # hist, bin_edges = np.histogram(mock_dict['LRG']['vz'], bins=np.linspace(-1000,1000,31),density=True)
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
