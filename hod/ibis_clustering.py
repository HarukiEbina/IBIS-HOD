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


import os,sys
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
print(sim_params)
# Get the metaparameters for the simulation.
meta = get_meta(sim_params['sim_name'],redshift=sim_params['z_mock'])
#
# additional parameter choices
want_zcv      = True
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
    std_hod_list.append( {'logM_cut': 11.25, 'logM1': 11.948970004336019, 'sigma': 0.66, 'kappa': 1, 'alpha': 0.33})
    # bright
    std_hod_list.append({'logM_cut': 11.0, 'logM1': 12.0, 'sigma': 0.66, 'kappa': 1, 'alpha': 0.5})
    
elif zbox==3.:
    # z=3.0
    # total
    std_hod_list.append({'logM_cut': 12.0, 'logM1': 12.698970004336019, 'sigma': 0.66, 'kappa': 1, 'alpha': 0.5})
    # bright
    std_hod_list.append({'logM_cut': 11.75, 'logM1': 12.448970004336019, 'sigma': 0.66, 'kappa': 1, 'alpha': 0.33})


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
outsim = outsim[outsim.find("_c0"):] + '_z%.1f_'%zbox
if want_zcv: outsuf += '_zcv'

outfn = "lae_clustering"+outsim+outsuf+".json"

if os.path.exists(outfn):
    with open(outfn, 'r') as file:
        d = json.load(file)
    
    dats = d['mocks']
    vecs = d['mock_vecs']
else:
    dats, vecs = [], []


# maxobj  = 1e21 # Should be an integer.
maxobj  = round(1e7) # Should be an integer.
# maxobj  = 8000000 # Should be an integer.
# maxobj  = 4000000 # Should be an integer.
# maxobj = round(1e-4 * 1e9) # just for downsampling
hodkeys = ['logM_cut','logM1','sigma','kappa','alpha',\
           # default values
           # 0 0 0 0
           'Acent','Bcent','Asat','Bsat',\
           # 0 1
           'alpha_c','alpha_s',\
           # 1 1 1 1
           's','s_v','s_p','s_r']
kk = 0

for std_dict in std_hod_list:
  lgMcut = std_dict['logM_cut']
  lgM1 = std_dict['logM1']
  sigm = std_dict['sigma']
  kappa = std_dict['kappa']
  alph = std_dict['alpha']
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
                        # if vec in vecs: continue
                      
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
                        print(mock_dict['LRG'].keys())
                        plt.hist(np.log10(mock_dict['LRG']['mass']),60,alpha=.5)
                        plt.savefig('plots/tmp.png')
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
