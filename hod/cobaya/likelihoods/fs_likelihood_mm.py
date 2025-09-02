import numpy as np

from cobaya.theory     import Theory
from cobaya.likelihood import Likelihood
from scipy.interpolate import InterpolatedUnivariateSpline as Spline

from linear_theory import *

# from velocileptors.LPT.lpt_rsd_fftw          import LPT_RSD
# from velocileptors.LPT.moment_expansion_fftw import MomentExpansion
from velocileptors.Utils.spherical_bessel_transform import SphericalBesselTransform as SBT
from velocileptors.EPT.ept_fullresum_fftw import REPT

from scipy.interpolate import interp1d
import os
import json

class FSLikelihood(Likelihood):
    """
    full-shape likelihood to fit for P(k,mu)
    assumes diagonal covariance

    This likelihood is for the matter-matter auto spectrum, so it's limited in nuisance parameters (only alpha0)
    """
    zfid: float
    # Hz_fid: float
    # chiz_fid: float
    
    datfn: str
    
    def initialize(self):
        """Sets up the class."""
        self.loadData()
        os.system('echo here')
        #

    def get_requirements(self):
        req = {'pt_pk_ell_mod': None,\
               'H0': None,\
               'sigma8': None,
               # 'b1': None,\
               # 'b2': None,\
               # 'bs': None,\
               'alpha0_mm': None,\
              }
        return(req)
    
    def logp(self,**params_values):
        """Return a log-likelihood."""
        
        fs_thy  = self.fs_predict()
        # fs_obs  = self.fs_observe(fs_thy)

        # obs = fs_obs
        obs = fs_thy
        MU = self.MU
        
        # obs = np.array([interp1d(obs[:,0],obs[:,i+1])(self.kvec) for i in range(3)])
        obs = np.array([interp1d(obs[:,0],obs[:,i+1],fill_value='extrapolate')(self.kvec) for i in range(3)])
        obs = np.repeat(obs[0],self.Nmu)+0.5*(3*MU**2-1)*np.repeat(obs[1],self.Nmu)#+0.125*(35*MU**4-30*MU**2+3)*np.repeat(obs[2],self.Nmu)
        
        self.last_thy = obs
        chi2 = np.einsum('i,i,i',self.dd-obs,self.cinv,self.dd-obs)
        self.last_chi2 = chi2
        return (-0.5*chi2)
        #
    def close(self):
        return 
        
        
    def loadData(self):
        """
        Loads the required data.
        
        Do this in two steps... first load full shape data then xirecon, concatenate after.
        
        The covariance is assumed to already be joint in the concatenated format.
        
        """
        
        ff = open(self.datfn)
        simdata = json.load(ff)
        dd = simdata['d_pkmu']
        cinv = simdata['Cinv']
        self.kvec = np.array(simdata['k'])
        Nmu=50
        k = np.repeat(self.kvec,Nmu)
        
        self.dd = np.array(dd)
        self.cinv = np.array(cinv)
        
        self.Nk = len(self.kvec)
        self.Nmu = 50
        mu = np.linspace(0.,1.,self.Nmu)
        MU = np.tile(mu,self.Nk)
        self.MU = MU                

    def fs_predict(self):
        """Use the PT model to compute P_ell, given biases etc."""
        
        pp   = self.provider
        modPT= pp.get_result('pt_pk_ell_mod')
        hub = pp.get_param('H0') / 100.
        sig8 = pp.get_param('sigma8')
        OmM = pp.get_param('omegam')
        hub  = pp.get_Hubble(0)[0]/100.
        omnuh2 = pp.get_param('omnuh2')
        fnu =  omnuh2/ hub**2 /OmM
        zfid = self.zfid
        ff   = f_of_a(1/(1.+zfid), OmegaM=pp.get_param('omegam')) * (1 - 0.6 * fnu)
        #
        # b1   = pp.get_param('b1') - 1.0
        # # b1   = pp.get_param('bsig8')/sig8 - 1.0
        # b2   = pp.get_param('b2')
        # bs   = pp.get_param('bs')
        alp0 = pp.get_param('alpha0_mm')
        alp2 = 0
        sn0  = 0
        sn2  = 0
        b1   = 1
        b2   = 0
        bs   = 0
        
        bias = [b1, b2, bs, 0.]
        cterm = [alp0,alp2,0,0]
        stoch = [sn0, sn2, 0]
        bvec = bias + cterm + stoch
        
        # bvec in EPT space
        kv, p0,p2,p4 = modPT.compute_redshift_space_power_multipoles(bvec,ff)
        # kv, p0, p2, p4 = modPT.combine_bias_terms_pkell(bvec)

        # p0 = np.repeat(p0,self.Nmu) 
        # p2 = np.repeat(p2,self.Nmu) 
        # p4 = np.repeat(p4,self.Nmu)
        # pkmu = p0+0.5*(3*MU**2-1)*p2+0.125*(35*MU**4-30*MU**2+3)*p4
        # return pkmu
        
        # # Put a point at k=0 to anchor the low-k part of the Spline.
        kv,p0 = np.append([0.,],kv),np.append([0.0,],p0)
        p2,p4 = np.append([0.,],p2),np.append([0.0,],p4)
        tt    = np.array([kv,p0,p2,p4]).T
        return(tt)
        #
        
#     def fs_observe(self,tt):
#         """Apply the window function matrix to get the binned prediction."""
        
#         # Have to stack ell=0, 2 & 4 in bins of 0.001h/Mpc from 0-0.4h/Mpc.
#         kv  = np.linspace(0.0,0.4,400,endpoint=False) + 0.0005
#         thy =                     Spline(tt[:,0],tt[:,1],ext=3)(kv)
#         thy = np.concatenate([thy,Spline(tt[:,0],tt[:,2],ext=3)(kv)])
#         thy = np.concatenate([thy,Spline(tt[:,0],tt[:,3],ext=3)(kv)])
        
#         # wide angle
#         expanded_model = np.matmul(self.matM, thy )
#         # Convolve with window (true) −> (conv) see eq. 2.18
#         convolved_model = np.matmul(self.matW, expanded_model )
#         # keep only the monopole and quadrupole
#         convolved_model = convolved_model[self.fitiis]
    
#         return convolved_model
    

class PT_pk_theory(Theory):
    """A class to return a PT P_ell module."""
    # From yaml file.
    zfid:     float
    # chiz_fid: float
    # Hz_fid:   float
    # datfn: str
    #
    def initialize(self):
        """Sets up the class."""
        # Don't need to do anything.
        # self.kvec = np.loadtxt(datfn)[0]
        # pass
    def get_requirements(self):
        """What we need in order to provide P_ell."""
        zg  = np.linspace(0,self.zfid,8,endpoint=True)
        # Don't need sigma8_z, fsigma8 or radial distance
        # here, but want them up in likelihood and they
        # only depend on cosmological things (not biases).
        req = {\
               'omegam': None,\
               'omnuh2': None,\
               'H0': None,\
               'Pk_interpolator': {'k_max': 30,'z': zg,\
                                   'nonlinear': False,\
                                   'vars_pairs': [['delta_nonu','delta_nonu']]},\
               'Hubble':   {'z': [0.0,self.zfid]},\
               'sigma8_z': {'z': [0.0,self.zfid]},\
               'fsigma8':  {'z': [self.zfid]},\
               'comoving_radial_distance': {'z': [self.zfid]}\
              }
        return(req)
    def get_can_provide(self):
        """What do we provide: a PT class that can compute xi_ell."""
        return ['pt_pk_ell_mod']
    
    def calculate(self, state, want_derived=True, **params_values_dict):
        """Create and initialize the PT class."""
        # Make shorter names.
        pp   = self.provider
        zfid = self.zfid
        
        # Get cosmological parameters
        OmM = pp.get_param('omegam')
        hub  = pp.get_Hubble(0)[0]/100.
        omnuh2 = pp.get_param('omnuh2')
        fnu =  omnuh2/ hub**2 /OmM
        
        
        ff   = f_of_a(1/(1.+zfid), OmegaM=pp.get_param('omegam')) * (1 - 0.6 * fnu)
        
        kmin = 5e-4
        kmax = 1.
        # kmax = 5.
        # kmin = .005
        # kmax = .5
        Nk = 500
        # Nk = 1000
        ki = np.logspace(np.log10(kmin),np.log10(kmax),Nk)
        # ki   = np.logspace(-3.0,1.0,200)
        pi   = pp.get_Pk_interpolator(nonlinear=False,var_pair=['delta_nonu','delta_nonu'])
        pi   = pi.P(self.zfid,ki*hub)*hub**3
        
        # # Work out the A-P scaling to the fiducial cosmology.        
        # Hz   = pp.get_Hubble(self.zfid)[0]/pp.get_Hubble(0)[0]
        # chiz = pp.get_comoving_radial_distance(self.zfid)[0]*hub
        # apar,aperp = self.Hz_fid/Hz,chiz/self.chiz_fid
        
        ept = REPT(ki,pi,rbao=110,sbao=None,)
        ept.compute_redshift_space_power_multipoles_tables(ff)
        state['pt_pk_ell_mod'] = ept
        #