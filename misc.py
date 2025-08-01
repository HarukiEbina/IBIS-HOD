from classy import Class
import numpy as np

from scipy.interpolate import interp1d


def get_cosmo_P18():
    params = {'output': 'tCl lCl mPk',
              'l_max_scalars': 2000,
              'lensing': 'yes',
              'P_k_max_h/Mpc': 50.,
              'non linear':'halofit', 
              'z_pk': '0.0,1087',
              'A_s': 2.10732e-9,
              'n_s': 0.96824,
              'alpha_s': 0.,
              'h': 0.6736,
              'N_ur': 2.0328,
              'N_ncdm': 1,
              'm_ncdm': 0.06,
              # 'm_ncdm': '0.01,0.05',
              'tau_reio': 0.0544,
              'omega_b': 0.02237,
              'omega_cdm': 0.1200,
              'Omega_k': 0.}

    cosmo = Class() 
    cosmo.set(params) 
    cosmo.compute() 
    return cosmo


def angular_distance(ra1_arr, dec1_arr, ra2_arr, dec2_arr):
    # Convert degrees to radians
    ra1_arr = np.radians(ra1_arr)
    dec1_arr = np.radians(dec1_arr)
    ra2_arr = np.radians(ra2_arr)
    dec2_arr = np.radians(dec2_arr)
    
    # Spherical law of cosines
    cos_theta = np.sin(dec1_arr) * np.sin(dec2_arr) + np.cos(dec1_arr) * np.cos(dec2_arr) * np.cos(ra1_arr - ra2_arr)
    
    # Compute the angle in radians
    theta = np.arccos(cos_theta)
    
    # Convert radians back to degrees
    theta_degrees = np.degrees(theta)
    
    return theta_degrees

def n2N(n_func,zedges,fsky,cosmo):
    h = cosmo.get_current_derived_parameters(['h'])['h']
    N_arr = np.zeros(len(zedges)-1)
    for i in range(len(zedges)-1):
        # idxs = (zz>zedges[i]) & (zz<zedges[i+1])
        zz = np.linspace(zedges[i],zedges[i+1],200)
        nz = n_func(zz)
        chi_zz = np.array([cosmo.comoving_distance(z) for z in zz])*h
        N_arr[i] += np.trapz(nz*chi_zz**2,chi_zz)
    N_arr *= 4*np.pi*fsky
    return N_arr

def N2n(N,zedges,fsky,cosmo):
    h = cosmo.get_current_derived_parameters(['h'])['h']
    V_bin = np.zeros(len(zedges)-1)
    chi = lambda z: cosmo.comoving_distance(z)*h
    for i in range(len(zedges)-1):
        V_bin[i] += 4*np.pi*fsky/3 * (chi(zedges[i+1])**3 - chi(zedges[i])**3)
    n_bin = N/V_bin
    def n_func(z):
        for i in range(len(zedges)-1):
            if z>=zedges[i] and z<=zedges[i+1]:
                return n_bin[i]
    zz = np.linspace(zedges[0],zedges[-1],200*len(zedges))
    n_bin_fine = np.array([n_func(z) for z in zz])
    return interp1d(zz,n_bin_fine,fill_value=0.,bounds_error=False)
    # return interp1d(zz,n_bin_fine,fill_value='extrapolate')


def pick_selection(ss,d):
    d.meta['NAMES'] = 'AR_LAE_M411,AR_LAE_M438,AR_LBG_M438,AR_LAE_M464,AR_LBG_M464,AR_LAE_M490,AR_LBG_M490,AR_LAE_M517,AR_LBG_M517,DJS_LAE_M411,DJS_LBG_M411,DJS_LAE_M438,DJS_LBG_M438,DJS_LAE_M464,DJS_LBG_M464,DJS_LAE_M490,DJS_LBG_M490,DJS_LAE_M517,DJS_LAE_6,HE_M411,HE_M438,HE_M464,HE_M490,HE_M517'
    d.meta['BITS'] = '0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23'
    d['IBIS_TARGET'] = [int(d['IBIS_TARGET'][i]) for i in range(len(d))]
    names = (d.meta["NAMES"].split(","))
    bits = np.array(d.meta["BITS"].split(",")).astype(int)
    
    bit = bits[names.index(ss)]
    return d[(d['IBIS_TARGET'] & 2**bit)>0]
    
def selection_idx(ss,d):
    d.meta['NAMES'] = 'AR_LAE_M411,AR_LAE_M438,AR_LBG_M438,AR_LAE_M464,AR_LBG_M464,AR_LAE_M490,AR_LBG_M490,AR_LAE_M517,AR_LBG_M517,DJS_LAE_M411,DJS_LBG_M411,DJS_LAE_M438,DJS_LBG_M438,DJS_LAE_M464,DJS_LBG_M464,DJS_LAE_M490,DJS_LBG_M490,DJS_LAE_M517,DJS_LAE_6,HE_M411,HE_M438,HE_M464,HE_M490,HE_M517'
    d.meta['BITS'] = '0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23'
    d['IBIS_TARGET'] = [int(d['IBIS_TARGET'][i]) for i in range(len(d))]
    names = (d.meta["NAMES"].split(","))
    bits = np.array(d.meta["BITS"].split(",")).astype(int)

    bit = bits[names.index(ss)]
    return (d['IBIS_TARGET'] & 2**bit)>0

def selection_idx_M(mm,d):
    d.meta['NAMES'] = 'AR_LAE_M411,AR_LAE_M438,AR_LBG_M438,AR_LAE_M464,AR_LBG_M464,AR_LAE_M490,AR_LBG_M490,AR_LAE_M517,AR_LBG_M517,DJS_LAE_M411,DJS_LBG_M411,DJS_LAE_M438,DJS_LBG_M438,DJS_LAE_M464,DJS_LBG_M464,DJS_LAE_M490,DJS_LBG_M490,DJS_LAE_M517,DJS_LAE_6,HE_M411,HE_M438,HE_M464,HE_M490,HE_M517'
    d.meta['BITS'] = '0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23'
    d['IBIS_TARGET'] = [int(d['IBIS_TARGET'][i]) for i in range(len(d))]
    names = (d.meta["NAMES"].split(","))
    bits = np.array(d.meta["BITS"].split(",")).astype(int)

    idx = np.zeros(len(d),dtype = bool)
    for sb in ['AR_LAE_','AR_LBG_','HE_','DJS_LBG_','DJS_LAE_']:
        ss = sb+mm
        try: 
            bit = bits[names.index(ss)]
            # print(ss)
        except: bit = -1
        if bit>-1:idx |= selection_idx(ss,d)
    return idx

def selection_idx_P(s_base,d):
    d.meta['NAMES'] = 'AR_LAE_M411,AR_LAE_M438,AR_LBG_M438,AR_LAE_M464,AR_LBG_M464,AR_LAE_M490,AR_LBG_M490,AR_LAE_M517,AR_LBG_M517,DJS_LAE_M411,DJS_LBG_M411,DJS_LAE_M438,DJS_LBG_M438,DJS_LAE_M464,DJS_LBG_M464,DJS_LAE_M490,DJS_LBG_M490,DJS_LAE_M517,DJS_LAE_6,HE_M411,HE_M438,HE_M464,HE_M490,HE_M517'
    d.meta['BITS'] = '0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23'
    d['IBIS_TARGET'] = [int(d['IBIS_TARGET'][i]) for i in range(len(d))]
    names = (d.meta["NAMES"].split(","))
    bits = np.array(d.meta["BITS"].split(",")).astype(int)
    
    idx = np.zeros(len(d),dtype = bool)
    for band in (m_bands+ ['6']):
        ss = s_base + '_' + band
        try: 
            bit = bits[names.index(ss)]
            # print(ss)
        except: bit = -1
        if bit>-1:idx |= selection_idx(ss,d)
    return idx