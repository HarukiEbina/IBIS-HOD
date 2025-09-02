import numpy as np

import json
from   astropy.table        import Table, Column, vstack
# from   scipy.interpolate    import InterpolatedUnivariateSpline as Spline
from scipy.interpolate import make_smoothing_spline
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from classy import Class

import sys, os
sys.path.append('..')

from rotate_to import rotate_to
# from   fake_lbg import MockLBG
from   calc_xi  import calc_xi
from   calc_wt  import calc_wt, calc_wt_cross
from   make_mask import SurveyMask

from forecast_utils import *

# m_bands = ['M411','M438','M464','M490','M517']
band_z = {
    'M411':[2.26,2.56],
    'M438':[2.47,2.77],
    'M464':[2.68,2.98],
    'M490':[2.89,3.2],
    'M517':[3.1,3.41]
}
band_chimid = {'M411': 3945.5839397567484, 'M438': 4112.149560017418, 'M464': 4265.362496255113, 'M490': 4409.927293886262, 'M517': 4540.894419987146}
# chi chimid=chibar = (chi1+chi2)/2

class mc_forecast(object):    


    #########
    # SETUP #
    #########
    def __init__(self, name, hod_params, basedir='.',zbox=3,sim_name='AbacusSummit_high_c000_ph100'):
        # if there is config, don't go through the work of generating it again
        # input of only filepath assumes that it's already generated
        self.basedir = os.path.join(basedir,sim_name,'z%.1f'%zbox,name,get_hod_string(hod_params))
        config_fn = os.path.join(self.basedir,'configs.json'.format(name))
        self.config_fn = os.path.abspath(config_fn)
        if os.path.exists(config_fn):
            self.read_config(config_fn)
        else:
            print('Config file does not exist\nGenerate config file with __init__')

    @classmethod
    def initialize(cls, name, cosmo, d,d_spec,maskfn,cra, cdc, diam, zmin,zmax,hod_params,use_mask,
                 area=-1,ntarg_ratio=-1,nspec_ratio=-1,custom_ntarg=-1,custom_nspec=-1,nbox=128, fsamp=1,overwrite=False, basedir='.',
                 m_bands = ['M411','M438','M464','M490','M517'],zbox=3,sim_name='AbacusSummit_high_c000_ph100',gal_type='LRG'):
        # generate config file if it doesn't exist

        init_basedir = basedir + ''
        basedir = os.path.join(basedir,sim_name,'z%.1f'%zbox,name,get_hod_string(hod_params))
        config_fn = os.path.join(basedir,'configs.json'.format(name))
        config_fn = os.path.abspath(config_fn)
        if not os.path.exists(config_fn) or overwrite:
            cls.gen_config(config_fn, name, cosmo, d,d_spec,maskfn,cra, cdc, diam, zmin,zmax,hod_params,use_mask,
                            area=area,ntarg_ratio=ntarg_ratio,nspec_ratio=nspec_ratio,custom_ntarg=custom_ntarg,custom_nspec=custom_nspec,nbox=nbox,fsamp=fsamp,
                            m_bands=m_bands,zbox=zbox,sim_name=sim_name,gal_type=gal_type)
        # self.read_config(config_fn)

        forecast = cls(name,hod_params,basedir=init_basedir,zbox=zbox,sim_name=sim_name)
        return forecast
        
    @staticmethod
    def gen_config(config_fn, name, cosmo, d, d_spec, maskfn, cra, cdc, diam, zmin, zmax, hod_params, use_mask, 
                   area=-1,ntarg_ratio=-1,nspec_ratio=-1,custom_ntarg=-1,custom_nspec=-1,nbox=128,fsamp=1,
                   m_bands = ['M411','M438','M464','M490','M517'],zbox=3,sim_name='AbacusSummit_high_c000_ph100',gal_type='LRG'):
        # copy of generate_config.ipynb
        # with extra functionalities (use_mask, area, ntarg_ratio)
        
        mask = SurveyMask(maskfn)
        h = cosmo.h()
        def chi(zz):
            if zz>0:return (cosmo.comoving_distance(zz)*h)
            else: return (0)
        chimin, chimax = chi(zmin), chi(zmax)
        chi0 = (chimin+chimax)/2
        chirange = (chimax-chimin)/2
        # ndensity = d['RA'].size/mask.area()
            
        opt_list = []
        pes_list = []
        for band in m_bands:        
            # define min and max of redshift ranges
            if band in band_z.keys():
                zmin, zmax = band_z[band][0], band_z[band][1]
                dtmp = d[d['SEL_%s'%band]]
            else:
                zmin, zmax = band_z['M411'][0], band_z['M517'][1]
            # chimin_tmp, chimax_tmp = chi(zmin),chi(zmax)
    
            # dtmp['EFFTIME'] = 12.15 * dtmp['TSNR2_LRG']
            
            print(d['RA'].size)
            # restrict to objects with redshift and 2hrs observation
            dtmp = dtmp[dtmp['EFFTIME']>7200]
            n_tot = dtmp['RA'].size
            dtmp = dtmp[dtmp['RRJM_DELTACHI2']>30]
            n_rr = dtmp['RA'].size
            # restrict to objects with desired redshift range (2.26,3.41)
            dtmp = dtmp[(dtmp['RRJM_Z']>zmin)&(dtmp['RRJM_Z']<zmax)]
            # dat = dat[(dat['RRJM_Z']>2.)&(dat['RRJM_Z']<zmax)]
            n_good = dtmp['RA'].size
        
            print(band)
            print(n_tot,'&',n_rr,'&',n_good)
            print('pessimistic %.3f'%(1-n_good/n_tot))
            print('optimistic %.3f'%(1-n_good/n_rr))
        
            opt_list.append(1-n_good/n_rr)
            pes_list.append(1-n_good/n_tot)
    
        fint_list = np.mean([opt_list,pes_list],axis=0)

        ndensity_band = []
        dNdchi_band = []

        bin_width = 10
        chivec = np.linspace(chimin,chimax,100000)
        for i, band in enumerate(m_bands):
            print(band)
            dtmp_spec = d_spec[d_spec['SEL_%s'%band]]
            dtmp = d[d['SEL_%s'%band]]
            ndensity = dtmp['RA'].size/mask.area()
            
            print('area is %.2f, ndensity is %.1f'%(mask.area(),ndensity))
            ndensity_band.append(ndensity)
            
            ### dNdchi
            zmin_tmp, zmax_tmp = band_z[band][0], band_z[band][1]
            chimin_tmp, chimax_tmp = chi(zmin_tmp), chi(zmax_tmp)
            print(chimin_tmp,chimax_tmp)
    
            bins = np.arange(chimin,chimax, bin_width) # spaced out by chi=10
            hchi = np.histogram(dtmp_spec['CHI'],bins = bins) # normalize
            x = (hchi[1][1:]+hchi[1][:-1])/2
            y = hchi[0]
            spl = make_smoothing_spline(x, y, lam=1e3)
            probabilities = spl(chivec)
            probabilities[probabilities<0] = 0

            # recent change
            probabilities[chivec<chimin_tmp] *= 0
            probabilities[chivec>chimax_tmp] *= 0
            
            probabilities *= d['RA'].size/d_spec['RA'].size/mask.area() * (1-fint_list[i])
            dNdchi_band.append(probabilities)
            plt.plot(chivec,probabilities,c='C%d'%i)
            plt.hist(dtmp_spec['CHI'], bins = bins, 
                     weights = np.ones(dtmp_spec['RA'].size)*d['RA'].size/d_spec['RA'].size/mask.area()* (1-fint_list[i]),color='C%d'%i,alpha=.5,)
        plt.plot([],[],c='k',label='fit')
        plt.xlabel(r'$\chi$')
        plt.legend()
        plt.savefig('plots/downsample_dNdz_band.png')
        plt.show()
        plt.clf()

        dNdchi_band = np.array(dNdchi_band)
        
        config = {
            "ndensity_band": ndensity_band,
            "chi_vector":chivec.tolist(),
            'dNdchi_band':dNdchi_band.tolist(),
            "bin_width": bin_width,
            "maskfn": os.path.abspath(maskfn),
            "cra": cra,
            "cdc": cdc,
            'diam': diam,
            "zmin": zmin,
            "zmax": zmax,
            "chi0": chi0,
            "chirange": chirange,
            "chimin": chimin,
            "chimax": chimax,
            "hod_params": hod_params,
            "cosmo_params":cosmo.pars,
            "interloper_band":fint_list.tolist(),
            'use_mask':use_mask,
            'ntarg_ratio':ntarg_ratio,
            'nspec_ratio':nspec_ratio,
            'custom_ntarg':custom_ntarg,
            'custom_nspec':custom_nspec,
            'fsamp':fsamp,
            'name':name,
            'nbox':nbox,
            'bands':m_bands,
            'zbox':zbox,
            'sim_name':sim_name,
            'gal_type':gal_type,
        }
        
        # Write the dictionary to a JSON file.
        try:
            os.makedirs(os.path.dirname(config_fn),exist_ok=True)
            with open(config_fn, "w") as f:
                json.dump(config, f, indent=4)
            print(f"Configuration successfully written to {config_fn}")
        except Exception as e:
            print(f"An error occurred while writing to {config_fn}: {e}")
        
    def read_config(self,config_fn,):
        # read config file
        self.config = json.load(open(config_fn,'r'))

        '''
        self.ndensity_band = config['ndensity_band']
        self.chi_vector = config['chi_vector']
        self.dNdchi_band = config['dNdchi_band']
        self.bin_width = config['bin_width']
        self.maskfn = config['maskfn']
        self.zmin = config['zmin']
        self.zmax = config['zmax']
        self.chi0 = config['chi0']
        self.chirange = config['chirange']
        self.chimin = config['chimin']
        self.chimax = config['chimax']
        self.hod_params = config['hod_params']
        self.cosmo_params = config['cosmo_params']
        self.interloper_band = config['interloper_band']
        self.use_mask = config['use_mask']
        self.ntarg_ratio = config['ntarg_ratio']
        self.nspec_ratio = config['nspec_ratio']
        self.custom_ntarg = config['custom_ntarg']
        self.custom_nspec = config['custom_nspec']
        '''
        
    #########
    # UTILS #
    #########

    def check_hod_exists(self,):
        name = self.config['name']
        nbox = self.config['nbox']
        hod_string = get_hod_string(self.config['hod_params'])
        bb = True
        for i in range(nbox):
            bb &= os.path.exists(os.path.join(self.basedir,'inputcats/dat%d.fits'%i))
        return bb

    def check_cat_exists(self,):
        name = self.config['name']
        nbox = self.config['nbox']
        hod_string = get_hod_string(self.config['hod_params'])
        bb = True
        for i in range(nbox):
            bb &= os.path.exists(os.path.join(self.basedir,'inputcats/tot%d.fits'%i)) # includes interlopers
        return bb

        
    def get_hod_spec(self,n):
        name = self.config['name']
        hod_string = get_hod_string(self.config['hod_params'])
        datfn = os.path.join(self.basedir,'inputcats/dat%d.fits'%n)
        dat = Table.read(datfn)
        return dat
        
    def get_hod_targ(self,n):
        name = self.config['name']
        hod_string = get_hod_string(self.config['hod_params'])
        datfn = os.path.join(self.basedir,'inputcats/tot%d.fits'%n) # includes interlopers
        dat = Table.read(datfn)
        return dat

    def get_hod_ran(self):
        name = self.config['name']
        hod_string = get_hod_string(self.config['hod_params'])
        datfn = os.path.join(self.basedir,'inputcats/ran.fits') 
        dat = Table.read(datfn)
        return dat

    
    ############################
    # GENERATING MOCK CATALOGS #
    ############################
    def make_hod(self,sbatch=True,overwrite=False,jobtype='debug'):
        if self.check_hod_exists() and not overwrite: return 

        if not sbatch:
            homedir = '/pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/mc'
            content=f"""
module load python
conda activate myenv_abacus

cd {homedir}
export OMP_NUM_THREADS=16
# argv[1] is path to config file
echo {self.config_fn}
python save_hod.py {self.config_fn}
            """
            sbatch_filename = "forecast_auto_2nd.sh"
            with open(sbatch_filename, "w") as f:
                f.write(content)

            
            bashCommand = "chmod +x ./{:s}".format(sbatch_filename)
            os.system(bashCommand)
            bashCommand = "./{:s}".format(sbatch_filename)
            os.system(bashCommand)
            return

            
        job_name = "forecast_sbatch"
        output_file = "slurm-%j.out"
        time = "00:15:00"  # HH:MM:SS
        script_to_run = "my_script.py"  # Replace with your actual job script
        homedir = '/pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/mc'

        sbatch_content = f"""#!/bin/bash -l
#SBATCH -J {job_name}
#SBATCH -N 1
#SBATCH -t {time}
#SBATCH -o clustering.out
#SBATCH -e clustering.err
#SBATCH -q {jobtype}
#SBATCH -C cpu
#SBATCH -A m68

module load python
conda activate myenv_abacus

cd {homedir}

# argv[1] is path to config file
date
echo {self.config_fn}
python save_hod.py {self.config_fn}
date
"""
        
        sbatch_filename = "forecast_auto_sbatch.sh"
        with open(sbatch_filename, "w") as f:
            f.write(sbatch_content)
            
        bashCommand = "sbatch {:s}".format(sbatch_filename)
        os.system(bashCommand)

    def make_cat(self,overwrite=False):
        '''
        Make the targ catalog, as a combination of hod + interlopers
        '''
        assert self.check_hod_exists()

        if self.check_cat_exists() and not overwrite: return
        # os.system('module load python')
        # os.system('conda activate myenv')
        # os.system('python save_interloper.py {:s}'.format(self.config_fn))

        self.save_interloper()
        return 
    

    def save_interloper(self,):
        config = self.config
        hod_params = config['hod_params']
        ndensity_band = config['ndensity_band'] # target density
        ndensity_band = np.array(ndensity_band) 
        name = config['name']
        nbox = config['nbox']
        maskfn = config['maskfn']
        use_mask = config['use_mask']
    
        ntarg_ratio = config['ntarg_ratio'] # use this for now
        # nspec_ratio = config['nspec_ratio']
        custom_ntarg = config['custom_ntarg'] # use this for now
        # custom_nspec = config['custom_nspec']
    
        cra = config['cra']
        cdc = config['cdc']
        diam = config['diam']

        m_bands = config['bands']
    
        if ntarg_ratio == -1 and custom_ntarg == -1: ntarg_ratio = 1.
        elif ntarg_ratio != -1: assert custom_ntarg == -1
        else: ntarg_ratio = custom_ntarg / np.sum(ndensity_band)
    
        hod_string = get_hod_string(self.config['hod_params'])
        
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
            
        outdir = os.path.join(self.basedir,'inputcats')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        
        intf_band = config['interloper_band'] # interloper fraction

        nints_band = [[] for i in range(len(m_bands))]
        for i in range(nbox):
            ran_tables = []
            for j, band in enumerate(m_bands):
                rng    = np.random.default_rng(i*len(m_bands)+j) # avoid all the bands having the same interloper positions
                mask = SurveyMask(maskfn)
                
                nran = 100000
                ran  = {}

                diam_deg = diam * 180/np.pi
                
                ran['RA' ] = rng.uniform(low=-diam_deg/2,high=diam_deg/2 ,size=nran)
                ran['DEC' ] = rng.uniform(low=-diam_deg/2,high=diam_deg/2 ,size=nran)
                ww, nra, ndc = in_field(ran)
                ran['RA'] = nra
                ran['DEC'] = ndc
                
                ran['CHI'] = -1 * np.ones(ran['RA'].size)
                
                fsamp = ndensity_band[j] / (ran['RA'].size / area ) * intf_band[j] * ntarg_ratio
                rand = rng.uniform(low=0,high=1,size=ran['RA'].size)
                ww   = np.nonzero( rand<fsamp )[0]
        
                for key in ran.keys():
                    ran[key] = ran[key][ww]
                    
                for j2, band2 in enumerate(m_bands):
                    if j==j2: ran['SEL_%s'%band2] = np.ones(ran['RA'].size,dtype=bool)
                    else: ran['SEL_%s'%band2] = np.zeros(ran['RA'].size,dtype=bool)
                ran['SEL'] = np.ones(ran['RA'].size,dtype=bool)
        
                ran_table = Table()
                for key in ran.keys():
                    ran_table[key] = ran[key]
                # print(band,ran['RA'].size)
                nints_band[j].append(ran['RA'].size)
                ran_tables.append(ran_table)
            
            ran_table = vstack(ran_tables)
            ran_table.write(os.path.join(outdir,'int%i.fits'%i),overwrite=True)
                
            dat_table = Table.read(os.path.join(outdir,'dat%d.fits'%i))
            for key in dat_table.keys():
                if key not in ran_table.keys():
                    ran_table[key] = np.zeros(ran_table['RA'].shape)
    
            dat_table = vstack([dat_table,ran_table])
            dat_table.write(os.path.join(outdir,'tot%d.fits'%i), format='fits',overwrite=True)
            
            del dat_table, ran_table

        for j, band in enumerate(m_bands):
            print('for {:s}, interloper size is {:.2f}+-{:.2f}'.format(band,np.mean(nints_band,axis=1)[j],np.std(nints_band,axis=1)[j]))
            print('with angular density {:.2f}+-{:.2f}'.format(np.mean(nints_band,axis=1)[j]/area,np.std(nints_band,axis=1)[j]/area))


    ########################
    # COMPUTING CLUSTERING #
    ########################
    def compute_wp(self,overwrite=False):
        '''
        Compute the projected correlation function w_p(R) for each box and compute covariance
        '''
        config = self.config

        
        # downsample from targ to spec...?

        # both ran and dat have masks applied
        # ranfn = os.path.join(catdir,'ran.fits')
        # ran = Table.read(ranfn)
        ran = self.get_hod_ran()
        
        tval,wps ,ngals = None,[],[]
        nbox = config['nbox']
        maskfn = config['maskfn']
        use_mask = config['use_mask']
        diam = config['diam']
        name = config['name']
        fsamp = config['fsamp']
        hod_string = get_hod_string(self.config['hod_params'])
        m_bands = config['bands']
        # outdir = 'hod/{:s}/{:s}/clustering'.format(name,hod_string)
        outdir = os.path.join(self.basedir,'clustering')
        if os.path.exists(os.path.join(outdir,'nbox{:d}_wp.txt'.format(nbox))) and not overwrite: 
            print('wp already computed')
            return

        
        for i in range(nbox):
            # spec
            dat = self.get_hod_spec(i)
            # targ
            # targ = self.get_hod_targ(i)
            targ = self.get_hod_spec(i)
            
            ### IF I WANT TARG -> SPEC DOWNSAMPLING ###
            # rng    = np.random.default_rng(1)
            # rand = rng.uniform(low=0,high=1,size=dat['RA'].size)
            # ww   = np.nonzero( rand<fsamp )[0]
            # dat = dat[ww]
        
            ngals.append(dat['RA'].size)
            
            # r_bins = np.logspace(0,1.5,8)
            r_bins = np.geomspace(0.4,45,9)[:-1]

        
            band_chi = {
                'M411':[3821.6517, 4069.5162],
                'M438':[3998.3714, 4225.9277],
                'M464':[4160.4664, 4370.2586],
                'M490':[4309.7838, 4510.0708],
                'M517':[4447.8731, 4633.9157],
            }
            
            wwt_cross = []
            for j,band in enumerate(m_bands):
                chimin, chimax = band_chi[band]
                chi0 = (chimin+chimax)/2
                chirange = (chimax-chimin)/2
                
                theta_bins = r_bins/chi0* 180/np.pi
                
                twt, tmp = calc_wt_cross(dat[dat['SEL_%s'%band]],targ,ran,bins=theta_bins)
                # twt, tmp = calc_wt(dat[dat['SEL_%s'%band]],targ,ran,bins=theta_bins)
                wwt_cross.append(tmp)
            f_chi = [targ[targ['SEL_%s'%band]]['RA'].size for j, band in enumerate(m_bands)]
            f_chi = np.array(f_chi,dtype=np.float_)
            N_chi = f_chi + 0.
            # int_f_chi = 0
            # for j,band in enumerate(m_bands):
            #     int_f_chi += (band_chi[band][1]-band_chi[band][0])*f_chi[j]
            # f_chi /= int_f_chi

            ### new
            # f_chi = []
            # for j,band in enumerate(m_bands):
            #     int_f_chi = np.trapz(self.config['dNdchi_band'][j],self.config['chi_vector'])
            #     f_chi.append(interp1d(self.config['chi_vector'],self.config['dNdchi_band'][j])(np.average(band_chi[band]))/int_f_chi)
            # f_chi=np.array(f_chi)

            # new new
            y = 0
            for j,band in enumerate(m_bands): y+= np.array(self.config['dNdchi_band'][j])
            int_f_chi = np.trapz(y,self.config['chi_vector'])
            f_chi = []
            for j,band in enumerate(m_bands):
                # f_chi.append(np.trapz(self.config['dNdchi_band'][j],self.config['chi_vector'])/int_f_chi)
                f_chi.append(np.trapz(self.config['dNdchi_band'][j],self.config['chi_vector'])/int_f_chi/(band_chi[band][1]-band_chi[band][0]))
                # intrange = (np.array(self.config['chi_vector'])>band_chi[band][0])&(np.array(self.config['chi_vector'])<band_chi[band][1])
                # f_chi.append(np.trapz(y[intrange],np.array(self.config['chi_vector'])[intrange])/int_f_chi/(band_chi[band][1]-band_chi[band][0]))
                # f_chi.append(interp1d(self.config['chi_vector'],self.config['dNdchi_band'][j])(np.average(band_chi[band]))/int_f_chi)
            f_chi=np.array(f_chi)

            
            ci = N_chi * f_chi / np.sum(N_chi*f_chi**2)

            ci = []
            ci_dict = json.load(open('/pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/data_measurements/ci.json','r'))
            for j,band in enumerate(m_bands): ci.append(ci_dict[band])
            ci = np.array(ci)

            
            # ci = 1/f_chi
            # ci = 1/2 * np.ones(2)
            print(f_chi,N_chi,ci)
            rwx = r_bins
            wp_manual = np.einsum('i,ij->j',ci,wwt_cross)
            
            wps.append(wp_manual)
            del dat, targ
        del ran
            
        rval = np.sqrt(rwx[:-1]*rwx[1:])
        wps  = np.array(wps)
        wpavg = np.mean(wps,axis=0)
        wperr = np.std( wps,axis=0)
        wpcor = np.corrcoef(wps,rowvar=False)
        navg = np.mean(np.array(ngals,dtype='float'))
        nerr = np.std( np.array(ngals,dtype='float'))
        if use_mask: 
            assert os.path.exists(maskfn)
            mask = SurveyMask(maskfn)
            area = mask.area()
        else: 
            area = (diam * 180/np.pi /2)**2 * np.pi

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        if use_mask: 
            assert os.path.exists(maskfn)
            mask = SurveyMask(maskfn)
            area = mask.area()
        else: 
            area = (diam * 180/np.pi /2)**2 * np.pi

        np.savetxt(os.path.join(outdir,'nbox{:d}_wp.txt'.format(nbox)),wps)
        np.savetxt(os.path.join(outdir,'nbox{:d}_r.txt'.format(nbox)),rval)
        
        
        with open(os.path.join(outdir,'nbox{:d}_wp_summary.txt'.format(nbox)),"w") as fout:
            fout.write("# Monte-Carlo calculation of wt using {:d} mocks.\n".\
                       format(wps.shape[0]))
            fout.write("# Field is area {:.2f}deg2.\n".\
                       format(area))
            fout.write("# Have {:.1f}+/-{:.2f} LBGs/field.\n".\
                       format(navg,nerr))
            fout.write("# Correlation matrix is:\n")
            for i in range(round(rval.size)):
                outstr = "#"
                for j in range(round(rval.size)): outstr += " {:8.4f}".format(wpcor[i,j])
                fout.write(outstr + "\n")
            fout.write("# {:>8s} {:>15s} {:>15s}\n".\
                       format("R[Mpc/h]","wt","dwt"))
            for i in range(rval.size):
                outstr = "{:10.3f} {:15.5e} {:15.5e}".format(rval[i],wpavg[i],wperr[i])
                fout.write(outstr+"\n")

        cov = np.zeros_like(wpcor)
        for i in range(rval.size):
            for j in range(rval.size):
                cov[i,j]=wperr[i]*wpcor[i,j]*wperr[j]
        np.savetxt(os.path.join(outdir,'nbox{:d}_wp_cov.txt'.format(nbox)),cov)

    def compute_wt(self,overwrite=False):
        '''
        Compute the projected correlation function w(theta) (targ x targ) for each box and compute covariance
        '''
        config = self.config
        
        # both ran and dat have masks applied
        ran = self.get_hod_ran()
        
        tval,wts ,ngals = None,[],[]
        nbox = config['nbox']
        maskfn = config['maskfn']
        use_mask = config['use_mask']
        diam = config['diam']
        name = config['name']
        hod_string = get_hod_string(self.config['hod_params'])
        m_bands = config['bands']
        outdir = os.path.join(self.basedir,'clustering')
        if os.path.exists(os.path.join(outdir,'wt_{:s}.txt'.format(m_bands[0]))) and not overwrite:
            print('wt already computed')
            return
        
        bins = np.geomspace(5e-3,6e-1,9)
        
        for i in range(nbox):
            # targ
            # targ = self.get_hod_targ(i)
            targ = self.get_hod_spec(i)
            
            for j,band in enumerate(m_bands):      
                twt,tmp = calc_wt(targ[targ['SEL_%s'%band]],ran[::50],bins=bins)
                if i==0: 
                    wts.append([])
                    ngals.append([])
                wts[j].append(tmp)
                ngals[j].append(targ[targ['SEL_%s'%band]]['RA'].size)
            del targ
        del ran
        wts  = np.array(wts) # (num band, nbox, w(theta))
        wtavg = np.mean(wts,axis=1)  # average over boxes
        wterr = np.std( wts,axis=1)
        wtcor = np.array([np.corrcoef(wts[j],rowvar=False) for j, band in enumerate(m_bands)])
        cov = np.zeros_like(wtcor)
        for j,band in enumerate(m_bands):
            for i in range(twt.size):
                for j2 in range(twt.size):
                    cov[j,i,j2]=wterr[j,i]*wtcor[j,i,j2]*wterr[j,j2]
                
        navg = np.mean(np.array(ngals,dtype='float'),axis=1)
        nerr = np.std( np.array(ngals,dtype='float'),axis=1)
        if use_mask: 
            assert os.path.exists(maskfn)
            mask = SurveyMask(maskfn)
            area = mask.area()
        else: 
            area = (diam * 180/np.pi /2)**2 * np.pi

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        if use_mask: 
            assert os.path.exists(maskfn)
            mask = SurveyMask(maskfn)
            area = mask.area()
        else: 
            area = (diam * 180/np.pi /2)**2 * np.pi

        np.savetxt(os.path.join(outdir,'theta.txt'),twt)
        for j, band in enumerate(m_bands):
            np.savetxt(os.path.join(outdir,'wt_{:s}.txt'.format(band)),wts[j])
        
        
            with open(os.path.join(outdir,'wt_{:s}_summary.txt'.format(band)),"w") as fout:
                fout.write("# Monte-Carlo calculation of wt using {:d} mocks.\n".\
                           format(wts[j].shape[0]))
                fout.write("# Field is area {:.2f}deg2.\n".\
                           format(area))
                fout.write("# Have {:.1f}+/-{:.2f} LBGs/field.\n".\
                           format(navg[j],nerr[j]))
                fout.write("# Correlation matrix is:\n")
                for i in range(round(twt.size)):
                    outstr = "#"
                    for j2 in range(round(twt.size)): outstr += " {:8.4f}".format(wtcor[j,i,j2])
                    fout.write(outstr + "\n")
                fout.write("# {:>8s} {:>15s} {:>15s}\n".\
                           format("R[Mpc/h]","wt","dwt"))
                for i in range(twt.size):
                    outstr = "{:10.3f} {:15.5e} {:15.5e}".format(twt[i],wtavg[j,i],wterr[j,i])
                    fout.write(outstr+"\n")

            np.savetxt(os.path.join(outdir,'wt_{:s}_cov.txt'.format(band)),cov[j])

    def sum_inverse_variance(self,datvec,covvec):
        tot_cov = 0
        for cov in covvec: tot_cov += np.linalg.inv(cov)
        tot_cov = np.linalg.inv(tot_cov) # C^-1 = C1^-1 + C2^-1

        tot_dat = 0
        for i, (dat,cov) in enumerate(zip(datvec,covvec)):
            tot_dat += np.einsum('ij,j->i',np.linalg.inv(cov),dat)
        tot_dat = np.einsum('ij,j->i',tot_cov,tot_dat) # x = C[C1^-1 x1 + C2^-1 x2]

        return tot_dat, tot_cov
    
    def compute_wR(self,overwrite=False):
        '''
        Basically the same as compute_wt, but with fixed rbins for all z
        '''
        config = self.config
        
        # both ran and dat have masks applied
        ran = self.get_hod_ran()
        
        tval,wts ,ngals = None,[],[]
        nbox = config['nbox']
        maskfn = config['maskfn']
        use_mask = config['use_mask']
        diam = config['diam']
        name = config['name']
        hod_string = get_hod_string(self.config['hod_params'])
        m_bands = config['bands']
        outdir = os.path.join(self.basedir,'clustering')
        if os.path.exists(os.path.join(outdir,'wR_{:s}.txt'.format(m_bands[0]))) and not overwrite:
            print('wR already computed')
            return
        rbins = np.geomspace(0.4,45,9)[:-1]
        # bins = np.geomspace(5e-3,6e-1,9)
        
        for i in range(nbox):
            # targ
            # targ = self.get_hod_targ(i)
            targ = self.get_hod_spec(i)
            ran = Table.read('/pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/mc/AbacusSummit_high_c000_ph100/z3.0/ibis_tertiary44_clustering/logM_cut_11.75_logM1_12.75_sigma_0.66_kappa_1.00_alpha_0.66/inputcats/ran%d.fits'%i)
            
            for j,band in enumerate(m_bands):      
                bins = rbins / (band_chimid[band]) * 180/np.pi
                twt,tmp = calc_wt(targ[targ['SEL_%s'%band]],ran[j:][::5],bins=bins)
                if i==0: 
                    wts.append([])
                    ngals.append([])
                wts[j].append(tmp)
                ngals[j].append(targ[targ['SEL_%s'%band]]['RA'].size)
            del targ
        del ran
        wts  = np.array(wts) # (num band, nbox, w(theta))
        wtavg = np.mean(wts,axis=1)  # average over boxes
        wterr = np.std( wts,axis=1)
        wtcor = np.array([np.corrcoef(wts[j],rowvar=False) for j, band in enumerate(m_bands)])
        cov = np.zeros_like(wtcor)
        for j,band in enumerate(m_bands):
            for i in range(twt.size):
                for j2 in range(twt.size):
                    cov[j,i,j2]=wterr[j,i]*wtcor[j,i,j2]*wterr[j,j2]
                
        navg = np.mean(np.array(ngals,dtype='float'),axis=1)
        nerr = np.std( np.array(ngals,dtype='float'),axis=1)
        if use_mask: 
            assert os.path.exists(maskfn)
            mask = SurveyMask(maskfn)
            area = mask.area()
        else: 
            area = (diam * 180/np.pi /2)**2 * np.pi

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        if use_mask: 
            assert os.path.exists(maskfn)
            mask = SurveyMask(maskfn)
            area = mask.area()
        else: 
            area = (diam * 180/np.pi /2)**2 * np.pi

        RwR = np.sqrt(rbins[1:]*rbins[:-1])
        np.savetxt(os.path.join(outdir,'R.txt'),RwR)
        for j, band in enumerate(m_bands):
            np.savetxt(os.path.join(outdir,'wR_{:s}.txt'.format(band)),wts[j])
        
            with open(os.path.join(outdir,'wR_{:s}_summary.txt'.format(band)),"w") as fout:
                fout.write("# Monte-Carlo calculation of wt using {:d} mocks.\n".format(wts[j].shape[0]))
                fout.write("# Field is area {:.2f}deg2.\n".format(area))
                fout.write("# Have {:.1f}+/-{:.2f} LBGs/field.\n".format(navg[j],nerr[j]))
                fout.write("# Correlation matrix is:\n")
                for i in range(round(RwR.size)):
                    outstr = "#"
                    for j2 in range(round(RwR.size)): outstr += " {:8.4f}".format(wtcor[j,i,j2])
                    fout.write(outstr + "\n")
                fout.write("# {:>8s} {:>15s} {:>15s}\n".format("R[Mpc/h]","wR","dwt"))
                for i in range(twt.size):
                    outstr = "{:10.3f} {:15.5e} {:15.5e}".format(RwR[i],wtavg[j,i],wterr[j,i])
                    fout.write(outstr+"\n")

            np.savetxt(os.path.join(outdir,'wR_{:s}_cov.txt'.format(band)),cov[j])
        totwR, totcov = self.sum_inverse_variance(wtavg,cov)
        np.savetxt(os.path.join(outdir,'wR_tot.txt'),totwR)
        np.savetxt(os.path.join(outdir,'wR_tot_cov.txt'),totcov)
            
    def compute_xi(self,zbins=None,overwrite=False,s_bins=None):
        '''
        Compute the 3d correlation function xi_ell(s) for each box and compute covariance
        '''
        config = self.config
        m_bands = config['bands']
        if zbins == None: zbins = range(len(m_bands))
        zmin = list(band_z.values())[zbins[0]][0]
        zmax = list(band_z.values())[zbins[-1]][1]
        zmid = (zmin+zmax)/2
        print(m_bands,zmin,zmid,zmax)
        # downsample from targ to spec...?

        # both ran and dat have masks applied
        # ranfn = os.path.join(catdir,'ran.fits')
        # ran = Table.read(ranfn)
        ran = self.get_hod_ran()
        rans = []
        for j, band in enumerate(m_bands): 
            if j in zbins: rans.append(ran[ran['SEL_%s'%band]])
        ran = vstack(rans)
        del rans
        if np.min(ran['RA'])<0:
            ran['RA']+=180
        
        rval,xi0s,xi2s,xi4s,ngals = None,[],[],[],[]
        nbox = config['nbox']
        maskfn = config['maskfn']
        use_mask = config['use_mask']
        diam = config['diam']
        name = config['name']
        fsamp = config['fsamp']
        hod_string = get_hod_string(self.config['hod_params'])
        outdir = os.path.join(self.basedir, 'clustering')

        if os.path.exists(os.path.join(outdir,'z{:.3f}_xi.txt'.format(zmid))) and not overwrite:
            print('xi already computed')
            return
        if s_bins is None: r_bins = np.geomspace(1,50,11)
        else: r_bins = s_bins

            
        for i in range(nbox):
            # spec
            dat = self.get_hod_spec(i)
            print(dat['RA'].size)
            # # targ
            # targ = self.get_hod_targ(i)
            dats = []
            for j, band in enumerate(m_bands): 
                if j in zbins: dats.append(dat[dat['SEL_%s'%band]])
            dat = vstack(dats)
            del dats
            print(dat['RA'].size)
            
            
            ### IF I WANT TARG -> SPEC DOWNSAMPLING ###
            rng    = np.random.default_rng(1)
            rand = rng.uniform(low=0,high=1,size=dat['RA'].size)
            ww   = np.nonzero( rand<fsamp )[0]
            dat = dat[ww]
        
            ngals.append(dat['RA'].size)
            # r_bins = np.logspace(0,1.5,8)
            band_chi = {
                'M411':[3821.6517, 4069.5162],
                'M438':[3998.3714, 4225.9277],
                'M464':[4160.4664, 4370.2586],
                'M490':[4309.7838, 4510.0708],
                'M517':[4447.8731, 4633.9157],
            }

            # compute xi with all bands together... for now
            if np.min(dat['RA'])<0:
                dat['RA']+=180
            
            rxi, xi0, xi2, xi4 = calc_xi(dat,ran[::10],bins=r_bins)
                        
            xi0s.append(xi0)
            xi2s.append(xi2)
            xi4s.append(xi4)
            del dat
        del ran
            
        rval = np.sqrt(rxi[:-1]*rxi[1:])
        xi0s  = np.array(xi0s)
        xi2s  = np.array(xi2s)
        xi4s  = np.array(xi4s)
        xis = np.concatenate((xi0s,xi2s,xi4s),axis=1) # (nbox, rval) *3 --> (nbox,3*rval) 
        xi0avg = np.mean(xi0s,axis=0)
        xi2avg = np.mean(xi2s,axis=0)
        xi4avg = np.mean(xi4s,axis=0)
        xiavg = np.mean(xis,axis=0)
        xi0err = np.std( xi0s,axis=0)
        xi2err = np.std( xi2s,axis=0)
        xi4err = np.std( xi4s,axis=0)
        xierr = np.std( xis,axis=0)
        xicor = np.corrcoef(xis,rowvar=False)
        navg = np.mean(np.array(ngals,dtype='float'))
        nerr = np.std( np.array(ngals,dtype='float'))
        if use_mask: 
            assert os.path.exists(maskfn)
            mask = SurveyMask(maskfn)
            area = mask.area()
        else: 
            area = (diam * 180/np.pi /2)**2 * np.pi

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        if use_mask: 
            assert os.path.exists(maskfn)
            mask = SurveyMask(maskfn)
            area = mask.area()
        else: 
            area = (diam * 180/np.pi /2)**2 * np.pi

        np.savetxt(os.path.join(outdir,'z{:.3f}_xi.txt'.format(zmid)),xis)
        np.savetxt(os.path.join(outdir,'z{:.3f}_xi0.txt'.format(zmid)),xi0s)
        np.savetxt(os.path.join(outdir,'z{:.3f}_xi2.txt'.format(zmid)),xi2s)
        np.savetxt(os.path.join(outdir,'z{:.3f}_xi4.txt'.format(zmid)),xi4s)
        np.savetxt(os.path.join(outdir,'z{:.3f}_xi_r.txt'.format(zmid)),rval)
        np.savetxt(os.path.join(outdir,'z{:.3f}_xi_avg.txt'.format(zmid)),xiavg)

        
        with open(os.path.join(outdir,'z{:.3f}_xi_summary.txt'.format(zmid)),"w") as fout:
            fout.write("# Monte-Carlo calculation of wt using {:d} mocks.\n".\
                       format(xis.shape[0]))
            fout.write("# Field is area {:.2f}deg2.\n".\
                       format(area))
            fout.write("# Have {:.1f}+/-{:.2f} LBGs/field.\n".\
                       format(navg,nerr))
            fout.write("# Correlation matrix is:\n")
            for i in range(3*round(rval.size)):
                outstr = "#"
                for j in range(3*round(rval.size)): outstr += " {:8.4f}".format(xicor[i,j])
                fout.write(outstr + "\n")
            fout.write("# {:>8s} {:>15s} {:>15s}\n".\
                       format("R[Mpc/h]","wt","dwt"))
            for i in range(rval.size):
                outstr = "{:10.3f} {:15.5e} {:15.5e}".format(rval[i],xi0avg[i],xi0err[i])
                fout.write(outstr+"\n")
            for i in range(rval.size):
                outstr = "{:10.3f} {:15.5e} {:15.5e}".format(rval[i],xi2avg[i],xi2err[i])
                fout.write(outstr+"\n")
            for i in range(rval.size):
                outstr = "{:10.3f} {:15.5e} {:15.5e}".format(rval[i],xi4avg[i],xi4err[i])
                fout.write(outstr+"\n")

        cov = np.zeros_like(xicor)
        for i in range(3*rval.size):
            for j in range(3*rval.size):
                cov[i,j]=xierr[i]*xicor[i,j]*xierr[j]
        np.savetxt(os.path.join(outdir,'z{:.3f}_xi_cov.txt'.format(zmid)),cov)

    def compute_dNdz(self,):
        config = self.config
        nbox = config['nbox']
        chimin = self.config['chimin']
        chimax = self.config['chimax']
        bin_width = self.config['bin_width']
        m_bands = config['bands']
        cosmo = Class()
        cosmo.set(config['cosmo_params'])
        cosmo.compute()
        h = cosmo.h()
        chi = lambda z: cosmo.comoving_distance(z)*h if z>0 else -1
        bins = np.arange(chimin,chimax, bin_width) # spaced out by chi=10
        dNdchis = []
        for i in range(nbox):
            spec = self.get_hod_spec(i)
            for j, band in enumerate(m_bands):
                dtmp_spec = spec[spec['SEL_%s'%band]]
                zmin_tmp, zmax_tmp = band_z[band][0], band_z[band][1]
                chimin_tmp, chimax_tmp = chi(zmin_tmp), chi(zmax_tmp)
                
                if i==0: dNdchis.append([])
                hchi = np.histogram(dtmp_spec['CHI'],bins=bins)
                x = (hchi[1][1:]+hchi[1][:-1])/2
                y = hchi[0]
                y[x<chimin_tmp] *= 0
                y[x>chimax_tmp] *= 0
                dNdchis[j].append(y)
                
        dNdchis  = np.array(dNdchis) # (num band, nbox, w(theta))
        dNdchiavg = np.mean(dNdchis,axis=1)  # average over boxes
        dNdchierr = np.std( dNdchis,axis=1)
        outdir = os.path.join(self.basedir,'clustering')
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        np.savetxt(os.path.join(outdir,'chi.txt'),x)
        for j, band in enumerate(m_bands):
            np.savetxt(os.path.join(outdir,'dNdchi_avg_{:s}.txt'.format(band)),dNdchiavg[j])
            np.savetxt(os.path.join(outdir,'dNdchi_err_{:s}.txt'.format(band)),dNdchierr[j])

        return 

    #########
    # CHAIN #
    #########

    def rewrite_xi_yaml(self,yamlfn,zbins=None,linear=False,chain=False):
        m_bands = self.config['bands']
        if zbins == None: zbins = range(len(m_bands))

        print(list(band_z.values()))
        zmin = list(band_z.values())[zbins[0]][0]
        zmax = list(band_z.values())[zbins[-1]][1]

        zmid = (zmin+zmax)/2

        name = self.config['name']
        hod_string = get_hod_string(self.config['hod_params'])
        nbox = self.config['nbox']
        
        # outdir = 'hod/{:s}/{:s}/cobaya'.format(name,hod_string)
        # datdir = 'hod/{:s}/{:s}/clustering'.format(name,hod_string)
        outdir = os.path.join(self.basedir,'cobaya')
        datdir = os.path.join(self.basedir,'clustering')

        yamldir = os.path.dirname(yamlfn)
        # newyamlfn = os.path.join(yamldir,

        tmp = 'z{:.3f}'.format(zmid)
        if linear: tmp += '_linear'
        if chain: tmp += '_chain'

        newyamlfn = os.path.join(outdir,'yamls',tmp+'.yaml')
        
        try:
            from cobaya.yaml import yaml_load_file
            yy = yaml_load_file(yamlfn)
        except:
            import yaml
            with open(yamlfn) as stream:
                yy = (yaml.safe_load(stream))

        yy['theory']['likelihoods.xi_likelihood.PT_xi_theory']['zfid']=zmid
        yy['likelihood']['likelihoods.xi_likelihood.FSLikelihood']['zfid']=zmid
        
        yy['likelihood']['likelihoods.xi_likelihood.FSLikelihood']['xifn']=os.path.abspath(os.path.join(datdir,'z{:.3f}_xi_avg.txt'.format(zmid)))
        yy['likelihood']['likelihoods.xi_likelihood.FSLikelihood']['rfn']=os.path.abspath(os.path.join(datdir,'z{:.3f}_xi_r.txt'.format(zmid)))
        yy['likelihood']['likelihoods.xi_likelihood.FSLikelihood']['covfn']=os.path.abspath(os.path.join(datdir,'z{:.3f}_xi_cov.txt'.format(zmid)))

        if linear: 
            nonlinear_params = ['b2','bs','alpha0','alpha2','SN2']
            for pp in nonlinear_params:
                del yy['params'][pp]['prior']
                del yy['params'][pp]['ref']
                yy['params'][pp]['value'] = 0.0
        
        if linear: yy['output'] = os.path.abspath(os.path.join(outdir,'output','z{:.3f}_linear'.format(zmid)))
        else: yy['output'] = os.path.abspath(os.path.join(outdir,'output','z{:.3f}'.format(zmid)))

        if not os.path.exists(os.path.dirname(newyamlfn)): os.makedirs(os.path.dirname(newyamlfn))
        with open(newyamlfn, 'w') as outfile:
            yaml.dump(yy, outfile, default_flow_style=False, sort_keys=False)

        return newyamlfn

    def run_xi_minimize(self,zbins=None,jobtype='debug',linear=False,yamlfn='/pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/mc/cobaya/yamls/fs_minimize.yaml'):
        name = self.config['name']
        hod_string = get_hod_string(self.config['hod_params'])
        nbox = self.config['nbox']
        # outdir = 'hod/{:s}/{:s}/cobaya'.format(name,hod_string)
        outdir = os.path.join(self.basedir,'cobaya')

        newyamlfn = self.rewrite_xi_yaml(yamlfn,zbins=zbins,linear=linear)
        newyamlfn = os.path.abspath(newyamlfn)

        job_name = "xi_{:s}".format(name)
        output_file = 'xi_ell' # os.path.join(outdir,'xi_ell')
        time = "00:30:00"  
        workdir = '/pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/mc/cobaya'
        sbatch_content = f"""#!/bin/bash -l
#SBATCH -J {job_name}
#SBATCH -N 1
#SBATCH -t {time}
#SBATCH -o {output_file}.out
#SBATCH -e {output_file}.err
#SBATCH -q {jobtype}
#SBATCH -C cpu
#SBATCH -A m68

module load python
source activate cobaya_07292024

export PYTHONPATH={workdir}:$PYTHONPATH

cd {workdir}
export OMP_NUM_THREADS=8
# srun -n 16 -c $OMP_NUM_THREADS cobaya-run -f {newyamlfn}
srun -n 16 -c $OMP_NUM_THREADS cobaya-run -f {newyamlfn}
"""
        
        sbatch_filename = "xi_chain.sh"
        with open(sbatch_filename, "w") as f:
            f.write(sbatch_content)
            
        bashCommand = "sbatch {:s}".format(sbatch_filename)
        os.system(bashCommand)

    def run_xi_chain(self,zbins=None,jobtype='debug',linear=False,yamlfn='/pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/mc/cobaya/yamls/fs_chain.yaml'):
        name = self.config['name']
        hod_string = get_hod_string(self.config['hod_params'])
        nbox = self.config['nbox']
        # outdir = 'hod/{:s}/{:s}/cobaya'.format(name,hod_string)
        outdir = os.path.join(self.basedir,'cobaya')

        newyamlfn = self.rewrite_xi_yaml(yamlfn,zbins=zbins,linear=linear,chain=True)
        newyamlfn = os.path.abspath(newyamlfn)

        job_name = "xi_{:s}".format(name)
        output_file = 'xi_ell' # os.path.join(outdir,'xi_ell')
        time = "00:30:00"  
        workdir = '/pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/mc/cobaya'
        sbatch_content = f"""#!/bin/bash -l
#SBATCH -J {job_name}
#SBATCH -N 1
#SBATCH -t {time}
#SBATCH -o {output_file}.out
#SBATCH -e {output_file}.err
#SBATCH -q {jobtype}
#SBATCH -C cpu
#SBATCH -A m68

module load python
source activate cobaya_07292024

export PYTHONPATH={workdir}:$PYTHONPATH

cd {workdir}
export OMP_NUM_THREADS=8
# srun -n 16 -c $OMP_NUM_THREADS cobaya-run -f {newyamlfn}
srun -n 16 -c $OMP_NUM_THREADS cobaya-run -f {newyamlfn}
"""
        
        sbatch_filename = "xi_chain.sh"
        with open(sbatch_filename, "w") as f:
            f.write(sbatch_content)
            
        bashCommand = "sbatch {:s}".format(sbatch_filename)
        os.system(bashCommand)

    ########
    # CHI2 #
    ########

    def chi2_wt(self,datwt):
        config = self.config
        name = config['name']
        hod_string = get_hod_string(config['hod_params'])
        outdir = os.path.join(self.basedir,'clustering')
        chi2 = 0
        for j, band in enumerate(m_bands):
            fn = os.path.join(outdir,'wt_%s_summary.txt'%band)
            dd = np.loadtxt(fn)
            rr = dd[:,0]
            wt = dd[:,1]
            dwt = dd[:,2]
            chi2 += np.sum((wt-datwt[j])**2/dwt**2)
            
        return chi2

    def calc_f2(self,i):
        # calculate f2 = 1/L for hod i
        config = self.config
        name = config['name']
        m_bands = config['bands']
        dat = self.get_hod_spec(i)
        band_chi = {
            'M411':[3821.6517, 4069.5162],
            'M438':[3998.3714, 4225.9277],
            'M464':[4160.4664, 4370.2586],
            'M490':[4309.7838, 4510.0708],
            'M517':[4447.8731, 4633.9157],
        }
        zmin, zmax = band_z['M411'][0],band_z['M517'][1]
        chimin, chimax = band_chi['M411'][0],band_chi['M517'][1]
        
        bin_width = 10
        chivec = np.linspace(chimin,chimax,100000)
        
        f2_dict = {}
        for j, band in enumerate(m_bands):
            dtmp_spec = dat[dat['SEL_%s'%band]]
            # ndensity = dtmp['RA'].size/mask.area()
            
            ### dNdchi
            zmin_tmp, zmax_tmp = band_z[band][0], band_z[band][1]
            chimin_tmp, chimax_tmp = band_chi[band][0], band_chi[band][1]
        
            bins = np.arange(chimin,chimax, bin_width) # spaced out by chi=10
            hchi = np.histogram(dtmp_spec['CHI'],bins = bins) # normalize
            x = (hchi[1][1:]+hchi[1][:-1])/2
            y = hchi[0]
            spl = make_smoothing_spline(x, y, lam=1e3)
            probabilities = spl(chivec)
            probabilities[probabilities<0] = 0
        
            # recent change
            probabilities[chivec<chimin_tmp] *= 0
            probabilities[chivec>chimax_tmp] *= 0
            
            # probabilities *= d['RA'].size/d_spec['RA'].size/6.11 * (1-fint_list[i])
            # dNdchi_band.append(probabilities)
        
            integral = np.trapz(probabilities,chivec)
            probabilities /= integral
            f2_dict[band] = np.trapz(probabilities**2,chivec)
        return f2_dict


    def chi2_wR(self,datwR,f2_dict=None):
        # computes chi2 for wp, not wR
        
        config = self.config
        name = config['name']
        m_bands = config['bands']
        hod_string = get_hod_string(config['hod_params'])
        outdir = os.path.join(self.basedir,'clustering')
        datarr, covarr = [], []

        if f2_dict==None: f2_dict = json.load(open('/pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/data_measurements/f2.json','r'))
        
        for j, band in enumerate(m_bands):
            f2 = f2_dict[band]
            assert f2>0
            datarr.append(np.loadtxt(os.path.join(outdir,'wR_{:s}.txt'.format(band))) / f2)
            covarr.append(np.loadtxt(os.path.join(outdir,'wR_{:s}_cov.txt'.format(band))) / f2**2 )
        datarr = np.array(datarr); covarr = np.array(covarr)
        datarr = np.average(datarr,axis=1)
        totwR, totcov = self.sum_inverse_variance(datarr,covarr)
        # assert np.allclose(totwR,np.loadtxt(os.path.join(outdir,'wR_tot.txt'))) # DO NOT load wR_tot

        tmp_datwR = datwR + 0
        for j, band in enumerate(m_bands):
            f2 = f2_dict[band]
            assert f2>0
            tmp_datwR[j] /= f2
        data_totwR, _ = self.sum_inverse_variance(tmp_datwR,covarr)
        
        diff = totwR-data_totwR
        chi2 = np.einsum('i,ij,j',diff,np.linalg.inv(totcov),diff)
            
        return chi2, totwR, data_totwR, totcov

'''    
    def chi2_wR(self,datwR):
        # computes chi2 for wp, not wR
        
        config = self.config
        name = config['name']
        m_bands = config['bands']
        hod_string = get_hod_string(config['hod_params'])
        nbox = config['nbox']
        outdir = os.path.join(self.basedir,'clustering')
        datarr, covarr = [], []

        # f2_dict = json.load(open('/pscratch/sd/h/hebina/AbacusLBG/ibis_tertiary44/LAE_auto_v2/data_measurements/f2.json','r'))
        for j, band in enumerate(m_bands):
            # f2 = f2_dict[band]
            # assert f2>0
            datarr.append(np.loadtxt(os.path.join(outdir,'wR_{:s}.txt'.format(band))))
            # covarr.append(np.loadtxt(os.path.join(outdir,'wR_{:s}_cov.txt'.format(band))))
        datarr = np.array(datarr)#; covarr = np.array(covarr)
        for i in range(nbox):
            f2_dict = self.calc_f2(i)
            for j, band in enumerate(m_bands):
                f2 = f2_dict[band]
                datarr[j,i] /= f2
                # covarr[j,i] /= f2**2
        covarr = [np.cov(datarr[j].T) for j, band in enumerate(m_bands)]
        covarr = np.array(covarr)
        datarr = np.average(datarr,axis=1)
        totwR, totcov = self.sum_inverse_variance(datarr,covarr)
        # assert np.allclose(totwR,np.loadtxt(os.path.join(outdir,'wR_tot.txt'))) # DO NOT load wR_tot

        tmp_datwR = datwR + 0
        for j, band in enumerate(m_bands):
            f2 = f2_dict[band]
            assert f2>0
            tmp_datwR[j] /= f2
        data_totwR, _ = self.sum_inverse_variance(tmp_datwR,covarr)
        
        diff = totwR-data_totwR
        chi2 = np.einsum('i,ij,j',diff,np.linalg.inv(totcov),diff)
            
        return chi2, totwR, data_totwR, totcov
'''



    ########
    # PLOT #
    ########

    # def plot_wp():
        






