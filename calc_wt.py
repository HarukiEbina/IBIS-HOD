#!/usr/bin/env python3
#
# Uses Corrfunc to compute the angular correlation function
# of a sample (provided data and randoms).
#
#
import numpy as np
import os

from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks





def calc_wt(dat,ran,bins=None):
    """Computes the (L-S) angular correlation function"""
    # Get number of threads, datas and randoms.
    nt   = int(os.getenv('OMP_NUM_THREADS',1))
    pp   = 'pair_product'
    # Bin edges are specified in Mpc/h, if nothing
    # is passed in, do log-spaced bins.
    if bins is None:
        Nbin = 15
        bb   = np.geomspace(0.005,0.5,Nbin+1)
    else:
        bb   = bins
    # RA and DEC should be in degrees, and all arrays should
    # be the same type, and it seems they need to be float
    # and not e.g. float32.  Ensure this now.
    dra = np.ascontiguousarray(dat['RA' ]).astype('float')
    ddc = np.ascontiguousarray(dat['DEC']).astype('float')
    rra = np.ascontiguousarray(ran['RA' ]).astype('float')
    rdc = np.ascontiguousarray(ran['DEC']).astype('float')
    if 'WT' in dat.keys():
        dwt = np.ascontiguousarray(dat['WT']).astype('float')
    else:
        dwt = np.ones_like(dat['RA']).astype('float')
    if 'WT' in ran.keys():
        rwt = np.ascontiguousarray(ran['WT']).astype('float')
    else:
        rwt = np.ones_like(ran['RA']).astype('float')
    # Get the weight sums and the pair counts.
    nd,nr= np.sum(dwt),np.sum(rwt)
    DD_c = DDtheta_mocks(1,nt,bb,RA1=dra,DEC1=ddc,weights1=dwt,weight_type=pp)
    RR_c = DDtheta_mocks(1,nt,bb,RA1=rra,DEC1=rdc,weights1=rwt,weight_type=pp)
    DR_c = DDtheta_mocks(0,nt,bb,RA1=dra,DEC1=ddc,weights1=dwt,\
                                 RA2=rra,DEC2=rdc,weights2=rwt,weight_type=pp)
    # Compute the normalized, weighted pair counts.
    dd   = DD_c['npairs']*DD_c['weightavg']/nd/nd
    dr   = DR_c['npairs']*DR_c['weightavg']/nd/nr
    rr   = RR_c['npairs']*RR_c['weightavg']/nr/nr
    # and hence the angular correlation function.
    tt,wt= np.sqrt(bb[:-1]*bb[1:]),(dd-2*dr+rr)/(rr+1e-30)
    return((tt,wt))
    #

def calc_wt_cross(dat1,dat2,ran,bins=None):
    """Computes the (L-S) angular correlation function"""
    # formally we need ran1 and ran2, but ignore for now
    
    # Get number of threads, datas and randoms.
    nt   = int(os.getenv('OMP_NUM_THREADS',1))
    pp   = 'pair_product'
    # Bin edges are specified in Mpc/h, if nothing
    # is passed in, do log-spaced bins.
    if bins is None:
        Nbin = 15
        bb   = np.geomspace(0.005,0.5,Nbin+1)
    else:
        bb   = bins
    # RA and DEC should be in degrees, and all arrays should
    # be the same type, and it seems they need to be float
    # and not e.g. float32.  Ensure this now.
    d1ra = np.ascontiguousarray(dat1['RA' ]).astype('float')
    d1dc = np.ascontiguousarray(dat1['DEC']).astype('float')
    d2ra = np.ascontiguousarray(dat2['RA' ]).astype('float')
    d2dc = np.ascontiguousarray(dat2['DEC']).astype('float')
    rra = np.ascontiguousarray(ran['RA' ]).astype('float')
    rdc = np.ascontiguousarray(ran['DEC']).astype('float')
    if 'WT' in dat1.keys():
        d1wt = np.ascontiguousarray(dat1['WT']).astype('float')
    else:
        d1wt = np.ones_like(dat1['RA']).astype('float')
    if 'WT' in dat2.keys():
        d2wt = np.ascontiguousarray(dat2['WT']).astype('float')
    else:
        d2wt = np.ones_like(dat2['RA']).astype('float')
    if 'WT' in ran.keys():
        rwt = np.ascontiguousarray(ran['WT']).astype('float')
    else:
        rwt = np.ones_like(ran['RA']).astype('float')
    # Get the weight sums and the pair counts.
    nd1,nd2,nr= np.sum(d1wt),np.sum(d2wt),np.sum(rwt)
    DD_c = DDtheta_mocks(0,nt,bb,RA1=d1ra,DEC1=d1dc,weights1=d1wt,RA2=d2ra,DEC2=d2dc,weights2=d2wt,weight_type=pp)
    RR_c = DDtheta_mocks(1,nt,bb,RA1=rra,DEC1=rdc,weights1=rwt,weight_type=pp)
    D1R_c = DDtheta_mocks(0,nt,bb,RA1=d1ra,DEC1=d1dc,weights1=d1wt,\
                                 RA2=rra,DEC2=rdc,weights2=rwt,weight_type=pp)
    D2R_c = DDtheta_mocks(0,nt,bb,RA1=d2ra,DEC1=d2dc,weights1=d2wt,\
                                 RA2=rra,DEC2=rdc,weights2=rwt,weight_type=pp)
    # Compute the normalized, weighted pair counts.
    dd   = DD_c['npairs']*DD_c['weightavg']/nd1/nd2
    d1r   = D1R_c['npairs']*D1R_c['weightavg']/nd1/nr
    d2r   = D2R_c['npairs']*D2R_c['weightavg']/nd2/nr
    rr   = RR_c['npairs']*RR_c['weightavg']/nr/nr
    # and hence the angular correlation function.
    tt,wt= np.sqrt(bb[:-1]*bb[1:]),(dd-d1r-d2r+rr)/(rr+1e-30)
    return((tt,wt))
    #



