#!/usr/bin/env python
#
#
# Generate a survey mask (an array of "filled" Healpix pixels)
# from a large (pre-generated) random file.  Write the result
# to a FITS file.
#
# Also contains a class that uses such a mask to return True
# for objects within the survey footprint/mask.
#
# Extinction correction coefficients for the ODIN NB filters
# a419 = 4.3238  # extinction correction for E(B-V)=1
# a501 = 3.54013 # extinction correction for E(B-V)=1
# a673 = 2.43846 # extinction correction for E(B-V)=1
# Magnitudes  ; reddening corrected mags
# m419=22.5-2.5*np.log10(od['flux_n419'])-a419*od['ebv']  # + 22.1069
# m501=22.5-2.5*np.log10(od['flux_n501'])-a501*od['ebv']  # + 22.1069
# m673=22.5-2.5*np.log10(od['flux_n673'])-a673*od['ebv']  # + 22.6585
#
#
import numpy  as np
import healpy as hp
import sys
import os


from   astropy.table import Table



def make_survey_mask(ran_fname,cut,nside=8192,is_nest=True):
    """Generate an astropy Table of the survey (inclusion)
    mask from the coordinates read from a file of randoms.
    We assume the file is large enough that Poisson fluctuations
    can be entirely neglected.
    An additional 'cut' is applied that all points lie within
    cut[2] degrees of RA/DEC=cut[0]/cut[1] (deg)."""
    # Read RA/DEC from the file, restrict to those with observations
    # in the desired filter and generate a list of Healpix pixel numbers.
    ran = Table.read(ran_fname) 
    '''
    Assume that circular cut (similar to data) has already been applied to the randoms
    in this case, radius = 1.4deg
    '''
    # ran = ran[ ran['MASKBITS']==0 ]
    # ran = ran[ ran['MASK']==0 ]

    # ran = ran[ ran['IN_ARJUN_MASK']==False ]
    # neededbands = ["N419","N501","N673"]
    # if filter_name=="N419": neededbands = ["N419","N501"]
    # if filter_name=="N501": neededbands = ["N501","N673"]
    # for band in neededbands:
    #     ran = ran[ ran['NOBS_'+band]>=10 ]
        # Now longer need ALLMASK_*, as this is included
        # in MASKBITS above.
        # ran = ran[ ran['ALLMASK_'+band]==0 ]
    # Apply the circle cut.
    if False:
        # A proper spherical cap cut centered on,
        # (cut[0],cut[1]) of radius cut[2] (all in degrees).
        cosmin= np.cos(np.radians(cut[2]))
        theta = np.radians(90-cut[1])
        phi   = np.radians(cut[0])
        cnhat = np.array([np.sin(theta)*np.cos(phi),\
                          np.sin(theta)*np.sin(phi),\
                          np.cos(theta)])
        theta = np.radians(90.-ran['DEC'])
        phi   = np.radians(ran['RA'])
        nhat  = np.zeros( (theta.size,3) )
        nhat[:,0] = np.sin(theta)*np.cos(phi)
        nhat[:,1] = np.sin(theta)*np.sin(phi)
        nhat[:,2] = np.cos(theta)
        cdist = np.dot(nhat,cnhat)
        theta = theta[ cdist>cosmin ]
        phi   = phi[   cdist>cosmin ]
    else:
        # Arjun actually applies a simpler "tangent" cut,
        # so implement that instead.
        rad2  = (ran['RA' ]-cut[0])**2 +\
                (ran['DEC']-cut[1])**2
        theta = np.radians(90.-ran['DEC'])
        phi   = np.radians(ran['RA'])
        print(theta.shape)
        if cut[2]<179.999:
            ww    = rad2<cut[2]**2
            theta = theta[ww]
            phi   = phi[ww]
        print(theta.shape)
    # and work out the HEALPIX pixel numbers:
    pixs  = hp.ang2pix(nside,theta,phi,nest=is_nest)
    # and generate a table containing the information to be returned.
    tt = Table({"HPXPIXEL":np.unique(pixs)})
    tt.meta["HPXNSID"] = nside
    tt.meta["HPXNEST"] = is_nest
    tt.meta["COMMENT"] = r'Mask lists HPX pixels included in survey'
    # We also want a random file -- this could be downsampled.
    rr = Table({'RA':ran['RA'][::10],'DEC':ran['DEC'][::10]})
    return( (tt,rr) )
    #






class SurveyMask:
    def __init__(self,maskfn,wtfn=None):
        """Reads a FITS file containing an inclusion mask and an
        optional 'weights' map."""
        mskd       = Table.read(maskfn)
        self.nside = mskd.meta["HPXNSID"]
        self.nest  = mskd.meta["HPXNEST"]
        self.pixs  = mskd["HPXPIXEL"]
        if 'MASK' in mskd.keys(): self.pixs = mskd["HPXPIXEL"][mskd['MASK']==0]
        
        if wtfn is None:
            self.wtmask = None
        else:
            self.wtmask = Table.read(wtfn)
    def area(self):
        """Returns the in_mask area in deg2."""
        apix = hp.nside2pixarea(self.nside,degrees=True)
        return(len(self.pixs)*apix)
    def weights(self,ras,decs):
        """Returns the weights for the points given by (RA,DEC) in degrees."""
        if self.wtmask is None:
            wt   = np.ones_like(ras)
        else:
            nside= self.wtmask.meta["HPXNSIDE"]
            nest = self.wtmask.meta["HPXNEST"]
            tt   = np.radians(90.-decs)
            pp   = np.radians(ras)
            pixs = hp.ang2pix(nside,tt,pp,nest=nest)
            wt   = np.zeros_like(ras)
            for pixnum,wval in \
                zip(self.wtmask["HPXPIXEL"],\
                    self.wtmask["FRACFLUX_N501_FRACSEL"]):
                k = np.nonzero(pixs==pixnum)[0]
                if k.size>0: wt[k] = wval
        return(wt)
    def __call__(self,ras,decs):
        """Returns a boolean array of whether the points pass the mask,
        with the points given by (RA,DEC) in degrees."""
        tt   = np.radians(90.-decs)
        pp   = np.radians(ras)
        pixs = hp.ang2pix(self.nside,tt,pp,nest=self.nest)
        return(np.in1d(pixs,self.pixs))
        #




if __name__=="__main__":
    # if len(sys.argv)!=3:
    #     raise RuntimeError("Usage: "+sys.argv[0]+" <ran-fname> <filter_name>")
    # ran_fname   = sys.argv[1]
    ran_fname = 'ibis-xmm.ran.fits'
    
    # filter_name = sys.argv[2].upper()
    #
    cut = [35.75, -4.75,180.]
    '''
    Don't want spherical cuts
    '''
    # if ran_fname.find('2band')>0:
    #     cut = [150.11,2.173,1.9]
    # if ran_fname.find('cosmos')>0:
    #     cut = [150.11,2.173,1.9]
    #
    tt,rr = make_survey_mask(ran_fname,cut)
    # tt,rr = make_survey_mask(ran_fname,cut,nside=4096)
    #
    tt.write('ibis_tertiary44_msk.fits',overwrite=True)
    rr.write("ibis_tertiary44_ran.fits",overwrite=True)
    #
