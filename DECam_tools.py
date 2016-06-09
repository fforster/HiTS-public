import numpy as np
from scipy.special import erf

class DECam_tools(object):

    def __init__(self, hitsdir):
        
        # open CCD numbers file
        self.hitsdir = hitsdir
        self.CCDn = {}
        (CCDstring, CCDnumber) = np.loadtxt("%s/etc/CCDnumbers.dat" % self.hitsdir, dtype = str).transpose()
        CCDnumber = np.array(CCDnumber, dtype = int)
        for i in range(len(CCDstring)):
            self.CCDn[CCDstring[i]] = CCDnumber[i]
            
        # open zero point file given the filtername
        self.IDzero = {}
        self.azero = {}
        self.e_azero = {}
        self.bzero = {}
        self.e_bzero = {}
        self.kzero = {}
        self.e_kzero = {}
        
        for filtername in ['u', 'g', 'r', 'i', 'z', 'Y']:
            (IDzeroi, filterzeroi, azeroi, e_azeroi, bzeroi, e_bzeroi, kzeroi, e_kzeroi) = np.loadtxt("%s/etc/zeropoints_%s.txt" % (self.hitsdir, filtername), dtype = str).transpose()
            self.IDzero[filtername] = np.array(IDzeroi, dtype = int)
            self.azero[filtername] = np.array(azeroi, dtype = float)
            self.e_azero[filtername] = np.array(e_azeroi, dtype = float)
            self.bzero[filtername] = np.array(bzeroi, dtype = float)
            self.e_bzero[filtername] = np.array(e_bzeroi, dtype = float)
            self.kzero[filtername] = np.array(kzeroi, dtype = float)
            self.e_kzero[filtername] = np.array(e_kzeroi, dtype = float)

    # function to convert fluxes into magnitudes given fluxes and errors in ADU, the CCD number, the exposure time and the airmass of the observation
    def ADU2mag(self, flux, e_flux, CCD, exptime, airmass, filtername):
        
        mag = np.ones(np.shape(flux)) * 30
        mag_1 = np.ones(np.shape(flux)) * 30
        mag_2 = np.ones(np.shape(flux)) * 30
        fluxp = flux + e_flux
        fluxm = flux - e_flux
        mflux = (flux > 0)
        mfluxp = (fluxp > 0)
        mfluxm = (fluxm > 0)
        mag[mflux] = np.array(-2.5 * np.log10(flux[mflux]) + 2.5 * np.log10(exptime) - self.azero[filtername][self.CCDn[CCD] - 1] - self.kzero[filtername][self.CCDn[CCD] - 1] * airmass)
        mag_1[mfluxp] = np.array(-2.5 * np.log10(fluxp[mfluxp]) + 2.5 * np.log10(exptime) - self.azero[filtername][self.CCDn[CCD] - 1] - self.kzero[filtername][self.CCDn[CCD] - 1] * airmass)
        mag_2[mfluxm] = np.array(-2.5 * np.log10(fluxm[mfluxm]) + 2.5 * np.log10(exptime) - self.azero[filtername][self.CCDn[CCD] - 1] - self.kzero[filtername][self.CCDn[CCD] - 1] * airmass)
        return (mag, mag - mag_1, mag_2 - mag)


    # function to convert fluxes into magnitudes for typical CCD
    def ADU2mag_avg(self, flux, exptime, airmass, filtername):
        
        mag = np.array(-2.5 * np.log10(flux) + 2.5 * np.log10(exptime) - np.average(self.azero[filtername]) - np.average(self.kzero[filtername]) * airmass)
        return mag

    # analytic form for efficiency vs ADU relation
    def efficiency_ADU(self, xs, offset, delta):
        return (1. + erf((log10(xs) - offset) / delta)) / 2.
        
    # analytic form for efficiency vs mag relation
    def efficiency_mag(self, xs, offset, delta):
        return (1. + erf(-(xs - offset) / delta)) / 2.

