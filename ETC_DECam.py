import numpy as np

# CC F. Forster, J. Martinez, J.C. Maureira, modified public DECam exposure time calculator
# -------------------------------------------------------------
# modified to include the effect of airmass and able to use empirical FWHM and/or sky brightness
# ETC_DECam object constructor: requires seeing in r band in arcsec (kw. seeing_r_arcsec, default if 0.75) and ETC version (5 or 6 accepted, default is 6, kw. vETC), it can also run in test LSST or HSC modes
# Available class functions are FWHM, SNR and findmag:
#   FWHM: compute FWHM given band and airmass
#   SNR: compute SNR given source magnitude (kw. mag), band (kw. mag), airmass (kw. airmass), exposure time (kw. exptime), sky properties and optional FWHM (fwhm)
#        skymode = mag: sky given in mag per arcsec2 at zenith (skymag), FWHM derived from seeing in r band, band and airmass
#        skymode = ADU: use empirical sky in ADU (skyADU), FWHM derived from seeing in r band, band and airmass
#        skymode = mag-FWHM: sky given in mag per arcsec2 at zenith (skymag), use empirical FWHM (fwhm)
#        skymode = ADU-FWHM: use empirical sky in ADU (skyADU), use empirical FWHM (fwhm)
#   findmag: find magnitude for given SNR given, same parameters as SNR, expect for source magnitude (not needed) and input SNR (SNRin)

class ETC_DECam(object):

    # class constructor given seeing in r band in arcsec and ETC version as keyword arguments 'seeing_r_arcsec' and 'vETC'
    def __init__(self, **kwargs):
        
        if 'seeing_r_arcsec' in kwargs.keys():
            self.seeing_r_arcsec = kwargs['seeing_r_arcsec']
        else:
            self.seeing_r_arcsec = 0.75 # median value at CTIO

        if 'vETC' in kwargs.keys():
            self.vETC = kwargs['vETC']
        else:
            self.vETC = 6 # recommended version

        # test parameters (default to false)
        self.testLSST = False
        self.testHSC = False
        if 'testLSST' in kwargs.keys():
            self.testLSST = kwargs['testLSST']
            if self.testLSST:
                print "ETC LSST mode ON (larger area and smaller pixel size)"
        if 'testHSC' in kwargs.keys():
            self.testHSC = kwargs['testHSC']
            if self.testHSC:
                print "ETC HSC mode ON (larger area and smaller pixel size)"
            if self.testHSC and self.testLSST:
                print "Both LSST and HSC test modes cannot be ON, setting HSC to OFF"
                self.testHSC = False

    # compute FWHM given band and airmass as keyword arguments 'band' and 'airmass'
    def FWHM(self, **kwargs):

        band = kwargs['band']
        
        # DECam filters central wavelengths
        cwl_u_nm = 375 # nm
        cwl_g_nm = 473.5 # nm
        cwl_r_nm = 638.5 # nm
        cwl_i_nm = 775.5 # nm
        cwl_z_nm = 922.5 # nm
        cwl_Y_nm = 995 # nm
        
        # scaling between seeing at different filters
        seeing_g_arcsec = self.seeing_r_arcsec * (cwl_r_nm / cwl_g_nm)**0.2
        seeing_u_arcsec = 0.2 + seeing_g_arcsec * (cwl_g_nm / cwl_u_nm)**0.2
        seeing_i_arcsec = self.seeing_r_arcsec * (cwl_r_nm / cwl_i_nm)**0.2
        seeing_z_arcsec = seeing_i_arcsec * (cwl_i_nm / cwl_z_nm)**0.2
        seeing_Y_arcsec = seeing_z_arcsec * (cwl_z_nm / cwl_Y_nm)**0.2
    
        # select seeing depending on band
        if band == 'u':
            seeing_arcsec = seeing_u_arcsec
        elif band == 'g':
            seeing_arcsec = seeing_g_arcsec
        elif band == 'r':
            seeing_arcsec = self.seeing_r_arcsec
        elif band == 'i':
            seeing_arcsec = seeing_i_arcsec
        elif band == 'z':
            seeing_arcsec = seeing_z_arcsec
        elif band == 'Y':
            seeing_arcsec = seeing_Y_arcsec
        
        fseeing_airmass = 1. / np.cos(np.arccos(1. / kwargs['airmass']))**(3./5.) # seeing correction factor due to airmass

        FWHM_arcsec = np.sqrt((seeing_arcsec * fseeing_airmass)**2 + 0.63**2)
        
        return FWHM_arcsec

    # compute SNR given band, source magnitude, exposure time and airmass
    # different modes for sky/FWHM calculation
    # skymode = mag: sky given in mag per arcsec2 at zenith (skymag), FWHM derived from seeing in r band, band and airmass
    # skymode = ADU: use empirical sky in ADU (skyADU), FWHM derived from seeing in r band, band and airmass
    # skymode = mag-FWHM: sky given in mag per arcsec2 at zenith (skymag), use empirical FWHM (fwhm)
    # skymode = ADU-FWHM: use empirical sky in ADU (skyADU), use empirical FWHM (fwhm)
    def SNR(self, **kwargs):#, mag, exptime, airmass, skymode, skylevel):
    
        # non-optional arguments
        band = kwargs['band']
        mag = kwargs['mag']
        exptime = kwargs['exptime']
        airmass = kwargs['airmass']
        skymode = kwargs['skymode']
        if (skymode == 'mag' or skymode == 'mag-FWHM') and not 'skymag' in kwargs.keys():
            raise ValueError("Missing sky magnitude (skymag)")
        elif (skymode == 'ADU' or skymode == 'ADU-FWHM')  and not 'skyADU' in kwargs.keys():
            raise ValueError("Missing sky level in ADU (skyADU)")
        elif (skymode == 'ADU-FWHM' or skymode == 'mag-FWHM') and not 'fwhm' in kwargs.keys():
            raise ValueError("Missing FWHM (fwhm)")

        # optional arguments: skymag or skyADU, fwhm
        
        # speed of light
        cspeed = 2.99792458e17 # nm
    
        # DECam filters central wavelengths
        cwl_u_nm = 375 # nm
        cwl_g_nm = 473.5 # nm
        cwl_r_nm = 638.5 # nm
        cwl_i_nm = 775.5 # nm
        cwl_z_nm = 922.5 # nm
        cwl_Y_nm = 995 # nm
        
        # scaling between seeing at different filters
        seeing_g_arcsec = self.seeing_r_arcsec * (cwl_r_nm / cwl_g_nm)**0.2
        seeing_u_arcsec = 0.2 + seeing_g_arcsec * (cwl_g_nm / cwl_u_nm)**0.2
        seeing_i_arcsec = self.seeing_r_arcsec * (cwl_r_nm / cwl_i_nm)**0.2
        seeing_z_arcsec = seeing_i_arcsec * (cwl_i_nm / cwl_z_nm)**0.2
        seeing_Y_arcsec = seeing_z_arcsec * (cwl_z_nm / cwl_Y_nm)**0.2
    
        # select seeing depending on band
        if band == 'u':
            bandpass_nm = 50 # nm
            cwl_nm = cwl_u_nm
            seeing_arcsec = seeing_u_arcsec
        elif band == 'g':
            bandpass_nm = 147 # nm
            cwl_nm = cwl_g_nm
            seeing_arcsec = seeing_g_arcsec
        elif band == 'r':
            bandpass_nm = 141 # nm
            cwl_nm = cwl_r_nm
            seeing_arcsec = self.seeing_r_arcsec
        elif band == 'i':
            bandpass_nm = 147 # nm
            cwl_nm = cwl_i_nm
            seeing_arcsec = seeing_i_arcsec
        elif band == 'z':
            bandpass_nm = 147 # nm
            cwl_nm = cwl_z_nm
            seeing_arcsec = seeing_z_arcsec
        elif band == 'Y':
            bandpass_nm = 50 # nm
            cwl_nm = cwl_Y_nm
            seeing_arcsec = seeing_Y_arcsec
    
        # filter width in Hz
        bandpass_Hz = cspeed / (cwl_nm - bandpass_nm / 2.) - cspeed / (cwl_nm + bandpass_nm / 2.)  # Hz
    
        # Planck's constant in MKS
        hPlanck_MKS = 6.62606957e-34 # m2 kg / s
        # hc / lambda in MKS (note cspeed and cwl_nm are in nm)
        hc_lambda = hPlanck_MKS * cspeed / cwl_nm
    
        # magnitude scale zero point flux
        zero_photrate = 3.631e-23 * bandpass_Hz / hc_lambda # photons / sec / m2  #  taken from ETC
    
        # primary mirror effective area
        area = 9.7 # m2 
        if self.testLSST:
            area = 35.04
        elif self.testHSC:
            area = 53.0

        # photons from zero magnitude source
        zero_phot = zero_photrate * area * exptime  # photons
    
        # primary mirror reflectivity
        if self.vETC == 6:
            aperture_eff = 2.04 
            if band == 'u':
                CCD_eff = 0.25
                primary_refl = 0.89
            elif band == 'g':
                CCD_eff = 0.59
                primary_refl = 0.89
            elif band == 'r':
                CCD_eff = 0.75
                primary_refl = 0.88
            elif band == 'i':
                CCD_eff = 0.85
                primary_refl = 0.87
            elif band == 'z':
                CCD_eff = 0.85
                primary_refl = 0.88
            elif band == 'Y':
                CCD_eff = 0.5
                primary_refl = 0.9
        elif self.vETC == 5:
            aperture_eff = 1.34
            if band == 'u':
                CCD_eff = 0.25
                primary_refl = 0.85
            elif band == 'g':
                CCD_eff = 0.7
                primary_refl = 0.93
            elif band == 'r':
                CCD_eff = 0.75
                primary_refl = 0.9
            elif band == 'i':
                CCD_eff = 0.85
                primary_refl = 0.87
            elif band == 'z':
                CCD_eff = 0.8
                primary_refl = 0.93
            elif band == 'Y':
                CCD_eff = 0.3
                primary_refl = 0.95
        else:
            raise ValueError("Please define an ETC version before calling SNR. EXIT")
    
        # vignetting
        vig = 1.0
    
        # filter, corrector and atmospheric transmission
        if band == 'u':
            filter_eff = 0.95
            corrector_eff = 0.75
            atmosph_eff = 0.7 * np.exp(1.) / np.exp(airmass)
        elif band == 'g':
            filter_eff = 0.9
            corrector_eff = 0.86
            atmosph_eff = 0.8 * np.exp(1.) / np.exp(airmass)
        elif band == 'r':
            filter_eff = 0.9
            corrector_eff = 0.86
            atmosph_eff = 0.9 * np.exp(1.) / np.exp(airmass)
        elif band == 'i':
            filter_eff = 0.9
            corrector_eff = 0.86
            atmosph_eff = 0.9 * np.exp(1.) / np.exp(airmass)
        elif band == 'z':
            filter_eff = 0.9
            corrector_eff = 0.86
            atmosph_eff = 0.9 * np.exp(1.) / np.exp(airmass)
        elif band == 'Y':
            filter_eff = 0.9
            corrector_eff = 0.75
            atmosph_eff = 0.95 * np.exp(1.) / np.exp(airmass)
    
            
        # derive FWHM from seeing and airmass or use empirical FWHM
        if skymode == 'mag' or skymode == 'ADU':
            fseeing_airmass = 1. / np.cos(np.arccos(1. / airmass))**(3./5.) # seeing correction factor due to airmass
            # airmass and optics effect on FWHM
            FWHM_arcsec = np.sqrt((seeing_arcsec * fseeing_airmass)**2 + 0.63**2)
        elif skymode == 'mag-FWHM' or skymode == 'ADU-FWHM':
            # if FWHM is provided do not scale seeing
            FWHM_arcsec = kwargs['fwhm']
        else:
            raise ValueError("SNR: skymode %s not recognized" % skymode)
    
        # aperture
        aperture = np.pi * ((aperture_eff * FWHM_arcsec) / 2.)**2  # aperture in arcsec2
    
        # final throughput
        throughput = CCD_eff * filter_eff * corrector_eff * primary_refl * vig * atmosph_eff

        # signal from zero magnitude source
        zero_signal = zero_phot * throughput # electrons from a zero mag source

        # electrons from source
        source_electrons = zero_signal * 10**(mag / -2.5) # electrons from a source of the given magnitude
    
        # DECam pixel scale
        DECam_pix = 0.263 # arcsec  0.2637 in the centre, 0.2626 at the edges
        if self.testLSST:
            DECam_pix = 0.2 # arcsec
        if self.testHSC:
            DECam_pix = 0.16 # arcsec

        # DECam typical gain
        gain = 4.0
    
        # electrons from the sky: 1) sky_ADU is given (use gain, aperture, pixel scale), 2) sky_mag is given (use zero_signal, aperture and airmass)
        if skymode == 'ADU' or skymode == 'ADU-FWHM':
            sky_electrons = kwargs['skyADU'] * gain * aperture / DECam_pix**2 # electrons from the sky per aperture given the empirical sky per pixel in ADU
        elif skymode == 'mag' or skymode == 'mag-FWHM':
            sky_electrons = 10**(-kwargs['skymag'] / 2.5) * zero_signal * aperture * airmass # electrons from the sky per aperture
        else:
            raise ValueError("SNR: skymode %s not recognized" % skymode)
    
        # readout noise per pixel
        RON_pix = 10. # electrons

        # readout noise per aperture
        RON_aper = RON_pix * np.sqrt(aperture / DECam_pix**2) # electrons
    
        # surce signal to noise ratio
        SNRout = source_electrons / np.sqrt(source_electrons + sky_electrons + RON_aper)
    
        return SNRout
    
    # find magnitude given target SNR (very stupid search for now)
    def findmag(self, **kwargs):
        
        magsarray = np.linspace(15., 27., 100000)
        if kwargs['skymode'] == 'mag':
            SNRs = self.SNR(band=kwargs['band'], mag=magsarray, exptime=kwargs['exptime'], airmass=kwargs['airmass'], skymode=kwargs['skymode'], skymag=kwargs['skymag'])
        elif kwargs['skymode'] == 'mag-FWHM':
            SNRs = self.SNR(band=kwargs['band'], mag=magsarray, exptime=kwargs['exptime'], airmass=kwargs['airmass'], skymode=kwargs['skymode'], skymag=kwargs['skymag'], fwhm=kwargs['fwhm'])
        elif kwargs['skymode'] == 'ADU':
            SNRs = self.SNR(band=kwargs['band'], mag=magsarray, exptime=kwargs['exptime'], airmass=kwargs['airmass'], skymode=kwargs['skymode'], skyADU=kwargs['skyADU'])
        elif kwargs['skymode'] == 'ADU-FWHM':
            SNRs = self.SNR(band=kwargs['band'], mag=magsarray, exptime=kwargs['exptime'], airmass=kwargs['airmass'], skymode=kwargs['skymode'], skyADU=kwargs['skyADU'], fwhm=kwargs['fwhm'])
        else:
            raise ValueError("findmag: wrong keyword arguments")

        if 'SNRin' not in kwargs.keys():
            raise ValueError("findmag: missing input SNR")
        else:
            return magsarray[np.argmin((SNRs - kwargs['SNRin'])**2)]
    
if __name__ == "__main__":

    print "\nInitializing ETC_DECam object with seeing_r=of 0.75\", version 6..."
    ETC = ETC_DECam(seeing_r_arcsec=0.75, vETC=6)

    print "\nTesting FWHM for an airmass vector..."
    print "u band, airmass between 1.0 and 1.6", ETC.FWHM(band='u', airmass=np.linspace(1.0, 1.6, 10))
    print "g band, airmass between 1.2 and 1.9", ETC.FWHM(band='g', airmass=np.linspace(1.2, 1.9, 10))

    print "\nTesting SNR and findmag with all skymodes with a 20 mag source, 173 sec exposure time in g band, airmass of 1.0..."
    SNRtest = ETC.SNR(band='g', mag=20, exptime=173, airmass=1.0, skymode='mag', skymag=22)
    magtest = ETC.findmag(band='g', SNRin=SNRtest, exptime=173, airmass=1.0, skymode='mag', skymag=22.0)
    print "mag (skymag=22): SNR %f <-> mag %f" % (SNRtest, magtest)
    if np.abs(20. - magtest) > 1e-4:
        print "mag not OK"
    else:
        print "OK"
    SNRtest = ETC.SNR(band='g', mag=20, exptime=173, airmass=1.0, skymode='mag-FWHM', skymag=22, fwhm=1.0)
    magtest = ETC.findmag(band='g', SNRin=SNRtest, exptime=173, airmass=1.0, skymode='mag-FWHM', skymag=22, fwhm=1.0)
    print "mag-FWHM (skymag=22, fwhm=1): SNR %f <-> mag %f"  % (SNRtest, magtest)
    if np.abs(20. - magtest) > 1e-4:
        print "mag-FWHM not OK"
    else:
        print "OK"
    SNRtest = ETC.SNR(band='g', mag=20, exptime=173, airmass=1.0, skymode='ADU', skyADU=120)
    magtest = ETC.findmag(band='g', SNRin=SNRtest, exptime=173, airmass=1.0, skymode='ADU', skyADU=120)
    print "ADU (skyADU=120): SNR %f <-> mag %f" % (SNRtest, magtest)
    if np.abs(20. - magtest) > 1e-4:
        print "ADU not OK"
    else:
        print "OK"
    SNRtest = ETC.SNR(band='g', mag=20, exptime=173, airmass=1.0, skymode='ADU-FWHM', skyADU=120, fwhm=1.0)
    magtest = ETC.findmag(band='g', SNRin=SNRtest, exptime=173, airmass=1.0, skymode='ADU-FWHM', skyADU=120, fwhm=1.0)
    print "ADU-FHWM (skyADU=120, fwhm=1): SNR %f <-> mag %f" % (SNRtest, magtest)
    if np.abs(20. - magtest) > 1e-4:
        print "ADU-FWHM not OK"
    else:
        print "OK"

    print "\nTesting findmag for input SNR of 5..."
    print "u band, mags at airmasses of 1.0, 1.2, 1.4, 1.6, 1.8, 2.0", map(lambda x: ETC.findmag(band='u', SNRin=5, exptime=173, airmass=x, skymode='mag', skymag=22.0), [1., 1.2, 1.4, 1.6, 1.8, 2.0])
    print "g band, mags at airmasses of 1.0, 1.2, 1.4, 1.6, 1.8, 2.0", map(lambda x: ETC.findmag(band='g', SNRin=5, exptime=173, airmass=x, skymode='mag', skymag=22.0), [1., 1.2, 1.4, 1.6, 1.8, 2.0])
    print "r band, mags at airmasses of 1.0, 1.2, 1.4, 1.6, 1.8, 2.0", map(lambda x: ETC.findmag(band='r', SNRin=5, exptime=173, airmass=x, skymode='mag', skymag=22.0), [1., 1.2, 1.4, 1.6, 1.8, 2.0])
    print "i band, mags at airmasses of 1.0, 1.2, 1.4, 1.6, 1.8, 2.0", map(lambda x: ETC.findmag(band='i', SNRin=5, exptime=173, airmass=x, skymode='mag', skymag=22.0), [1., 1.2, 1.4, 1.6, 1.8, 2.0])
    print "z band, mags at airmasses of 1.0, 1.2, 1.4, 1.6, 1.8, 2.0", map(lambda x: ETC.findmag(band='z', SNRin=5, exptime=173, airmass=x, skymode='mag', skymag=22.0), [1., 1.2, 1.4, 1.6, 1.8, 2.0])
    print "Y band, mags at airmasses of 1.0, 1.2, 1.4, 1.6, 1.8, 2.0", map(lambda x: ETC.findmag(band='Y', SNRin=5, exptime=173, airmass=x, skymode='mag', skymag=22.0), [1., 1.2, 1.4, 1.6, 1.8, 2.0])

