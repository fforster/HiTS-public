import numpy as np
import sys


# UNDER CONSTRUCTION

class ETC_DECam(object):


    def __init__(self, **kwargs):
        
        self.vETC = kwargs['vETC']
        self.seeing_r_arcsec = kwargs['seeing_r_arcsec']
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
                print "Both LSST and HSC test modes cannot be ON, exit."
                sys.exit()

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

    # compute SNR given magnitude (vETC is ETC version, e.g. 5 or 6)
    def SNR(self, **kwargs):#, mag, exptime, airmass, skymode, skylevel):
    
        # compulsory arguments
        band = kwargs['band']
        mag = kwargs['mag']
        exptime = kwargs['exptime']
        airmass = kwargs['airmass']
        skymode = kwargs['skymode']
        if (skymode == 'mag' or skymode == 'mag-FWHM') and not 'skymag' in kwargs.keys():
            print "Missing sky magnitude (skymag)"
            sys.exit()
        elif (skymode == 'ADU' or skymode == 'mag-FWHM') and not 'fwhm' in kwargs.keys():
            print "Missing fwhm (fwhm)"
            sys.exit()
        elif skymode == 'ADU' and not 'skyADU' in kwargs.keys():
            print "Missing sky level in ADU (skyADU)"
            sys.exit()

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
    
        zero_photrate = 3.631e-23 * bandpass_Hz / hc_lambda # photons / sec / m2  #  where does the number come from?
    
        # primary mirror effective area
        area = 9.7 # m2 35.04 LSST
        # LSST test
        if self.testLSST:
            area = 35.04
        elif self.testHSC:
            area = 53.0
    
        zero_phot = zero_photrate * area * exptime  # photons
    
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
            print "Please define an ETC version before calling SNR. EXIT"
            sys.exit()
    
        vig = 1.0
    
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
    
            
        # how to scale seeing with airmass
        if skymode == 'mag':
            fseeing_airmass = 1. / np.cos(np.arccos(1. / airmass))**(3./5.) # seeing correction factor due to airmass
            # airmass and optics effect on FWHM
            FWHM_arcsec = np.sqrt((seeing_arcsec * fseeing_airmass)**2 + 0.63**2)
        else:
            # if FWHM is provided do not scale seeing
            FWHM_arcsec = seeing_arcsec
    
        # aperture
        aperture = np.pi * ((aperture_eff * FWHM_arcsec) / 2.)**2  # aperture in arcsec2
    
        throughput = CCD_eff * filter_eff * corrector_eff * primary_refl * vig * atmosph_eff
        zero_signal = zero_phot * throughput # electrons from a zero mag source
        source_electrons = zero_signal * 10**(mag / -2.5) # electrons from a source of the given magnitude
    
        DECam_pix = 0.263 # arcsec  0.2637 in the centre, 0.2626 at the edges
        gain = 4.0
        # test LSST
        if self.testLSST:
            DECam_pix = 0.2 # arcsec
        if self.testHSC:
            DECam_pix = 0.16 # arcsec
    
        # electrons from the sky: 1) sky_ADU is given (use gain, aperture, pixel scale), 2) sky_mag is given (use zero_signal, aperture and airmass)
        if skymode == 'ADU':
            sky_electrons = kwargs['skyADU'] * gain * aperture / DECam_pix**2 # electrons from the sky per aperture given the empirical sky per pixel in ADU
        elif skymode == 'mag' or skymode == 'mag-FWHM':
            sky_electrons = 10**(-kwargs['skymag'] / 2.5) * zero_signal * aperture * airmass # electrons from the sky per aperture
        else:
            print "No sky level given. EXIT"
            sys.exit()
    
        RON_pix = 10. # electrons
        RON_aper = RON_pix * np.sqrt(aperture / DECam_pix**2) # electrons
    
        SNRout = source_electrons / np.sqrt(source_electrons + sky_electrons + RON_aper)
    
        return SNRout
    
    # find magnitude given target SNR
    def findmag(self, **kwargs):
        
        magsarray = np.linspace(15., 27., 100000)
        if kwargs['skymode'] == 'mag':
            SNRs = self.SNR(band=kwargs['band'], mag=magsarray, exptime=kwargs['exptime'], airmass=kwargs['airmass'], skymode=kwargs['skymode'], skymag=kwargs['skymag'])
        elif kwargs['skymode'] == 'mag-FWHM':
            SNRs = self.SNR(band=kwargs['band'], mag=magsarray, exptime=kwargs['exptime'], airmass=kwargs['airmass'], skymode=kwargs['skymode'], skymag=kwargs['skymag'], fwhm=kwargs['fwhm'])
        elif kwargs['skymode'] == 'ADU':
            SNRs = self.SNR(band=kwargs['band'], mag=magsarray, exptime=kwargs['exptime'], airmass=kwargs['airmass'], skymode=kwargs['skymode'], skyADU=kwargs['skyADU'], fwhm=kwargs['fwhm'])
        else:
            print "findmag: wrong keyword arguments"
            sys.exit()

        return magsarray[np.argmin((SNRs - kwargs['SNRin'])**2)]
    
if __name__ == "__main__":

    ETC = ETC_DECam(seeing_r_arcsec=0.75, vETC=6)
    print ETC.SNR(band='g', mag=20, exptime=173, airmass=1.0, skymode='ADU', skyADU=120, fwhm=1.0)
    print ETC.SNR(band='g', mag=20, exptime=173, airmass=1.0, skymode='mag', skymag=22)
    print ETC.SNR(band='g', mag=20, exptime=173, airmass=1.0, skymode='mag-FWHM', skymag=22, fwhm=1.0)
    print ETC.findmag(band='u', SNRin=5, exptime=173, airmass=1.6, skymode='mag', skymag=22.0) - ETC.findmag(band='u', SNRin=5, exptime=173, airmass=1.0, skymode='mag', skymag=22.0)

    ETC = ETC_DECam(seeing_r_arcsec=0.75, vETC=6)
    print ETC.FWHM(band='u', airmass=np.linspace(1.0, 1.6, 10))
    print ETC.FWHM(band='g', airmass=np.linspace(1.2, 1.9, 10))
    ETC = ETC_DECam(seeing_r_arcsec=1.0, vETC=6)
    print ETC.FWHM(band='u', airmass=np.linspace(1.0, 1.6, 10))
    print ETC.FWHM(band='g', airmass=np.linspace(1.2, 1.9, 10))
    #print ETC.findmag('g', 5, 160, 2.0, 'mag', 22.0) - ETC.findmag('g', 5, 160, 1.1, 'mag', 22.0)
    #print ETC.findmag('g', 5, 87, 1.8, 'mag', 22.0) - ETC.findmag('g', 5, 87, 1.1, 'mag', 22.0)
