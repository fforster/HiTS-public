import numpy as np
import sys, getopt

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
            self.vETC = 7 # recommended version

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
        if "nread" in kwargs.keys():
            nread = kwargs['nread']
        else:
            nread = 1
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
        if self.vETC == 7:
            aperture_eff = 2.04 
            RON_pix = 7. # readout noise per pixel, electrons
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
        elif self.vETC == 6:
            aperture_eff = 2.04 
            RON_pix = 10. # readout noise per pixel, electrons
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
            RON_pix = 10. # readout noise per pixel, electrons
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
        DECam_pix = 0.264 # arcsec  0.2637 in the centre, 0.2626 at the edges
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
    
        # readout noise per aperture
        if self.vETC == 5 or self.vETC == 6:
            RON_aper = RON_pix * np.sqrt(aperture / DECam_pix**2) # electrons^2/pixel^2
        elif self.vETC == 7:
            RON_aper = RON_pix**2 * (aperture / DECam_pix**2) # electrons^2/pixel^2
            
        # surce signal to noise ratio
        SNRout = source_electrons / np.sqrt(source_electrons + sky_electrons + nread * RON_aper)
    
        return SNRout
    
    # find magnitude given target SNR (very stupid search for now)
    def findmag(self, **kwargs):
        
        magsarray = np.linspace(15., 27., 100000)
        if kwargs['skymode'] == 'mag':
            SNRs = self.SNR(band=kwargs['band'], mag=magsarray, exptime=kwargs['exptime'], nread=kwargs['nread'], airmass=kwargs['airmass'], skymode=kwargs['skymode'], skymag=kwargs['skymag'])
        elif kwargs['skymode'] == 'mag-FWHM':
            SNRs = self.SNR(band=kwargs['band'], mag=magsarray, exptime=kwargs['exptime'], nread=kwargs['nread'], airmass=kwargs['airmass'], skymode=kwargs['skymode'], skymag=kwargs['skymag'], fwhm=kwargs['fwhm'])
        elif kwargs['skymode'] == 'ADU':
            SNRs = self.SNR(band=kwargs['band'], mag=magsarray, exptime=kwargs['exptime'], nread=kwargs['nread'], airmass=kwargs['airmass'], skymode=kwargs['skymode'], skyADU=kwargs['skyADU'])
        elif kwargs['skymode'] == 'ADU-FWHM':
            SNRs = self.SNR(band=kwargs['band'], mag=magsarray, exptime=kwargs['exptime'], nread=kwargs['nread'], airmass=kwargs['airmass'], skymode=kwargs['skymode'], skyADU=kwargs['skyADU'], fwhm=kwargs['fwhm'])
        else:
            raise ValueError("findmag: wrong keyword arguments")

        if 'SNRin' not in kwargs.keys():
            raise ValueError("findmag: missing input SNR")
        else:
            return magsarray[np.argmin((SNRs - kwargs['SNRin'])**2)]
    
if __name__ == "__main__":
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hts:e:o:b:t:n:m:a:f:M:A:S", ["help", "test", "seeing=", "ETC=", "skymode=", "band=", "exptime=", "nread=", "mag=", "airmass=", "fwhm=", "skymag=", "skyADU=", "SNR="])
    except getopt.GetoptError:
        print 'python ETC_DECam.py --help'
        sys.exit(1)
    for opt, arg in opts:

        if opt in ('-h', '--help'):
            print "\nCC F. Forster, J. Martinez, J.C. Maureira, modified DECam ETC.\n---------------------------------------------------\n"
            print "Please report any problems to francisco.forster@gmail.com"
            print "Options: help, test, seeing= [arcsec], ETC=, skymode=, band=, exptime= [sec], mag=, airmass=, fwhm= [arcsec], skymag= [mag/arcsec2], skyADU= [ADU/pix], SNR=\n"
            print "Code contains special class to call the ETC from inside python code, but can also be called from the command line."
            
            print "\nExample inside python code:\n"
            
            print """   # import module
   from ETC_DECam import *

   # Initializing ETC_DECam object with seeing_r=0.75\", version 7...
   ETC = ETC_DECam(seeing_r_arcsec=0.75, vETC=7)

   # Testing FWHM for an airmass vector...
   print "   u band, airmass between 1.0 and 1.6", ETC.FWHM(band='u', airmass=np.linspace(1.0, 1.6, 10))
   print "   g band, airmass between 1.2 and 1.9", ETC.FWHM(band='g', airmass=np.linspace(1.2, 1.9, 10))

   # Testing SNR and findmag with all skymodes with a 20 mag source, 173 sec exposure time in g band, airmass of 1.0...
   print ETC.SNR(band='g', mag=20, exptime=173, nread=1, airmass=1.0, skymode='mag', skymag=22)
   print ETC.findmag(band='g', SNRin=SNRtest, exptime=173, nread=1, airmass=1.0, skymode='mag', skymag=22.0)\n\n"""

            print "Command line has two basic modes: giving an input magnitude to get an input signal to noise ratio (SNR) and viceversa, i.e. getting a limiting magnitude given a SNR.\n"
            print "Then, the skymode variabla controls how the environmental variables are defined:"
            print "mag: sky is given in mag/arcsec2 and the FWHM is derived from the airmass and seeing at zenith in r band (default is 0.75\" at zenith)"
            print "ADU: sky is given in ADU/pixel and the FWHM is derived from the airmass and seeing at zenith in r band (default is 0.75\" at zenith)"
            print "mag-FWHM: sky is given in mag/arcsec2 and FWHM is manually input in arcsec (airmass is also needed to compute extra extinction and atmosphere emission)"
            print "ADU-FWHM: sky is given in ADU/pixeland FWHM is manually input in arcsec (airmass is also needed to compute extra extinction and atmosphere emission)\n"
            
            print "Command line examples:\n"
            print """   python ETC_DECam.py --skymode mag --mag 20 --band g --exptime 173 --nread 1 --airmass 1.0 --skymag 22 
   SNR(mag=20.000000): 248.947527
   python ETC_DECam.py --skymode ADU --mag 20 --band g --exptime 173 --nread 1 --airmass 1.0 --skyADU 200 
   SNR(mag=20.000000): 259.685555
   python ETC_DECam.py --skymode mag-FWHM --mag 20 --band g --exptime 173 --nread 1 --airmass 1.0 --skymag 22 --fwhm 2
   SNR(mag=20.000000): 175.916219
   python ETC_DECam.py --skymode ADU-FWHM --mag 20 --band g --exptime 173 --nread 1 --airmass 1.0 --skyADU 200 --fwhm 2
   SNR(mag=20.000000): 191.589959
   
   python ETC_DECam.py --skymode mag --SNR 5 --band g --exptime 173 --nread 1 --airmass 1.0 --skymag 22 
   Magnitude(SNR=5.000000): 24.803138
   python ETC_DECam.py --skymode ADU --SNR 5 --band g --exptime 173 --nread 1 --airmass 1.0 --skyADU 200 
   Magnitude(SNR=5.000000): 24.945099
   python ETC_DECam.py --skymode mag-FWHM --SNR 5 --band g --exptime 173 --nread 1 --airmass 1.0 --skymag 22 --fwhm 2
   Magnitude(SNR=5.000000): 24.073411
   python ETC_DECam.py --skymode ADU-FWHM --SNR 5 --band g --exptime 173 --nread 1 --airmass 1.0 --skyADU 200 --fwhm 2
   Magnitude(SNR=5.000000): 24.216212"""
            print "\n\n\n"
            

            sys.exit(1)

        elif opt in ('-t', '--test'):

            print "\nInitializing ETC_DECam object with seeing_r=of 0.75\", version 7..."
            
            ETC = ETC_DECam(seeing_r_arcsec=0.75, vETC=7)
            
            print "\nTesting FWHM for an airmass vector..."
            print "   u band, airmass between 1.0 and 1.6", ETC.FWHM(band='u', airmass=np.linspace(1.0, 1.6, 10))
            print "   g band, airmass between 1.2 and 1.9", ETC.FWHM(band='g', airmass=np.linspace(1.2, 1.9, 10))
            print "   r band, airmass between 1.2 and 1.9", ETC.FWHM(band='r', airmass=np.linspace(1.2, 1.9, 10))
            print "   i band, airmass between 1.2 and 1.9", ETC.FWHM(band='i', airmass=np.linspace(1.2, 1.9, 10))
            print "   z band, airmass between 1.2 and 1.9", ETC.FWHM(band='z', airmass=np.linspace(1.2, 1.9, 10))
            print "   Y band, airmass between 1.2 and 1.9", ETC.FWHM(band='Y', airmass=np.linspace(1.2, 1.9, 10))
            
            print "\nTesting SNR and findmag with all skymodes with a 20 mag source, 173 sec exposure time in g band, airmass of 1.0..."
            SNRtest = ETC.SNR(band='g', mag=20, exptime=173, nread=1, airmass=1.0, skymode='mag', skymag=22)
            magtest = ETC.findmag(band='g', SNRin=SNRtest, exptime=173, nread=1, airmass=1.0, skymode='mag', skymag=22.0)
            print "   mag (skymag=22): SNR %f <-> mag %f" % (SNRtest, magtest)
            if np.abs(20. - magtest) > 1e-4:
                print "   mag not OK"
            else:
                print "   OK"
            SNRtest = ETC.SNR(band='g', mag=20, exptime=173, nread=1, airmass=1.0, skymode='mag-FWHM', skymag=22, fwhm=1.0)
            magtest = ETC.findmag(band='g', SNRin=SNRtest, exptime=173, nread=1, airmass=1.0, skymode='mag-FWHM', skymag=22, fwhm=1.0)
            print "   mag-FWHM (skymag=22, fwhm=1): SNR %f <-> mag %f"  % (SNRtest, magtest)
            if np.abs(20. - magtest) > 1e-4:
                print "   mag-FWHM not OK"
            else:
                print "   OK"
            SNRtest = ETC.SNR(band='g', mag=20, exptime=173, nread=1, airmass=1.0, skymode='ADU', skyADU=120)
            magtest = ETC.findmag(band='g', SNRin=SNRtest, exptime=173, nread=1, airmass=1.0, skymode='ADU', skyADU=120)
            print "   ADU (skyADU=120): SNR %f <-> mag %f" % (SNRtest, magtest)
            if np.abs(20. - magtest) > 1e-4:
                print "   ADU not OK"
            else:
                print "   OK"
            SNRtest = ETC.SNR(band='g', mag=20, exptime=173, nread=1, airmass=1.0, skymode='ADU-FWHM', skyADU=120, fwhm=1.0)
            magtest = ETC.findmag(band='g', SNRin=SNRtest, exptime=173, nread=1, airmass=1.0, skymode='ADU-FWHM', skyADU=120, fwhm=1.0)
            print "   ADU-FHWM (skyADU=120, fwhm=1): SNR %f <-> mag %f" % (SNRtest, magtest)
            if np.abs(20. - magtest) > 1e-4:
                print "   ADU-FWHM not OK"
            else:
                print "   OK"

            print "\nTesting findmag for input SNR of 5 and exposure time of 173 sec, assuming sky magnitudes of 22.8, 22.1, 21.1, 20.1, 18.7 and 18 in ugrizY..."
            print "   u band, mags at airmasses of 1.0, 1.2, 1.4, 1.6, 1.8, 2.0", map(lambda x: ETC.findmag(band='u', SNRin=5, exptime=173, nread=1, airmass=x, skymode='mag', skymag=22.8), [1., 1.2, 1.4, 1.6, 1.8, 2.0])
            print "   g band, mags at airmasses of 1.0, 1.2, 1.4, 1.6, 1.8, 2.0", map(lambda x: ETC.findmag(band='g', SNRin=5, exptime=173, nread=1, airmass=x, skymode='mag', skymag=22.1), [1., 1.2, 1.4, 1.6, 1.8, 2.0])
            print "   r band, mags at airmasses of 1.0, 1.2, 1.4, 1.6, 1.8, 2.0", map(lambda x: ETC.findmag(band='r', SNRin=5, exptime=173, nread=1, airmass=x, skymode='mag', skymag=21.1), [1., 1.2, 1.4, 1.6, 1.8, 2.0])
            print "   i band, mags at airmasses of 1.0, 1.2, 1.4, 1.6, 1.8, 2.0", map(lambda x: ETC.findmag(band='i', SNRin=5, exptime=173, nread=1, airmass=x, skymode='mag', skymag=20.1), [1., 1.2, 1.4, 1.6, 1.8, 2.0])
            print "   z band, mags at airmasses of 1.0, 1.2, 1.4, 1.6, 1.8, 2.0", map(lambda x: ETC.findmag(band='z', SNRin=5, exptime=173, nread=1, airmass=x, skymode='mag', skymag=18.7), [1., 1.2, 1.4, 1.6, 1.8, 2.0])
            print "   Y band, mags at airmasses of 1.0, 1.2, 1.4, 1.6, 1.8, 2.0", map(lambda x: ETC.findmag(band='Y', SNRin=5, exptime=173, nread=1, airmass=x, skymode='mag', skymag=18), [1., 1.2, 1.4, 1.6, 1.8, 2.0])

            print "\n\n\n\n"
            sys.exit()
            

        elif opt in ('-s', '--seeing'):
            seeing_r_arcsec = float(arg)
        elif opt in ('-e', '--ETC'):
            vETC = int(arg)
        elif opt in ('-o', '--skymode'):
            skymode = arg
        elif opt in ('-b', '--band'):
            band = arg
        elif opt in ('-t', '--exptime'):
            exptime = float(arg)
        elif opt in ('-n', '--nread'):
            nread = float(arg)
        elif opt in ('-m', '--mag'):
            mag = float(arg)
        elif opt in ('a', '--airmass'):
            airmass = float(arg)
        elif opt in ('-f', '--fwhm'):
            fwhm = float(arg)
        elif opt in ('-M', '--skymag'):
            skymag = float(arg)
        elif opt in ('-A', '--skyADU'):
            skyADU = float(arg)
        elif opt in ('-S', '--SNR'):
            SNR = float(arg)

    if 'seeing_r_arcsec' in locals() and 'vETC' in locals():
        ETC = ETC_DECam(seeing_r_arcsec=seeing_r_arcsec, vETC=vETC)
    elif 'seeing_r_arcsec' in locals():
        ETC = ETC_DECam(seeing_r_arcsec=seeing_r_arcsec)
    elif 'vETC' in locals():
        ETC = ETC_DECam(vETC=vETC)
    else:
        ETC = ETC_DECam()

    # manage exceptions, missing lots of them!
    if 'skymode' in locals():
        
        if 'SNR' in locals():

            if skymode == 'mag':
                print "Magnitude(SNR=%f): %f" % (SNR, ETC.findmag(band=band, SNRin=SNR, exptime=exptime, nread=nread, airmass=airmass, skymode='mag', skymag=skymag))
            elif skymode == 'mag-FWHM':
                print "Magnitude(SNR=%f): %f" % (SNR, ETC.findmag(band=band, SNRin=SNR, exptime=exptime, nread=nread, airmass=airmass, skymode='mag-FWHM', skymag=skymag, fwhm=fwhm))
            elif skymode == 'ADU':
                print "Magnitude(SNR=%f): %f" % (SNR, ETC.findmag(band=band, SNRin=SNR, exptime=exptime, nread=nread, airmass=airmass, skymode='ADU', skyADU=skyADU))
            elif skymode == 'ADU-FWHM':
                print "Magnitude(SNR=%f): %f" % (SNR, ETC.findmag(band=band, SNRin=SNR, exptime=exptime, nread=nread, airmass=airmass, skymode='ADU-FWHM', skyADU=skyADU, fwhm=fwhm))
        
        else:
            
            if skymode == 'mag':
                print "SNR(mag=%f): %f" % (mag, ETC.SNR(band=band, mag=mag, exptime=exptime, nread=nread, airmass=airmass, skymode='mag', skymag=skymag))
            elif skymode == 'mag-FWHM':
                print "SNR(mag=%f): %f" % (mag, ETC.SNR(band=band, mag=mag, exptime=exptime, nread=nread, airmass=airmass, skymode='mag-FWHM', skymag=skymag, fwhm=fwhm))
            elif skymode == 'ADU':
                print "SNR(mag=%f): %f" % (mag, ETC.SNR(band=band, mag=mag, exptime=exptime, nread=nread, airmass=airmass, skymode='ADU', skyADU=skyADU))
            elif skymode == 'ADU-FWHM':
                print "SNR(mag=%f): %f" % (mag, ETC.SNR(band=band, mag=mag, exptime=exptime, nread=nread, airmass=airmass, skymode='ADU-FWHM', skyADU=skyADU, fwhm=fwhm))

    else:
        
        print "Define skymode (mag, mag-FWHM, ADU, ADU-FWHM)"
        sys.exit()
        
    

