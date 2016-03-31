# Public files for the HiTS survey

## ALTAZ animations:

>  altaz2014a.avi

>  altaz2015a.avi

## RADEC animations:

> 2014A_anim.gif

> 2015A_anim.gif

## Modified DECam ETC

> ETC_DECam.py

CC F. Forster, J. Martinez, J.C. Maureira, modified DECam ETC.
---------------------------------------------------

Please report any problems to francisco.forster@gmail.com
Options: help, test, seeing= [arcsec], ETC=, skymode=, band=, exptime= [sec], mag=, airmass=, fwhm= [arcsec], skymag= [mag/arcsec2], skyADU= [ADU/pix], SNR=

Code contains special class to call the ETC from inside python code, but can also be called from the command line.

Example inside python code:
```
   # import module
   from ETC_DECam import *

   # Initializing ETC_DECam object with seeing_r=of 0.75", version 6...
   ETC = ETC_DECam(seeing_r_arcsec=0.75, vETC=6)

   # Testing FWHM for an airmass vector...
   print "   u band, airmass between 1.0 and 1.6", ETC.FWHM(band='u', airmass=np.linspace(1.0, 1.6, 10))
   print "   g band, airmass between 1.2 and 1.9", ETC.FWHM(band='g', airmass=np.linspace(1.2, 1.9, 10))

   # Testing SNR and findmag with all skymodes with a 20 mag source, 173 sec exposure time in g band, airmass of 1.0...
   print ETC.SNR(band='g', mag=20, exptime=173, airmass=1.0, skymode='mag', skymag=22)
   print ETC.findmag(band='g', SNRin=SNRtest, exptime=173, airmass=1.0, skymode='mag', skymag=22.0)
```

Command line has two basic modes: giving an input magnitude to get an input signal to noise ratio (SNR) and viceversa, i.e. getting a limiting magnitude given a SNR.

Then, the skymode variabla controls how the environmental variables are defined:

mag: sky is given in mag/arcsec2 and the FWHM is derived from the airmass and seeing at zenith in r band (default is 0.75" at zenith)

ADU: sky is given in ADU/pixel and the FWHM is derived from the airmass and seeing at zenith in r band (default is 0.75" at zenith)

mag-FWHM: sky is given in mag/arcsec2 and FWHM is manually input in arcsec (airmass is also needed to compute extra extinction and atmosphere emission)

ADU-FWHM: sky is given in ADU/pixeland FWHM is manually input in arcsec (airmass is also needed to compute extra extinction and atmosphere emission)

Command line examples:
```
   python ETC_DECam.py --skymode mag --mag 20 --band g --exptime 173 --airmass 1.0 --skymag 22 
   SNR(mag=20.000000): 248.947527
   
   python ETC_DECam.py --skymode ADU --mag 20 --band g --exptime 173 --airmass 1.0 --skyADU 200 
   SNR(mag=20.000000): 259.685555
   
   python ETC_DECam.py --skymode mag-FWHM --mag 20 --band g --exptime 173 --airmass 1.0 --skymag 22 --fwhm 2
   SNR(mag=20.000000): 175.916219
   
   python ETC_DECam.py --skymode ADU-FWHM --mag 20 --band g --exptime 173 --airmass 1.0 --skyADU 200 --fwhm 2
   SNR(mag=20.000000): 191.589959
   
   python ETC_DECam.py --skymode mag --SNR 5 --band g --exptime 173 --airmass 1.0 --skymag 22 
   Magnitude(SNR=5.000000): 24.803138
   
   python ETC_DECam.py --skymode ADU --SNR 5 --band g --exptime 173 --airmass 1.0 --skyADU 200 
   Magnitude(SNR=5.000000): 24.945099
   
   python ETC_DECam.py --skymode mag-FWHM --SNR 5 --band g --exptime 173 --airmass 1.0 --skymag 22 --fwhm 2
   Magnitude(SNR=5.000000): 24.073411
   
   python ETC_DECam.py --skymode ADU-FWHM --SNR 5 --band g --exptime 173 --airmass 1.0 --skyADU 200 --fwhm 2
   Magnitude(SNR=5.000000): 24.216212
```
