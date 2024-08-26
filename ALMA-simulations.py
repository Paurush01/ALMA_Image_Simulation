import casatasks
from casatasks import importfits
from casatasks import simobserve
from casatasks import concat
from casatasks import simanalyze
from casatasks import simalma
from casatasks import exportfits
from casatasks import tclean
import cv2
from PIL import Image
from astropy.io import fits
import numpy as np
import os
import shutil


def preinducingnoise(imageloc,imgsize):
    #preinduces SNR=10 noise in the images which leads to various noise levels (SNR=3-15) in the final images
     image=(fits.getdata(imageloc, ext=0))
     snp = 10.0
     sky_sigma = np.max(image)/ snp
     image += np.random.normal(0,np.max(image)/snp,size=(imgsize, imgsize))
     header01=fits.getheader(imageloc)
     hdu=fits.PrimaryHDU(image,header=header01)
     hdu.writeto("noisy_"+imageloc,overwrite=True)

     

def ALMAsim(imageloc,imgsize,perpixelres,obsfreq,inwid,config,delfiles):
    #generates mock ALMA imaging for different configurations
    #imageloc (string) = file location of the frequency extracted fits file.
    #imgsize (integer) = number of pixels on any axis. for FLARES SKIRT output the value is 400
    #perpixelres (float) = per pixel resolution in arcseconds
    #obsfreq (float) = observed frame frequency of observation in Ghz
    #inwid (float) =  bandwidth of observation in Ghz
    #config (integer) = ALMA configuration from Cycle 10 (1-10). Change the hardcoded 10 to revert to previous cycles
    #delfiles (boolean) = Keep True to delete the CASA files created to simulate ALMA and only retain the final fits (just saves data)

    preinducingnoise(imageloc,imgsize)

    importfits(fitsimage="noisy_"+imageloc,imagename="ALMAsims.casa",whichrep=0,whichhdu=-1,zeroblanks=True,overwrite=True,defaultaxes=False,defaultaxesvalues=[],beam=[])
    freq=str(obsfreq)+"GHz"
    configname=""
    try:
        val = int(config)
        if(val<0) or (val>10):
            return("Invalid configuration number:")
        else :
            configname="alma.cycle10."+str(val)+".cfg"
    except ValueError:
        return ("Invalid configuration number")

    simobserve(project="ALMAsims",skymodel="ALMAsims.casa/",inbright="",indirection="J2000 19h00m00 -40d00m00",incell=str(perpixelres)+"arcsec",incenter=freq,inwidth=str(inwid)+"GHz",complist="",compwidth=str(inwid)+"GHz",setpointings=True,ptgfile="$project.ptg.txt",integration="20s",direction="",mapsize=['', ''],maptype="ALMA",pointingspacing="",caldirection="",calflux="1Jy",obsmode="int",refdate="2014/05/21",hourangle="transit",totaltime="1000s",antennalist=configname,sdantlist="aca.tp.cfg",sdant=0,thermalnoise="tsys-atm",user_pwv=5,t_ground=269.0,t_sky=260.0,tau0=0.1,seed=11111,leakage=0.0,graphics="both",verbose=False,overwrite=True)

    tclean(vis="./ALMAsims/ALMAsims."+"alma.cycle10."+str(val)+".ms", imagename="ALMAsims_tclean", phasecenter="J2000 19h00m00 -40d00m00",imsize=imgsize,cell=str(perpixelres)+"arcsec", specmode='mfs',deconvolver='hogbom', gridder='standard', weighting='natural', niter= 10000,minpsffraction=0.05, maxpsffraction=0.8)
    
    exportfits("ALMAsims_tclean.image",fitsimage="alma_c"+str(val)+"_"+imageloc,velocity=False,optical=False,bitpix=-32,minpix=0,maxpix=-1,overwrite=True,dropstokes=False,stokeslast=True,history=True,dropdeg=False)

    if (delfiles):
        shutil.rmtree("ALMAsims", ignore_errors=True)
        shutil.rmtree("ALMAsims.casa", ignore_errors=True)
        shutil.rmtree("ALMAsims.casa", ignore_errors=True)
        shutil.rmtree("ALMAsims_tclean.image", ignore_errors=True)
        shutil.rmtree("ALMAsims_tclean.mask", ignore_errors=True)
        shutil.rmtree("ALMAsims_tclean.model", ignore_errors=True)
        shutil.rmtree("ALMAsims_tclean.pb", ignore_errors=True)
        shutil.rmtree("ALMAsims_tclean.psf", ignore_errors=True)
        shutil.rmtree("ALMAsims_tclean.residual", ignore_errors=True)
        shutil.rmtree("ALMAsims_tclean.sumwt", ignore_errors=True)
        os.remove("noisy_"+imageloc)

    return ("FITS image saved as: "+ "alma_c"+str(val)+"_"+imageloc)


def FLARES_SKIRT_perpixelres(z):
    #just a shortcut for knowing pixel resolutions at different redshift snapshots for FLARES SKIRT files
    pixelrez=["0.023","0.025","0.028","0.030","0.032","0.035"]
    return (pixelrez[z-5])

def rest_to_obs(freq,z):
    #rest frame frequency should be entered in microns(um)
     lamda=freq*1e-6*(z+1)
     c=299792458.00
     f=(c/lamda)/1e9
     freq=round(f,2)
     return (freq)

if __name__ == '__main__':
    loc="z_6_flares_08_gal_054.fits"
    print(ALMAsim(loc,400,FLARES_SKIRT_perpixelres(6),rest_to_obs(158,6),7.5,3,True))
