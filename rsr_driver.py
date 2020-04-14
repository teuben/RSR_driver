#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" A script to easily configure the Redshift Search Receiver (RSR) d
    data reduction
 (C) 2019 Large Millimeter Telescope (LMT), Mexico
 Released under GNU Public License (GPL v3)
 author: David O. Sánchez Argüelles (domars@inaoep.mx)
"""
try:
    import dreampy3 as dreampy
    from dreampy3.redshift.netcdf import RedshiftNetCDFFile
except ImportError:
    import dreampy
    from dreampy.redshift.netcdf import RedshiftNetCDFFile
from datetime import date
from scipy import signal
from collections import OrderedDict 
import sys
import numpy
import glob
import os.path
import argparse
import warnings

script_version ="0.2rc"

def rsrFileSearch (obsnum, chassis, root='/data_lmt/', full = True):
    """ Locate the RSR files for a given obsnum and chassis number
    
        Args:
          obsnum (int): Unique observation number for a LMT observation.
          chassis (int): Number of the RSR chassis. Must be in range [0-3].
          root (str): Absolute path where the LMT data is stored.
          full (bool): If true returns the absolute file path. If false return the relative file path. (Default is True)
        Returns:
          filename (str): Path to RSR data file if it exists or an empty string if the file is not found.
    """
    if chassis < 0 or chassis >3:
            print ("Error no chassis %d" % chassis)
            return ""

    chassisDir = "%s/RedshiftChassis%d/"%(root,chassis)
    filename = "*%06d*.nc" % (obsnum)
    
    fileList = glob.glob(chassisDir+filename)
    
    if len(fileList)>0:
        return fileList[0]

    return ""

def rsr_output_header(hdu, infodict):
    """ Returns a string with useful information for output spectrum
    
        Args:
           hdu (object): Header Data Unit processed from a RedshiftNetCDFFile
           infodict (dict): Dictionary containing additional information to verbose into the header. The elements of this dictionary must have been creater with the add_info function
        
        Returns:
           header(str): The header string with the data reduction details
        
        See:
            add_info
    """
    sigma = hdu.sigma.copy().squeeze()

    cdate = date.today()
    string = "------------Redshift Receiver Spectrum---------\n"
    string +="Telescope: Large Millimeter Telescope\n"
    string +="Source: %s\n" % hdu.header.SourceName
    string +="Source RA: %s\n" % hdu.header.RA
    string +="Source DEC: %s\n" % hdu.header.DEC
    string +="Pipeline version (DREAMPY): %s\n" % dreampy.version()
    string +="Driver script version: %s\n" % script_version 
    string +="Date of Reduction (YYYY-MM-DD): %s\n"%cdate.strftime("%Y-%b-%d")
    string +="Frequency Units: GHz\n"
    string +="Spectrum Units: K (T_A)\n"
    string +="Band intervals (GHz):"
    for i in range(sigma.shape[0]-1):
        string +="(%.3f-%.3f)," %(hdu.frequencies[i].min(), hdu.frequencies[i].max())
    string+= "(%.3f-%.3f)\n" %(hdu.frequencies[-1].min(), hdu.frequencies[-1].max())
    string +="Sigma per band: "
    for i in range(sigma.shape[0]-1):
        string +="%f," %sigma[i]
    string += "%f\n"%sigma[-1]
    for ikey in infodict.keys():
        string += "{}: {} {}\n".format(ikey, infodict[ikey]["value"], infodict[ikey]["units"])
    string +="------------------------------------------------\n"
    
    return string

def rebaseline (nc, fwind=33):
    """Performs a Savitzky-Golay filter on RSR band based data
    
       Args:
        nc (RedshiftNetCDFFile): Object containing the spectrum.
        fwind (int): Filter window. Number of samples to construct the filter on. Need to be an odd number. (Default 33).
       Note:
        Data in the nc.hdu.spectrum will be replaced with filtered data
    """
    
    if fwind %2 ==0:
        raise Exception ("You must provide an odd number for the Savitzky-Golay filter to work")
    nband = nc.hdu.spectrum.shape[1]
    data = nc.hdu.spectrum.copy().squeeze()
    freq = nc.hdu.frequencies
    
    for iband in range(nband):
        sec_start =[]
        sec_end =[]
        onScan=False
        for isample in range(data.shape[1]):
            if numpy.isfinite (data[iband, isample]) and data[iband, isample] != nc.hdu.blank and not onScan:
                onScan=True
                sec_start.append(isample)
                continue
            if onScan and ( not numpy.isfinite (data[iband, isample]) or data[iband, isample] == nc.hdu.blank):
                sec_end.append(isample)
                onScan = False
        nscans = len(sec_start)
        print ("Found %d scans on band %s on obs %s"%(nscans,iband, nc.hdu.header.obs_string()))
        for iscan in range(nscans):
            slen = sec_end[iscan]-sec_start[iscan]
            if fwind > slen-3:
                redfactor = 256/fwind
                awind = (slen)/redfactor
                if awind % 2 == 0:
                    awind+=1
                print("Requested SGF window %d is too large for valid data scan %d. Changing window to %d" %(fwind, slen, awind ) )
            else:
                awind = fwind
            data[iband,sec_start[iscan]:sec_end[iscan]] -= signal.savgol_filter(data[iband,sec_start[iscan]:sec_end[iscan]],awind,3)
    nc.hdu.spectrum[0,:,:] = data
    nc.hdu.baseline()
    
def freq_gaussian(freq,A,f_mu, f_width):
    """1D Gaussian Spectral Line 
    
        Args:
           freq (ndarray): Frequency axis
           A (float):   Line amplitude
           f_mu (float): Line central frequency
           f_width (float): Line width (Use FWHM).
        Return:
           gauss_line (ndarray): An array containing the spectral line
    """
    u = (freq-f_mu)
    sigma = 0.425*f_width
    u /= sigma
    return A*numpy.exp(-0.5*u**2)

def insert_sim(nc, A,f_mu, f_width):
    """ Insert a simulated gaussian line in the spectrum
    
        Args:
           nc (RedshiftNetCDFFile): NetCDF object with the processed  spectrum from a RedshiftNetCDFFile
           A (float):   Line amplitude
           f_mu (float): Line central frequency
           f_width (float): Line width (Use FWHM).
        Note:
            Updates the nc.hdu.spectrum variable
    """
    
    nreps = nc.hdu.spectrum.shape[0]
    nbands = nc.hdu.spectrum.shape [1]
    
    for iband in range(nbands):
        if f_mu > nc.hdu.frequencies[iband].min() and f_mu < nc.hdu.frequencies[iband].max():
            gmodel = freq_gaussian(nc.hdu.frequencies[iband], A, f_mu, f_width)
            #import ipdb; ipdb.set_trace()
            for irep in range(nreps):
                nc.hdu.spectrum[irep,iband]+=gmodel
    

def update_compspec_sigma(hdu):
    """ Updates the header sigma array with current standrard deviation values per band
    
        Args:
           hdu (object): Header Data Unit processed from a RedshiftNetCDFFile
        Note:
            Updates the hdu.sigma values
    """
    
    nbands = hdu.spectrum.shape[1]
    
    for iband in range(nbands):
        hdu.sigma[0,iband] = numpy.nanstd(hdu.spectrum[0,iband])
        
def add_info (infodict, label, value, units):
    """ Add information to a dictionary in order to be compatible with the rsr_output_header function
    
        Args:
            infodict (dict): Dictionary to be updated.
            label (str): Label (and key) to be added.
            value (any): Value to be written to the header.
            units (str): Unit of the value. Can be an empty string for unitless or adimentional units.
    """
    infodict[label]={"value":value, "units":units}
    

def vel2dfreq (freq, vel):
    """Convert a rotational velocity into the observed increment of frequency in a line f_width
    
        Args:
            freq (float):    Central Frequency, assumend in GHz
            vel (float):     Displacement velocity assumed in km/s
            
    """
    c =  299792.458
    beta = vel/c
    return freq*((1.+beta)/(1.-beta))**0.5-freq

def dfreq2vel(freq, dfreq):

    """Convert a increment of frequency in a velocity displacement
    
        Parameters:
            freq (float):   Central Frequency, assumend in GHz
            dfreq (float):   Increment in Frequency, assumed in GHz
            
    """
    c =  299792.458 
    
    nfreq = freq+dfreq
    return numpy.abs((freq**2-nfreq**2)/(nfreq**2+freq**2)*c)



def rsr_driver_start (clargs):
    
    """ Routine that calls the DREAMPY routines to reduce the data
    
        Args: 
            clargs (list): List of parameters from command line, use sys.argv[1:] or construct it using the shlex package
    """

    parser = argparse.ArgumentParser(description='Simple wrapper to process RSR spectra')

    parser.add_argument('obslist', help="Text file with obsnums to process, one obsnum per row")
    parser.add_argument('-p', dest ="doplot", action="store_true", help="Produce default plots")
    parser.add_argument('-t','--threshold',dest="cthresh", type=float, help="Thershold sigma value when coadding all observations")
    parser.add_argument('-o','--output', dest="output", default="", help="Output file name containing the spectrum")
    parser.add_argument('-f','--filter', dest="filter", default=0, help="Apply Savitzky-Golay filter (SGF) to reduce large scale trends in the spectrum. Must be an odd integer. This value represent the number of channels used to aproximate the baseline. Recomended values are larger than 21. Default is to not apply the SGF", type=int)
    parser.add_argument('-s','--smothing', dest ="smooth", default=0, type=int, help="Number of channels of a boxcar lowpass filter applied  to the coadded spectrum. Default is to no apply filter")
    parser.add_argument('-r','--repeat_thr', dest = "rthr", type = float,help="Thershold sigma value when averaging single observations repeats")  

    parser.add_argument('--simulate', nargs='+', help="Insert a simulated line into spectrum. The format is a list or a set of three elements Amplitude central_frequency line_velocity_width.", type = float)

    parser.add_argument('-d','--data_lmt_path', dest = "data_lmt",help="Path where the LMT data is located (defaults to /data_lmt", default="/data_lmt/")  


    args = parser.parse_args(clargs)

    hdulist = []
    windows = {}

    process_info=OrderedDict()

    tint = 0.0
    real_tint = 0.0
    source = None

    Obslist = numpy.loadtxt(args.obslist, unpack=True, dtype="int")

    if args.simulate:
        
        if len(args.simulate) != 3:
            raise Exception("Incorrect information provided to simulation function. Need 3 parameters separated by spaces")
        args.simulate[-1] = vel2dfreq(args.simulate[1], args.simulate[-1])

    for ObsNum in Obslist:
        for chassis in [0,1,2,3]:

            filename = rsrFileSearch(ObsNum, chassis, root=args.data_lmt)
            
            if filename == "":
                    print ('File not found for chassis %d ObsNumber: %s ' % (chassis,ObsNum))
                    continue
                    

            nc = RedshiftNetCDFFile(filename)
            if nc.hdu.header.ObsPgm != 'Bs':
                    nc.close()
                    continue
            if source is None:
                    source = nc.hdu.header.SourceName
        
            nc.hdu.process_scan(corr_linear=True) 
            
            if args.simulate:
                insert_sim(nc, *args.simulate)
            
            if chassis in(2,3): 
                    nc.hdu.blank_frequencies( {3: [(95.5,97.0),]} )
            nc.hdu.baseline(order = 1, subtract=True, windows=windows)
            
            if args.rthr:
                nc.hdu.average_all_repeats(weight='sigma', threshold_sigma=float(args.rthr))
            else:
                nc.hdu.average_all_repeats(weight='sigma')
            if args.filter >0:
                rebaseline(nc, args.filter)


            integ = 2*int(nc.hdu.header.get('Bs.NumRepeats'))*int(nc.hdu.header.get('Bs.TSamp'))
            tint += integ
            real_tint += (nc.hdu.data.AccSamples/48124.).mean(axis=1).sum()        
            hdulist.append(nc.hdu)
            nc.sync()
            nc.close()

    if len(hdulist) == 0:
        raise FileNotFoundError("No RSR data files were found. Check the -d parameter or the DATA_LMT environment variable")

    add_info(process_info, "Integration Time", real_tint, "s")

    if args.filter >0:
        freqd = args.filter*(nc.hdu.frequencies[1,0]-nc.hdu.frequencies[1,1])
        add_info(process_info, "SGF Window Size", args.filter," channels")
        add_info(process_info, "SGF Frequency HP Cut", freqd , "GHz")
        add_info(process_info, "SGF Velocity HP Cut", dfreq2vel(nc.hdu.frequencies[1,0],freqd) ,"km/s")

    if args.rthr:
        add_info(process_info, "Repeat Sigma Threshold", args.rthr, "K")

    if args.simulate:
        add_info(process_info, "Simulated: ", "{} K, {} GHz, {} GHz".format(args.simulate[0], args.simulate[1], args.simulate[2]), "") 

    hdu = hdulist[0]  # get the first observation
    if args.cthresh:
        hdu.average_scans(hdulist[1:],threshold_sigma=args.cthresh)
    else:
        hdu.average_scans(hdulist[1:])
    add_info(process_info,"Coadd Sigma threshold", args.cthresh, "K")

    if args.smooth > 0:
        add_info(process_info, "Smoothing Channels", args.smooth, "")
        hdu.smooth(nchan=args.smooth)
    hdu.make_composite_scan()
    
    update_compspec_sigma(hdu)

    if args.output == "":
        outfile = "%s_rsr_spectrum.txt"%source
    else:
        outfile = args.output
    out_array = numpy.array ([hdu.compfreq,hdu.compspectrum.squeeze()]).T
    numpy.savetxt(outfile, out_array, header = rsr_output_header(hdu,process_info))

    if args.doplot:
        try:
            from dreampy3.redshift.plots import RedshiftPlot
        except ImportError:
            from dreampy.redshift.plots import RedshiftPlot
        pl1 = RedshiftPlot()
        pl = RedshiftPlot()
        pl1.plot_spectra(hdu)
        pl.clear()
        pl.plot(hdu.compfreq, hdu.compspectrum[0,:], linestyle = 'steps-mid')
        pl.set_xlim(72.5, 111.5)
        pl.set_xlabel('Frequency (GHz)')
        pl.set_ylabel('TA* (K)')
        pl.set_subplot_title("%s Tint=%f hrs " %(hdu.header.SourceName, real_tint/3600.0))

if __name__ == "__main__":
    """ Simple wrapper to process RSR spectra
    
        See:
            rsr_driver -h
    
    """
    rsr_driver_start(sys.argv[1:])
