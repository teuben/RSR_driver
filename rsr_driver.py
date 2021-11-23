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
    from dreampy3.redshift.netcdf.redshift_scan import RedshiftScan
    from dreampy3.redshift.utils.spectrum_utils import blank as spec_blank_value
except ImportError:
    import dreampy
    from dreampy.redshift.netcdf import RedshiftNetCDFFile
    from dreampy.redshift.netcdf.redshift_scan import RedshiftScan
    from dreampy.redshift.utils.spectrum_utils import blank as spec_blank_value
from datetime import date
from scipy import signal
from scipy import interpolate
from scipy import stats
from collections import OrderedDict 
import sys
import numpy
import glob
import os.path
import argparse
import warnings
import copy
import shlex

script_version ="0.6.0-pjt"


def rsrFileSearch (obsnum, chassis, root='/data_lmt/', full = True):
    """ Locate the RSR files for a given obsnum and chassis number
    
        Args:
          obsnum (int): Unique observation number for a LMT observation.
          chassis (int): Number of the RSR chassis. Must be in range [0-3].
          root (str): Absolute path where the LMT data is stored.
          full (bool): If true returns the absolute file path. If false return 
          the relative file path. (Default is True)
        Returns:
          filename (str): Path to RSR data file if it exists or an empty string 
          if the file is not found.
    """
    if chassis < 0 or chassis >3:
            print("Error no chassis %d" % chassis)
            return ""

    chassisDir = "%s/RedshiftChassis%d/"%(root,chassis)
    filename = "*%06d*.nc" % (obsnum)
    
    fileList = glob.glob(chassisDir+filename)
    
    if len(fileList)>0:
        return fileList[0]

    return ""

def rsr_output_header(hdu, infodict, add_comment = False):
    """ Returns a string with useful information for output spectrum
    
        Args:
           hdu (object): Header Data Unit processed from a RedshiftNetCDFFile
           infodict (dict): Dictionary containing additional information to verbose
           into the header. The elements of this dictionary must have been creater 
           with the add_info function
        
        Returns:
           header(str): The header string with the data reduction details
        
        See:
            add_info
    """
    sigma = hdu.sigma.copy().squeeze()

    cdate = date.today()
    
    ctag = ""
    if add_comment:
        ctag = "# "
    
    
    string = ctag + "------------Redshift Receiver Spectrum---------\n"
    string += ctag + "Telescope: Large Millimeter Telescope\n"
    string += ctag + "Source: %s\n" % hdu.header.SourceName
    string += ctag + "Source RA: %s\n" % hdu.header.RA
    string += ctag + "Source DEC: %s\n" % hdu.header.DEC
    string += ctag + "Pipeline version (DREAMPY): %s\n" % dreampy.version()
    string += ctag + "Driver script version: %s\n" % script_version 
    string += ctag + "Date of Reduction (YYYY-MM-DD): %s\n"%cdate.strftime("%Y-%b-%d")
    string += ctag + "Frequency Units: GHz\n"
    string += ctag + "Spectrum Units: K (T_A)\n"
    string += ctag + "Band intervals (GHz):"
    for i in range(sigma.shape[0]):
        wvalid = numpy.where(numpy.isfinite(hdu.spectrum[0,i]))
        if len(wvalid[0]) > 0:
            string +="(%.3f-%.3f)," %(hdu.frequencies[i,wvalid].min(), hdu.frequencies[i,wvalid].max())
    #string+= "(%.3f-%.3f)\n" %(hdu.frequencies[-1].min(), hdu.frequencies[-1].max())
    string = string[:-1]+"\n"
    string += ctag + "Sigma per band: "
    for i in range(sigma.shape[0]):
        string +="%f," %sigma[i]
    #string += "%f\n"%sigma[-1]
    string=string[:-1]+"\n"
    for ikey in infodict.keys():
        string += ctag + "{}: {} {}\n".format(ikey, infodict[ikey]["value"], infodict[ikey]["units"])
    string +=ctag +  "------------------------------------------------\n"
    
    return string

   

def rsr_savgolfilter(data_in, freq_in, fwind, exclude=None):
        
    """ Apply Savitzky-Golay filter (SGF) to the input data. This is a high-pass version of the filter.
    
        Args:
            freq_in (float ndarray): RSR Frequency
            data_in (float ndarray): RSR Spectrum
            fwind (float):           Window parameter to SGF
        
        Note: The RSR data (spectra) is mean subtracted before any calculation.
    """
    if  fwind %2 ==0:
        raise Exception ("You must provide an odd number for the Savitzky-Golay filter to work")
    
    slen = len (data_in)


    if fwind > slen-3:
        redfactor = 256//fwind
        awind = (slen)//redfactor
        if awind % 2 == 0:
            awind+=1
        print("Requested SGF window %d is too large for valid data scan %d. Changing window to %d" %(fwind, slen, awind ) )
    else:
        awind = fwind
    #Now we check if the user wants a section to be ignored from the SGF calculations

    data_filter = data_in.copy()
    freq_filter = freq_in.copy()
    spline_weight = numpy.ones(slen)

    if not exclude is None:
        for ifreq, iwidth in zip(exclude['freqs'], exclude['widths']):
            if ifreq >= freq_filter.min() and ifreq<=freq_filter.max():
                wex = numpy.where(numpy.logical_and(freq_filter >= ifreq-iwidth,
                                                    freq_filter <= ifreq+iwidth))
                wnex = numpy.where(numpy.logical_or(freq_filter < ifreq-iwidth,
                                                    freq_filter > ifreq+iwidth))
                spline_weight[wex] = 0.

                if freq_filter[-1] < freq_filter[0]:
                    dinterp = interpolate.UnivariateSpline(freq_filter[wnex][::-1], data_filter[wnex][::-1])#,s =len(data_filter)/4)
                else:
                    dinterp = interpolate.UnivariateSpline(freq_filter[wnex], data_filter[wnex])#,s =len(data_filter)/4)
                
                data_filter[wex] = dinterp(freq_filter[wex])
                
                print("SGF: Removing interval {} to {} GHz".format(ifreq-iwidth,ifreq+iwidth))
        
    data_filter = data_in - signal.savgol_filter(data_filter,awind,3)
    return data_filter

def rsr_notchfilter(freq_in, data_in, sigma):
            
    """ Implement a naive notch filter to RSR spectra in order to remove oscillations in the baseline. It work by analizing the RSR data
        PSD. If a peak in the PSD exceed the sigma thershold it is removed.
    
        Args:
            freq_in (float ndarray): RSR Frequency
            data_in (float ndarray): RSR Spectrum
            sigma (float):           Sigma threshold level to consider a baseline oscilation to be removed.
        
        Note: The RSR data (spectra) is mean subtracted before any calculation.
    """

    dspec = data_in.copy()
    fspec = freq_in.copy()
    
    dspec -=dspec.mean()
    
    fft_spec = numpy.fft.rfft(dspec)
    fft_freq = numpy.fft.rfftfreq(len(dspec), d=numpy.abs(fspec[1]-fspec[0]))
    spec_spec = numpy.abs(fft_spec)**2
    thspec,mins, maxs = stats.sigmaclip(spec_spec[1:], high= 4)
    w = numpy.where(spec_spec > sigma*thspec.std())
    fft_spec[0]=0.0
    fft_spec[w]=0.0
    dspec=numpy.fft.irfft(fft_spec)
    
    return dspec

def rebaseline (nc, farg=None, exclude = None, notch = False):
    """Performs a Savitzky-Golay (SG) or a Notch filter on RSR band based data. SG filter is suitable for low order trends in the data. Notch filter is sutable for    oscilations in the baseline. 
    
       Args:
        nc (RedshiftNetCDFFile): Object containing the spectrum.
        farg (int): Filter argument window. In case of the SG filter is the window lenght. For the notch filter is the sigma thershold for cutting a frequency
        exclude (dict): A dictionary with the frequencies to exclude from baseline calculations.
        notch (bool): If false (default) apply Savitzky-Golay filter. If True apply notch filter.
       Note:
        Data in the nc.hdu.spectrum will be replaced with filtered data
    """
    
    if farg is None:
        if not notch:
            farg = 55
        else:
            farg = 10
    
    if type(nc) is RedshiftNetCDFFile:
        hdu = nc.hdu
    elif type(nc) is RedshiftScan:
        hdu = nc
    else:
        raise TypeError("Argument passed to rebaseline must be a Redshift NetCDF File or Redshift Scan")
    nband = hdu.spectrum.shape[1]
    data = hdu.spectrum.copy().squeeze()
    freq = hdu.frequencies
    
    for iband in range(nband):
        sec_start =[]
        sec_end =[]
        onScan=False
        for isample in range(data.shape[1]):
            if numpy.isfinite (data[iband, isample]) and data[iband, isample] != hdu.blank and not onScan:
                onScan=True
                sec_start.append(isample)
                continue
            if onScan and ( not numpy.isfinite (data[iband, isample]) or data[iband, isample] == hdu.blank):
                sec_end.append(isample)
                onScan = False
        nscans = len(sec_start)
        print("Found %d scans on band %s on obs %s"%(nscans,iband, hdu.header.obs_string()))
        for iscan in range(nscans):
            
            if not notch:
                data_filter = rsr_savgolfilter(data[iband,sec_start[iscan]:sec_end[iscan]], \
                                               freq[iband,sec_start[iscan]:sec_end[iscan]], \
                                               int(farg),exclude=exclude)
            else:
                data_filter = rsr_notchfilter(freq[iband,sec_start[iscan]:sec_end[iscan]], \
                                              data[iband,sec_start[iscan]:sec_end[iscan]], \
                                              sigma=farg)
            
            data[iband,sec_start[iscan]:sec_end[iscan]] =data_filter
           
    hdu.spectrum[0,:,:] = data
    hdu.baseline()
    
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
    sigma = f_width/(2*(2*numpy.log(2))**0.5)
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
            for irep in range(nreps):
                nc.hdu.spectrum[irep,iband]+=gmodel
    

def update_compspec_sigma(hdu):
    """ Updates the header sigma array with current standard deviation values per band
    
        Args:
           hdu (object): Header Data Unit processed from a RedshiftNetCDFFile
        Note:
            Updates the hdu.sigma values
    """
    
    nbands = hdu.spectrum.shape[1]
    
    spectrum = hdu.spectrum.copy()
    spectrum[hdu.spectrum == spec_blank_value] = numpy.nan
    if numpy.ma.is_masked(spectrum):
        spectrum.mask[hdu.spectrum == spec_blank_value] = True
    
    hdu_sigma = hdu.sigma.copy()
    
    for iband in range(nbands):
        if numpy.ma.is_masked(spectrum):
            if spectrum[0,iband].mask.all():
                hdu_sigma[0,iband] = numpy.nan
            else:
                hdu_sigma[0,iband] = numpy.std(spectrum[0,iband])
        else:
            hdu_sigma[0,iband] = numpy.nanstd(spectrum[0,iband])
            
    #Now create a sigma array with the frequencies
    
    compsigma = numpy.zeros(hdu.compfreq.shape)
    bandlimits = setup_default_windows(hdu)
    
    
    nbands = bandlimits.keys()
    for i in range(hdu.compfreq.shape[0]):
        sigma_vals = []
        for ikey in nbands:
            if hdu.compfreq[i]>= bandlimits[ikey][0][0] and hdu.compfreq[i]<= bandlimits[ikey][0][1]:
                mindis = numpy.abs(hdu.frequencies[ikey] - hdu.compfreq[i])
                w = numpy.where(mindis == mindis.min())
                if numpy.isfinite(spectrum[0,ikey,w]):
                    sigma_vals.append(hdu_sigma[0,ikey])
        if len(sigma_vals) > 0:
            sigma_vals = numpy.array(sigma_vals)
            compsigma[i] = numpy.sqrt(numpy.sum(sigma_vals**2)/sigma_vals.shape[0]) 
        else:
            compsigma[i] = numpy.nan
    
    return compsigma
        
def add_info (infodict, label, value, units=""):
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

def setup_default_windows (nc):
    """ Generate a the default window values for baseline estimation
    
        Parameters:
            nc (RedshiftNetCDFFile):    A redshift raw data file to get the frequency values
            
        Returns:
            windows (dict): A dictionary with each band maximum a minimum values
    """
    
    if type(nc) is RedshiftNetCDFFile:
        hdu = nc.hdu
    elif type(nc) is RedshiftScan:
        hdu = nc
    else:
        raise TypeError("Argument passed to setup_default_windows must be a Redshift NetCDF File or Redshift Scan")
    
    nbands = hdu.frequencies.shape[0]
    windows ={}
    
    for i in range(nbands):
        windows[i] = [(hdu.frequencies[i].min(),hdu.frequencies[i].max() )]
    
    print("Generating default baseline windows values to include all data", windows)
    
    return windows

def update_windows (windows, limits, exclude):
    """ Modify the windows dictionaty to exlucde the interval freqs-widths to freqs+widths from the baseline estimation
    
        Parameters:
            windows (dict): Dictionary with the current values of windows
            limits (dict): Limits of each RSR band, usually a defaul windows dict
            exclude (dict): A dictionary with the freqs and widths to exclude 
    """
    freqs = exclude['freqs']
    widths = exclude['widths']
    
    sorti = numpy.argsort (freqs)
    freqs = freqs[sorti]
    widths = widths[sorti]
    
    exclude['freqs'] = freqs
    exclude['widths'] = widths
    
    
    for i in range(len(freqs)):
        for iband in range (len(limits)):
            if freqs[i]>=limits[iband][0][0] and freqs[i]<=limits[iband][0][-1]:
                newtuples = []
                for ituple in windows[iband]:
                    if freqs[i]>= ituple[0] and freqs[i] <= ituple[-1]:
                        newtuples.append((ituple[0],freqs[i]-widths[i]))
                        newtuples.append((freqs[i]+widths[i], ituple[1]))
                    else:
                        newtuples.append(ituple)
                windows[iband]= newtuples
                

def load_obsnum_file(ifile):
    """ Get a list ob observation number from the input file
    
        The raw LMT observation files are identified by a unique number (Observation Number or ObsNum). 
        The input file could contain a single observation number per line or two observation numbers
        separated by an hyphen, this notation is used to add to the observation list a range of files 
    
        Args: 
            ifile (sting): The input file containing the ObsNum
        
        Returns:
            A list with the ObsNums parsed from files
    """
    
    obsnum_list = []
    nline = 0
    with open (ifile) as dfile:
        for iline in dfile.readlines():
            nline += 1
            sline = iline.strip().split("-")

            if sline[0][0] == "#":
                continue
            if len(sline)==1:
                onum = int(sline[0])
                obsnum_list.append(onum)
            elif len(sline)==2:
                snum = int(sline[0])
                enum = int(sline[-1])
                obsnum_list.extend(range(snum,enum+1))
            else:
                print("Input file on line %d does not contain a valid ObsNum or ObsNume range. Ignoring")

    if len(obsnum_list) == 0:
        raise ValueError ("Input file does not contain any valid lines")
    return obsnum_list


def waterfall_plot (hdu, fig=None, thresh=None):
    
    from matplotlib import pyplot as plt
    
    nreps = hdu.spectrum.shape[0]
    
    if fig is None and nreps == 1:
        raise RuntimeError ("Averaged spectrum need a previously created figure")
    
    plot_data = hdu.spectrum.copy()
    w = numpy.where(plot_data == spec_blank_value)
    plot_data[w] = numpy.nan
    
    clip_data,_,_ = stats.sigmaclip(plot_data[numpy.isfinite(plot_data)],3)
    data_std = numpy.std(clip_data)
    
    alphas = numpy.ones(plot_data.shape[:2])
    if thresh:
        w = numpy.where(hdu.sigma >thresh)
        alphas[w] = 0.3
    
    onm = hdu.header.ObsNum
    try:
        chassis = hdu.header.ChassisNumber[0]
    except:
        chassis = hdu.header.ChassisNumber
    if fig is None:
        fig, axes = plt.subplots (nrows=nreps+1, ncols=1, sharex = True)
        for irep in range(nreps):
             for iband in range(6):
                axes[irep].plot(hdu.frequencies[iband], plot_data[irep,iband], alpha = alphas[irep,iband])
                axes[irep].set_ylim(-5*data_std, 5*data_std)
             axes[irep].text(hdu.frequencies.mean(),3*data_std,"ObsNum %d, Chassis %d, Repeat %d" %(onm,chassis,irep), ha="center")
    else:
        axes = fig.axes
        for iband in range(6):
            axes[-1].plot(hdu.frequencies[iband], plot_data[0,iband], alpha = alphas[0,iband])
            axes[-1].set_ylim(axes[0].get_ylim())
        axes[-1].text(hdu.frequencies.mean(),3.0/5.0*axes[0].get_ylim()[1],"ObsNum %d, Chassis %d, Averaged" %(onm,chassis), ha="center")
    
    return fig
    

def rsr_driver_start (clargs):
    
    """ Routine that calls the DREAMPY package to reduce the data
    
        Args: 
            clargs (list): List of parameters from command line, use sys.argv[1:] or construct it using the shlex package
    """

    parser = argparse.ArgumentParser(description='Simple wrapper to process RSR spectra')

    parser.add_argument('obslist', help="Text file with obsnums to process. Either one obsnum per row or a range of observation numbers separated by hyphens.")
    parser.add_argument('-p', dest ="doplot", action="store_true", help="Produce default plots")
    parser.add_argument('-t','--threshold',dest="cthresh", type=float, help="Thershold sigma value when coadding all observations")
    parser.add_argument('-o','--output', dest="output", default="", help="Output file name containing the spectrum")
    parser.add_argument('-f','--filter', dest="filter", default=0, help="Apply Savitzky-Golay filter (SGF) to reduce large scale trends in the spectrum. Must be an odd integer. This value represent the number of channels used to aproximate the baseline. Recomended values are larger than 21. Default is to not apply the SGF", type=int)
    parser.add_argument('-s','--smothing', dest ="smooth", default=0, type=int, help="Number of channels of a boxcar lowpass filter applied  to the coadded spectrum. Default is to no apply filter")
    parser.add_argument('-r','--repeat_thr', dest = "rthr", type = float,help="Thershold sigma value when averaging single observations repeats")  
    parser.add_argument('-n', '--notch_sigma', dest = "notch", help="Sigma cut for notch filter to eliminate large frecuency oscillations in spectrum. Needs to be run with -f option.", type =float)

    parser.add_argument('--simulate', nargs='+', help="Insert a simulated line into spectrum. The format is a list or a set of three elements Amplitude central_frequency line_velocity_width.", type = float)

    parser.add_argument('-d','--data_lmt_path', dest = "data_lmt",help="Path where the LMT data is located (default is to look for the DATA_LMT environment variable or the /data_lmt folder")  
    
    parser.add_argument('-b', dest="baseline_order", default = 1, help="Baseline calculation order", type=int)
    
    parser.add_argument('--exclude', nargs='+', dest ="exclude", help="A set of frequencies to exclude from baseline calculations. Format is central frequenciy width. Ej --exclude 76.0 0.2 96.0 0.3 excludes the 75.8-76.2 GHz and the 95.7-96.3 intervals from the baselien calculations.")
    
    parser.add_argument ('-j', dest="jacknife", action="store_true",help ="Perform jacknife simulation")


    parser.add_argument('-c', dest= "chassis", nargs='+', help = "List of chassis to use in reduction. Default is the four chassis")

    parser.add_argument('-B', '--badlags', help="A bad lags file with list of (chassis,board,channel) tuples as produced by seek_bad_channels")
    
    parser.add_argument('-R', '--rfile', help="A file with information of band data to ignore from analysis. \
                        The file must include the obsnum, chassis and band number to exclude separated by comas. One band per row")
    
    parser.add_argument('-w', '--waterfall-file', dest ="waterfall", help= "Request the driver to produce waterfall plot for each input file", type=str)

    parser.add_argument('--no-baseline-sub', dest="nosub", action="store_false", help="Disable subtraction of \
                        polinomial baseline. NOT RECOMMENDED.")

    args = parser.parse_args(clargs)

    hdulist = []
    windows = None
    exclude = None
    
    data_lmt_path = "/data_lmt/"
    res_data_lmt_path = ""
    if not args.data_lmt:
        if "DATA_LMT" in os.environ.keys():
            data_lmt_path = os.environ["DATA_LMT"]
    else:
        data_lmt_path=args.data_lmt
        if "DATA_LMT" in os.environ.keys():
            res_data_lmt_path = os.environ["DATA_LMT"]
    os.environ["DATA_LMT"]=data_lmt_path

    process_info=OrderedDict()

    tint = 0.0
    real_tint = 0.0
    source = None
    
    
    Obslist = load_obsnum_file(args.obslist)
    
    if numpy.ndim(Obslist)==0:
        tmplist = []
        tmplist.append(Obslist)
        Obslist = tmplist

    if args.simulate:
        if len(args.simulate) != 3:
            raise Exception("Incorrect information provided to simulation function. Need 3 parameters separated by spaces")
        args.simulate[-1] = vel2dfreq(args.simulate[1], args.simulate[-1])
    
    if args.baseline_order < 0 or args.baseline_order > 3:
        raise ValueError ("Requested a polynomial baseline subtraction of order %d. DREAMPY supports  \
                          order in range from 0 to 3. Check the input of -b parameter" % args.baseline_order)
    else:
        add_info(process_info, "Polynomial Baseline Order", args.baseline_order)
    
    if args.chassis:
        if len(args.chassis) > 4:
            raise ValueError ("Incorrect number of chassis supplied. RSR only have four chassis")
        
        use_chassis = numpy.array(args.chassis, dtype=int)
        use_chassis = tuple (numpy.unique(use_chassis))
    else:
        use_chassis = tuple  ((0,1,2,3))
        
    remove_keys = None
    rpattern = "O%06d_c%02d"
    if args.rfile:
        print("Processing rfile %s" % args.rfile)
        remove_keys ={}
        with open (args.rfile) as rfile:
            for iline in rfile.readlines():
                if iline[0] == '#':
                    continue
                if iline.isspace():
                    continue
                #  warning:   although called band, it's really board
                ronum, rchassis, rband = iline.split("#")[0].split(",")
                rkey = rpattern %(int (ronum),int(rchassis))
                if not rkey in remove_keys.keys():
                    remove_keys[rkey] = []
                remove_keys[rkey].append(int(rband))

    # PJT hack (either way, this will reset the lags if you don't supply them)
    if args.badlags:
        dreampy.badlags(args.badlags)
    else:
        dreampy.badlags(None)

    alltau = []
    waterpdf = None
    backend = None
    if args.waterfall:
        import matplotlib
        backend= matplotlib.get_backend()
        matplotlib.use("pdf")
        from matplotlib.backends.backend_pdf import PdfPages
        waterpdf = PdfPages(args.waterfall)
        

    info_obsnum=""

    for ObsNum in Obslist:
        nfiles_per_onum = 0
        tint_per_onum = 0.
        tau_onum = 0.
        ntau_onum = 0
        for chassis in use_chassis:

            filename = rsrFileSearch(ObsNum, chassis, root=data_lmt_path)
            
            if filename == "":
                    print('File not found for chassis %d ObsNumber: %s ' % (chassis,ObsNum))
                    continue
                    

            nc = RedshiftNetCDFFile(filename)
            if nc.hdu.header.ObsPgm != 'Bs':
                    nc.close()
                    continue
            if source is None:
                    source = nc.hdu.header.SourceName
                            
        
            try: 
            	nc.hdu.process_scan(corr_linear=True) 
            except:
                print("Cannot proces file %s. This is probably an incomplete observation" % filename)
                nc.close()
                continue
            
            info_obsnum += "%d_%d," % (ObsNum, chassis) 
            
            if windows is None:
                windows = setup_default_windows(nc)
                
                if args.exclude:
                    if len(args.exclude) % 2 !=0:
                        raise Exception ("Incorrect information provided in keyword exclude")
                    limits = copy.deepcopy(windows)
                    exclude = {}
                    exclude['freqs'] = numpy.array(args.exclude[0::2], dtype=float)
                    exclude['widths'] = numpy.array(args.exclude[1::2], dtype=float)
                    
                    update_windows(windows,limits,exclude)
                                   
            
            if args.simulate:
                insert_sim(nc, *args.simulate)
            
            if chassis in(2,3): 
                    nc.hdu.blank_frequencies( {3: [(95.5,97.0),]} )
            
            if not remove_keys is None:
                rk = rpattern % (ObsNum, chassis)
                if rk in remove_keys.keys():
                    for iband in remove_keys[rk]:
                        nc.hdu.blank_frequencies ({iband:[(windows[iband][0][0],windows[iband][0][-1]),]})
                        print("Remove board %d from ObsNum %d Chassis %d" % (iband, ObsNum, chassis))

            
            nc.hdu.baseline(order = args.baseline_order, subtract=True, windows=windows)
            
            if not waterpdf is None:
                wfig = waterfall_plot(nc.hdu,thresh=args.rthr)
            
            if args.rthr:
                nc.hdu.average_all_repeats(weight='sigma', threshold_sigma=args.rthr)
            else:
                nc.hdu.average_all_repeats(weight='sigma')
            if args.filter:
                rebaseline(nc, farg= args.filter,exclude=exclude)
            
            if not waterpdf is None:
                waterfall_plot(nc.hdu, wfig, thresh=args.cthresh)
                waterpdf.savefig(wfig)
                



            integ = 2*int(nc.hdu.header.get('Bs.NumRepeats'))*int(nc.hdu.header.get('Bs.TSamp'))
            itau = nc.hdu.header.get('Radiometer.Tau')
            if itau > 0.0 and numpy.isfinite(itau):
                tau_onum+=itau
                ntau_onum+=1
            tint += integ
            tint_per_onum += (nc.hdu.data.AccSamples/48124.).mean(axis=1).sum()
            nfiles_per_onum +=1
            hdulist.append(nc.hdu)
            nc.sync()
            nc.close()
        if nfiles_per_onum > 0:
            real_tint += tint_per_onum/nfiles_per_onum
        if ntau_onum >0:
            alltau.append(tau_onum/ntau_onum)
    
    avgtau = numpy.mean(alltau)

    if len(hdulist) == 0:
        raise FileNotFoundError("No RSR data files were found. Check the -d parameter or the DATA_LMT environment variable")
    
    add_info (process_info, "Input Observations", info_obsnum)

    if not waterpdf is None:
        waterpdf.close()

    add_info(process_info, "Integration Time", real_tint, "s")
    add_info(process_info, "Average Opacity (220GHz)", numpy.around(avgtau,2))
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
    
    if args.jacknife:
        rand_signs = numpy.sign(numpy.random.uniform(-0.5,0.5,size = len (hdulist)))
        for ihdu, isign in zip(hdulist,rand_signs):
            wvalues = numpy.where(numpy.logical_and(numpy.isfinite(ihdu.spectrum), ihdu.spectrum != spec_blank_value))
            if isign == 0:
                isign = 1
            ihdu.spectrum[wvalues]*=isign
    
    if args.cthresh:
        hdu.average_scans(hdulist[1:],threshold_sigma=args.cthresh)
        add_info(process_info,"Coadd Sigma threshold", args.cthresh, "K")
    else:
        hdu.average_scans(hdulist[1:])
    

    if args.filter:
        rebaseline(hdu,farg=args.filter, exclude=exclude)

    if args.notch:
        rebaseline(hdu,farg=args.notch,notch=True)
        add_info(process_info, "Notch Filter Sigma", args.notch, "")
    if args.smooth > 0:
        add_info(process_info, "Smoothing Channels", args.smooth, "")
        hdu.smooth(nchan=args.smooth)
        
    hdu.make_composite_scan()
    compsigma = update_compspec_sigma(hdu)
    
    add_info(process_info, "RSR driver cmd", " ".join(clargs))
    if args.output == "":
        outfile = "%s_rsr_spectrum.txt"%source
    else:
        outfile = args.output
    out_array = numpy.array ([hdu.compfreq,hdu.compspectrum.squeeze(), compsigma]).T
    numpy.savetxt(outfile, out_array, header = rsr_output_header(hdu,process_info))
    bspec_outfilename = outfile.replace(".txt", "_bandspec.txt")
    hdu.write_spectra_to_ascii(bspec_outfilename)
    #Prepend the header to the bandspec output file
    bspec_file = open(bspec_outfilename)
    bspec_lines = bspec_file.readlines()
    bspec_lines.insert(0,rsr_output_header(hdu,process_info, add_comment=True))
    bspec_file.close()
    bspec_file = open(bspec_outfilename,"w")
    bspec_file.writelines(bspec_lines)
    bspec_file.close()
    

    if args.doplot:
        if not backend is None:
            print("Attempt to switch backend to %s"%backend)
            try:
                matplotlib.use(backend, warn=False , force=True)
            except TypeError:
                matplotlib.use(backend, force=True)
            
        from matplotlib import pyplot as plt
        
        pl = plt.figure()
        plt.ion()
        plt.step(hdu.compfreq, hdu.compspectrum[0,:], where="mid")
        plt.xlim(72.5, 111.5)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('TA* (K)')
        plt.suptitle("%s Tint=%f hrs " %(hdu.header.SourceName, real_tint/3600.0))
        plb = plt.figure()
        for i in range(hdu.spectrum.shape[1]):
            plt.step(hdu.frequencies[i], hdu.spectrum[0,i,:], where="mid", label = "Band %d"%i)
        plt.legend(loc="best")
        plt.xlim(72.5, 111.5)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('TA* (K)')
        plt.suptitle("%s Tint=%f hrs " %(hdu.header.SourceName, real_tint/3600.0))
        plt.show()
        
        
if __name__ == "__main__":
    """ Simple wrapper to process RSR spectra
    
        See:
            rsr_driver -h
    """

    rsr_driver_start(sys.argv[1:]) 
