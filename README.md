# RSR_driver
A wrapper script to work with the RSR DREAMPY pipeline

## Requirements
DREAMPY or DREAMPY3 needs to be installed for the script to work. Also you need to configure a directory with the LMT data using the DATA_LMT environment variable.

## Installation

To install start by cloning this repository:

```
$ git clone https://github.com/LMTdevs/RSR_driver.git
```

## Basic Usage

The Large Millimeter Telescope ([LMT](http://lmtgtm.org/)) raw observations are identified with an unique number known as **Observation Number** (ObsNum). To use the RSR_driver wrapper just write in a text file all the ObsNum you want to process, the easisest option one *ObsNum* per row.

Then simply run: 
```
$ python rsr_driver.py <obsnum_file>
```

You can get information of all the parameter available by typying:
```
$ python rsr_driver.py -h
```
## The ObsNum input file

Deep integration on astronomical sources usually require several hours of observations. A typicall RSR spectrum is constructed by adding several five minute observations. Usually you observe yout favourite source by couple of hour each night. To avoid writting down all the observations of a single observing night you can use the *range notation* in the ObsNum file:

```
#Contents of file obsnum.txt
11222        #This adds a single obsnum
11223        #This add another obsnum
12111-12131  #This add all observations of the second night
13011-13015  #This add all observations of the third night
```
Don't worry if there are calibration or pointing observation within the range, the wrapper will ignore those. Please be **carefull** with the range notation, you don't want to mix observations from different sources within the range. Also it is *STRONGLY* advised to avoid writting ranges which span within multiple observation nights. Now it is time to execute the wrapper
```
$ python rsr_driver.py obsnum.txt
```
If you want to visualize the output just add the -p flag to the command and it will be produce a couple of interactive plots.

## Description of the output

By default the wrapper uses the obsnum input filename and adds the suffixes **_rsr_spectum** and **_rsr_spectrum_bandspec** to store the output spectra. 
To understand the difference between files you need to know that the RSR 3mm band is divided into 6 bands. The data of each band is stored and processed independently by most of the steps of the DREAMPY pipeline.

The **_rsr_spectrm_bandspec** contains the average spectrum of al lthe processed raw file per band. The shape of the data is [256,12]. The odd columns contain the frequency values in GHz and the even columns contain the spectrum values in antenna temperature units (K).

The **rsr_spectrum** is considered the main output of the pipeline. In this case the frequency channels that overlap within adjacent band are averaged and a single frequency axis is defined. The shape of the data is [256,3]. The first column contains the frequency values in GHz, the second and the third contains the spectrum and the statistical error in antenna temperature units (K).
