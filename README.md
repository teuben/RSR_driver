# RSR_driver
A wrapper script to work with the RSR DREAMPY pipeline

## Requirements
DREAMPY or ([DREAMPY3](https://github.com/lmt-heterodyne/dreampy3)) needs to be installed for the script to work. Also you need to configure a directory with the LMT data using the DATA_LMT environment variable. 

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

A deep integration on astronomical sources usually require several hours of observations. A typicall RSR spectrum is constructed by adding several five minute observations. Usually you observe yout favourite source by couple of hours each night. To avoid writting down all the observations of a single observing night you can use the *range notation* in the ObsNum file:

```
#Contents of file obsnum.txt
11222        #This adds a single obsnum
11223        #This add another obsnum
12111-12131  #This add all observations of the second night from 12111 to 12131 
13011-13015  #This add all observations of the third night from 12011 to 13015
```
Don't worry if there are calibration or pointing observations within the range, the wrapper will ignore those. Please be **carefull** with the range notation, you do not want to mix observations from different sources. Also it is **STRONGLY** advised to avoid writting ranges which span within multiple observation nights. Now it is time to execute the wrapper
```
$ python rsr_driver.py obsnum.txt
```
If you want to visualize the output just add the -p flag to the command and it will be produce a couple of interactive plots.

## Description of the output

By default the wrapper uses the obsnum input filename and adds the suffixes **_rsr_spectum** and **_rsr_spectrum_bandspec** to store the output spectra. 
To understand the difference between these files you need to know that the RSR 3mm frequency range is divided into 6 bands. The data of each band is stored and processed independently by most of the steps in the DREAMPY pipeline.

The **_rsr_spectrm_bandspec** contains the average spectrum of all the processed raw files separated per band. The shape of the data is [256,12]. The odd columns contain the frequency values in GHz and the even columns contain the spectrum values in antenna temperature units (K). In the DREAMPY plots this spectrum shows each band with a different color.

The **rsr_spectrum** is considered the main output of the pipeline. In this case the frequency channels that overlap within adjacent bands are averaged and a single frequency axis is defined. The shape of the data is [256,3]. The first column contains the frequency values in GHz, the second and the third contains the spectrum and the statistical error in antenna temperature units (K).

## Dealing with noisy data

DREAMPY uses a weighted average to produce the output spectra. In this scheme, the standard deviation of each band is used to compute the weights. However, sometimes this is not enough to produce good quality output spectrum. The rsr_driver exports two parameters to remove noisy data with the flags **-r** and **-t**. Both define a thershold levels for remove a particular band from the reduction process.

The **-r** flags works on a single raw file on a *repeat* level. The RSR is dual beam, dual polarization receiver. A standrard observation contains five repeats where 30 s are integrated with the source on one beam and 30 s with the source in the complementary polarization beam. This accounts for a  minute integration per repeat. As an example, to remove all the repeats which standard deviation exceeds 10mK you can issue the following command:

```
python rsr_driver.py obsnum.txt -r 0.01
```

The **-t** flag works when averaing all the raw files. This removes from the final average all the bands exceeding the value sent to the RSR_driver via this parameter. A good practice is to use a equal or smaller value than the **-r**. Otherwise it might have no effect.

## Notes on baseline subtraction

The DREAMPY package implements a polynomial baseline subtraction in the band data. You can specify the order of the polynomial by the **-b** flag. Notice that the maximum order implemented is 3. You can disable the baseline subtraction with the **--no-baseline-sub**.

Additionaly, you can use a high-pass Savitzky-Golay filter (SGF) to remove large variations in the baseline with the parameter **-f**. This parameter requires a number of samples used to calculate the sections of the filter. A larger value will consider more samples for a section of the filter hence it will remove larger scale variations. A recommended value is 55.


