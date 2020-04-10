# RSR_driver
A collection of script to work with the RSR DREAMPY pipeline

## Requirements
DREAMPY needs to be installed for the script to work. Also you need to configure a directory with the LMT data.

## Usage

In a text file write down all the Observation Numbers corresponding to the raw RSR NetCDF files.

Then simply run: 
```
$ python rsr_driver.py <obsnum_file>
```

You can get information of all the parameter available by typying:
```
$ python rsr_driver.py -h
```
