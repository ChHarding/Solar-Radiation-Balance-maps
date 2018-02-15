# Solar-Radiation-Balance-maps

![TSR](https://github.com/ChHarding/Solar-Radiation-Balance-maps/blob/master/TSR/TSR_monthly_rotates_60.gif)

- High resolution Solar Radiation Balance maps from ECMWF 40-year reanalysis model data
- Folders contain maps (pdfs and gifs), model data and python code to create the maps
- You are free to use the maps or run the code to generate your own.
- Used in a lab activity to teach about climate in Geology/Astronomy 106L at Iowa State University: 
  - The learning objective for this lab is for students to understand how short- and long-wave radiation budgets affect monthly and yearly temperature fluctuations.
  - Students work in pairs, and are given a map that shows either annual mean incoming solar radiation or annual mean outgoing thermal radiation. 
  - Students answer a series of questions about these maps, and then combine with another group in order to compare the differences between incoming and outgoing radiation values. 
  - By doing this, they develop an understanding of which parts of the world experience an energy surplus or deficit.  
  - Students are then asked to hypothesize why the energy deficit (surplus) areas don’t continue to get colder (warmer)

## Folder content

data/ECMWF40_moda_Sep1957_Aug2002_SSR_STR_TSR_TTR.nc  
- netcdf file from http://apps.ecmwf.int/datasets/data/era40-moda
- Monthly means of Daily means, moda, 6, 2.5°, 1957-09-01 to 2002-08-01, Forecast, 40 years reanalysis
- Variables: Surface net solar radiation, Surface net thermal radiation, Top net solar radiation, Top net thermal radiation  (SSR, STR, TSR, TTR)
- Radiation Quantities in the ECMWF model and MARS.pdf documents the variables and gives some example maps

src/Radiation_maps_from_ECMWF40_data_monthly.py
src/Radiation_maps_from_ECMWF40_data_yearly.py
- Python 2 scripts for generating monthly or yearly maps from the Variables
- reads in the netcdf file and generates a bunch of geotiffs using GDAL so they could be used in a GIS
- uses numpy to store data in arrays
- data values (energy) is converted from J to W, duration is number of secs per month 
- for TTR and STR, the absolute values are used (otherwise they would be negative)
- calculates the yearly or monthly averages
- uses Matplotlib/Basemap to create the maps
- lots of plot parameters you could change for the maps:
  - variable (currently ttr but can be a list for batch processing)
  - area (currently global)
  - projection (currently Robinson)
  - visual center (median) of the map but can be a list for batch processing
  - grids (currently 30 degree spacings)
  - colormap/colorramp: currently CMRmap but can be a list for batch processing
  - contours (on/off)
  - for Radiation_maps_from_ECMWF40_data_yearly.py only one map is created, Radiation_maps_from_ECMWF40_data_yearly.py creates a map for each month
  - if batch processing is used (e.g. if a list of medians is given), a map is produced for each median, this can be combined with a list of colormaps 
  - monthly maps at different medians can be used as frames to create and animation showing the changes over the year, plus a slow rotation around the globe
  
geotiffs:
- rasters for monthly and annual means for ssr, str, ttr and tsr
- Radiation_data.mxd - ArcMap document for use in ArcGIS 10.3 or later
- shapefile for continents (continents.shp) and countries (world30.shp)
  
SSR: Surface net solar radiation maps (pdfs)

TSR: Top net solar radiation maps (pdfs and animated gif)

TTR: Top net thermal radiation (pdfs and animated gif)

( Surface net thermal radiation (STR) is in the netcdf file but no maps were created)

## Creating your own maps
- use ArcGIS with the rasters in the geotiffs folders (load Radiation_data.mxd to get started)

or:

- get python 2.7, install numpyt, netCDF4 osgeo (for gdal) matplotlib and Basemap
- decide if you want to plot yearly or monthly maps, edit the corresponding python script
- change this for the variable(s) to plot:

```
# variables to plot
#varnames = [ "ttr", "tsr", "str", "ssr"] 
varnames = [ "ttr"] 
```

- change this for the medians/meridians to plot:

```
    # plot at median      
    #medians = range(-90,181,90) # must be -180 to +180
    medians = [0]
```

- pick your colormap(s), here are some suggestions

```
	cmaps  = [plt.cm.CMRmap]
        """
                  plt.cm.gist_ncar,
                  plt.cm.gist_rainbow,
                  plt.cm.gist_earth,
                  plt.cm.gnuplot2,
                  plt.cm.CMRmap,
                  plt.cm.cubehelix,
                  plt.cm.BuPu_r,
                  plt.cm.YlGnBu_r,
        """
```
- fiddle with the basemap and matplotlob plot parameters
- run the code
- maps will end up in their variable folder (e.g. TTR),
- map filenames will contain the meridian (_mrd=XXX) and the colormap (_cmap=XXX) used.
- monthly maps will contain the month as number (_month=01 to 12)
- Radiation_maps_from_ECMWF40_data_monthly.py has some provisions for creating lower res jpgs which can be used to create animated gifs



