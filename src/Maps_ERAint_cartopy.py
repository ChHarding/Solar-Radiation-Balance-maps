# Global Radiation maps for Cinzia Cervato and Diane Thatcher (Geol 106)
# Chris Harding (charding@iastate.edu) Sep 2015

# given a netCDF file of 1 or more variables, calculate the mean and
# plot it with basemap (matplotlib) and save it as a geotif file

# Ch: Apr. 2020: converted to Python 3 

import datetime
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt  # matplotlib state machine interface
import netCDF4
import numpy as np
from mpl_toolkits.basemap import Basemap, addcyclic
from osgeo import gdal, gdal_array, osr


# helper function: get date based on number of hours since Jan.1, 1900 (the 0 for this netcdf file)
def get_datetime(nhours): 
    w = nhours / (24 * 7) 
    delta = datetime.timedelta(weeks=w)  
    date = datetime.datetime(1900,  1, 1) + delta # see ncheader.txt
    return date 
#
# Save numpy array as geotiff using GDAL (http://www.gdal.org/gdal_tutorial.html)
#
def save_as_geotiff(data, filename, geotransform, projection):
    """ saves numpy 2D array as a unprojected geotiff
    data - 2D numpy array (must be float)
    filename - name of tiff file
    geotransform - GDAL Geotransform 
                   6 numbers, e.g. (0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
                      top left x, w-e pixel resolution, rotation - 0 if image is "north up", 
		      top left y, rotation - 0 if image is "north up", n-s pixel resolution
    projection - GDAL projection: "" for unprojected 
                 Note: this is IGNORED for now! GDAL supports projections such as NAD83 or UTM, but I've not tested 
		               anything but unprojected with WGS84
    """
    # check if we have a float array
    t = str(data.dtype)
    if not "float" in t: 
        print("save_as_geotiff(): data array must be float!")
        return

    driver  = gdal.GetDriverByName('GTiff')
    nrows,ncols = data.shape
    num_bands = 1
    output_raster = driver.Create(filename ,ncols, nrows, num_bands, gdal.GDT_Float32)  # Open the file for 32 bit float
    output_raster.SetGeoTransform(geotransform)  
    #output_raster.SetProjection(projection)      # set projection  - ignored for now
    
    # Hack: set to +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs, so ArcGIS knows that its a GCS
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS("WGS84") 
    output_raster.SetProjection(srs.ExportToWkt())  
                                
    # Write data into file
    output_raster.GetRasterBand(1).WriteArray(data)
    output_raster.GetRasterBand(1).ComputeStatistics(False) # False -> write full stats, not just the overview  stats
    output_raster = None # close buffer
    
    print("Created geotiff", filename)

# helper function to find names of variables in a .nc file 
# original at http://schubert.atmos.colostate.edu/~cslocum/netcdf_example.html#code
def ncdump(nc_fid):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object

    '''
    def print_ncattr(key):
        try:
            print("\t\ttype:", repr(nc_fid.variables[key].dtype))
            for ncattr in nc_fid.variables[key].ncattrs():
                s = nc_fid.variables[key].getncattr(ncattr)
                print('\t\t%s:' % ncattr, repr(s))
        except KeyError:
            print("\t\tWARNING: %s does not contain variable attributes" % key)

    print("NetCDF Global Attributes:")
    nc_attrs = nc_fid.ncattrs()
    for nc_attr in nc_attrs:
        s = nc_fid.getncattr(nc_attr)
        s = s[:19] if len(s) > 20 else s
        print('\t%s:' % nc_attr, repr(s))

    print("NetCDF dimension information:")
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    for dim in nc_dims:
        print("\tName:", dim)
        print("\t\tsize:", len(nc_fid.dimensions[dim]))
        print_ncattr(dim)

    print("NetCDF variable information:")
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    for var in nc_vars:
        if var not in nc_dims:
            print('\tName:', var)
            print("\t\tdimensions:", nc_fid.variables[var].dimensions)
            print("\t\tsize:", nc_fid.variables[var].size)
            print_ncattr(var)

# ---------------------- MAIN -----------------------------
start_folder = os.getcwd()
print("in folder", start_folder)


# netcdf file, dl from http://apps.ecmwf.int/datasets/data/era40-moda
# Monthly means of Daily means, moda, 6, 2.5, 1957-09-01...2002-08-01, Forecast, 40 years reanalysis
# This file contains:  Surface net solar radiation, Surface net thermal radiation, Top net solar radiation, 
#  Top net thermal radiation  (SSR, STR, TSR, TTR)

# I run ncdump -h on this file and put the header info in ncheader.txt, here's a part of it:
#dimensions:
	#longitude = 144 ;    <- number of cells along east-west
	#latitude = 73 ;      <- along NS
	#time = UNLIMITED ; // (540 currently) <- each cell has 540 timesteps, with a timestamp
#variables:
	#float longitude(longitude) ;
		#longitude:units = "degrees_east" ;
		#longitude:long_name = "longitude" ;
	#float latitude(latitude) ;
		#latitude:units = "degrees_north" ;
		#latitude:long_name = "latitude" ;
	#int time(time) ;
		#time:units = "hours since 1900-01-01 00:00:0.0" ; <- timestamp encodes how many hours! have pass since
		#time:long_name = "time" ;
		#time:calendar = "gregorian" 
	#short tsr(time, latitude, longitude) ;  <- the value value array,  540 timeslices, each 73 by 144
		#tsr:scale_factor = 151.303801138358 ;  <- tricky: take the cell value, first multiply by this
		#tsr:add_offset = 4957620.34809943 ;    <-         and then add this
		#tsr:_FillValue = -32767s ;
		#tsr:missing_value = -32767s ;         
		#tsr:units = "J m**-2" ;                <- to get the actual J/m2 value
		#tsr:long_name = "Top net solar radiation" ;
		#tsr:standard_name = "toa_net_upward_shortwave_flux" ;

# note for ERAint lpc and cp: values are already unpacked, so scale and offset are meaningless and
# misleading!


# name of netcdf file on disk 
print(os.listdir("data"))
#fn = "data/ECMWF40_moda_Sep1957_Aug2002_SSR_STR_TSR_TTR.nc" # I copied it here from the data folder 
#fn = r"C:\tmp\ERAint_79_19_LSP_CP_unpacked.nc"
fn = "data/ERAint_79_19_totprecip_SurfSolRadDwd.nc"
ncattr = {} # dict for storing attributes

print("reading netcdf file", fn)
nc = netCDF4.Dataset(fn, diskless=False) 

#ncdump(nc) # prints out metadata, in case you need to see the variable names

# variables to plot
#varnames = ["ttr", "tsr", "str", "ssr"] 
#varnames = ["ssdr"] # Surface solar radiation downwards
varname = "tp" # total precip 


def get_mean(varname, nc, save_geotiff=False, folder="geotiffs"):
    '''from opened nc, calculate and return the sum of varname over all month
       if save_as_geotiff is True, save mean as geotiff in folder folder'''

    # from the nc file, grab some attribute values we will need later (as a dict)
    # the time array, the lat/lon arrays and the value array (as numpy arrays)

    v  = nc.variables 
    curr_var = v[varname] # the nc variable
    
    print("working on", varname)
    for i in zip(curr_var.dimensions, curr_var.shape): # print out number of cells for each dimension
        print("\t",i[0], i[1])
        
    # get attributes of current var as a dict (see ncheader), 
    # see http://home.strw.leidenuniv.nl/~sfinx/netcdf4.html
    for a in ("scale_factor", "add_offset", "missing_value", "units", "long_name"):
       try:
           ncattr[a] = curr_var.getncattr(a)
           print("\t", a, ncattr[a])
       except:
           print(a, "not found")
    
    # get the time array (1D), each with the time in secs from 1/1/1900 for each of the 540 time slices
    time_ = np.array(v["time"]) 
    total_hours = time_[-1] - time_[0]
    print("\ttime: from", get_datetime(time_[0]), "to", get_datetime(time_[-1]), "total hours:",total_hours)
    print()

    # get data value array (3D)
    vl = np.array(v[varname])
    print(varname, "data values", vl.shape, vl.dtype, vl.min(), vl.max())  # we will later scale  and offset them  
    
    mean = vl.mean(axis=0) 
    print(" min %f  mean %f  max %f" % (vl.min(), vl.mean(), vl.max()))
    print("mean"); b,v = np.histogram(mean);print(b);print(v.round(2))

    # make histograms - this helped me to determine that the data was already unpacked
    '''
    import seaborn as sns
    sns.distplot(mean.flatten(), kde=False, bins=100,
                 hist_kws={"histtype": "step", "linewidth": 1,"alpha": 1, "color": "g"})

    import matplotlib.pyplot as plt
    plt.show()
    print()
    #std = vl.std(axis=0)
    #print("std"); b,v = np.histogram(std);print(b);print(v)
    '''
    
    # for my vars, no conversion is needed, so I'm commenting it out ...
    '''
    # For the sanity check we just approximated the number of seconds per month but to properly convert the data to W/m2
    # we really need to know the exact number of secs for each timeslice. 
    # _time contains all month as timestamps in hours(!)
    # We get the duration of each month in hours by subtracting the timestamp value of the next month 
    # from the timestamp value of the current month (works for all except the very last month ...),
    # and multiplying this difference by 3600 to get the duration in seconds.
    
    duration = np.zeros(len(time_)) # will contain number of secs per month (duraction of each time slice)
    i = 0
    while i < len(time_)-1: # don't include last! 
           num_hours = time_[i+1] - time_[i] # timestamp of next month - current -> number of hours in current month
           duration[i] = num_hours * 3600.0
           #print i,d
           i += 1    
    duration[i] = duration[i-12] # last month: use same as a year ago  
    #print(i,duration[i])  
       
    # convert each slice from J/m**2 to Watts/m**2 by div. it by its duration (in secs)
    sc  = ncattr["scale_factor"]
    ofs = ncattr["add_offset"]    
    for i,curr_slice in enumerate(vl):  # go through all slices, for loop will iterate by time dimension (axis=0 [time, lat, lon])
        #print(curr_slice.shape) # check that we have a 2D slice 
        curr_slice = curr_slice * sc + ofs #... data are first scaled before the offset is added.  http://home.strw.leidenuniv.nl/~sfinx/netcdf4.html:
        curr_slice /= duration[i] # convert J to W, duration is number of secs per month 
        #print(i, duration[i], np.mean(curr_slice))
        vl[i] = curr_slice # overwrite with Watt/s array
   
    # For ttr, negate values (as per Diana's request so the 106 students don't have to deal with negative energy)
    if varname == "ttr" or varname == "str": vl = -vl
    
    # get the mean along then time axis
    print("global min %f  mean %f  max %f" % (vl.min(), vl.mean(), vl.max())) # mean of converted 3D array
  
    sc  = ncattr["scale_factor"]
    ofs = ncattr["add_offset"]

    proc_vl = vl * sc + ofs
    print(proc_vl.min(), proc_vl.max())
    
    # b/c of the scale+offset issue we can't simply take the mean with mean()
    # but have to convert each slice first
    for i,curr_slice in enumerate(vl):
        print(curr_slice.min(), curr_slice.mean(), curr_slice.max())
        curr_slice = (curr_slice * sc + ofs) 
        print(curr_slice.min(), curr_slice.mean(), curr_slice.max())

    mean = vl.mean(axis=0) 
    print(" min %f  mean %f  max %f" % (vl.min(), vl.mean(), vl.max()))
    #print("mean"); b,v = np.histogram(mean);print(b);print(v)
    #std = vl.std(axis=0)
    #print("std"); b,v = np.histogram(std);print(b);print(v)
    
    '''

    #
    # Export as geotiff using GDAL
    #
    if save_geotiff == True:
        # make geotransform (http://gis.stackexchange.com/questions/8392/how-do-i-convert-affine-coordinates-to-lat-lng)
        # get cell size in degrees (should be 2.5, according to ECMWF40 docs, but lat has 73 not 72, so it must be -90 ... 0 ... +90)
        xres = 360/float(vl.shape[2]) # lon    
        yres = (180)/float(vl.shape[1]-1) # lat   -1 as we have 73 not  72
        xmin = 0; ymin = 90
        
        # I assume the data is cell centered I so offset by 1/2 a cell
        xmin -= xres/2.0 
        ymin += yres/2.0 
        #print "geotransform:", xmin, xres, ymin, yres   
        geotransform=(xmin,xres,0,ymin,0, -yres)     #    
        
        fn = varname + "_ann_mean.tif" # name of geotiff file
        fn = folder + os.sep + fn # put all in geotiffs folder
        save_as_geotiff(mean, fn, geotransform, None)

    return mean

#
# MAIN
#

# get mean
mean = get_mean(varname, nc, save_geotiff=True, folder="geotiffs")
mean *= 30 * 12 * 1000
print("min %f  mean %f  max %f" % (mean.min(), mean.mean(), mean.max()))

long_name = "Total rainfall mm/year"


# get lat/lon arrays (1D)
lat = np.array(nc.variables["latitude"])
#print lat, lat.shape, len(lat)
lon = np.array(nc.variables["longitude"])     

#
# Plot mean via matplotlib/Basemap and save as pdf
#

""" # as we will use a world wide projection (robinson), we don't need to set a bounding box
# However, other projection may require it, so I left this in as an example
# extent of mp to be plotted, around a specific meridian
mymeridian = 60 # border of map
llcrnrlon = -(360-mymeridian) #lon.min()
llcrnrlat = -90 #lat.min()
ucrnrlon  = mymeridian #lon.max()
urcrnrlat = 90 #lat.max() 
    
# if the data is not global, get the bounding box from the coords
llcrnrlon = lon.min()
llcrnrlat = lat.min()
ucrnrlon  = lon.max()
urcrnrlat = lat.max() 
#print "bbox ll ur lon lat", llcrnrlon, llcrnrlat, ucrnrlon, urcrnrlat
"""
    

# add one column, otherwise we get a gap on leftmost lon with basemap
mean,lon = addcyclic(mean, lon) 
x, y = np.meshgrid(lon, lat)  # makes 2D arrays from 1D lat and 1D lon

# plot at median      
#medians = range(-90,181,90) # must be -180 to +180
medians = [-90]
for med in medians:        
    
    # colormap(s) used for pcolormesh 
    # https://jakevdp.github.io/blog/2014/10/16/how-bad-is-your-colormap/
    #cmaps  = [plt.cm.nipy_spectral, plt.cm.CMRmap]
    cmaps  = [plt.cm.viridis_r]
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
    for cmap in cmaps:
        print("colormap", cmap.name)        
    
        # clip colormap  https://matplotlib.org/3.2.1/tutorials/colors/colormap-manipulation.html
        from matplotlib.colors import ListedColormap, LinearSegmentedColormap
        viridisBig = plt.cm.get_cmap('viridis_r', 512)
        max_color=0.85
        cmap = ListedColormap(viridisBig(np.linspace(0.0, max_color, 8)))

        # draw a basemap and keep handle 
        m = Basemap(
            llcrnrlon=-180,llcrnrlat=-60, urcrnrlon=180, urcrnrlat=70, 
            resolution='c', # c,l,i,h,f 
            #area_thresh = 1000,
            #projection='robin', # robinson projection
            #projection='moll',
            #projection='kav7',
            #projection='cyl',
            #projection='hammer',
            #lon_0=med,    # sets the median for some projections         
            ) 

        m.drawcoastlines(linewidth=1.5, 
                        #color='Black'
                        color='grey',
                        #color="#080808",
                        )
        m.fillcontinents(color='grey',lake_color='grey', alpha=0.3)
        m.drawcountries(linewidth=0.5, 
                        color='grey',
                        #color="#080808",
                        )
        #m.drawrivers(linewidth=1, color='Blue')
        #m.drawstates(linewidth=1.5) # US states
    
        m.drawparallels(np.arange(-90.,90.,30), labels=[1,1,0,0], linewidth=1) # labels left, right, top or bottom of the plot
        m.drawmeridians(np.arange(0.,360.,30), labels=[0,0,0,1], linewidth=1)
        
        # draw equator
        elo = [-180, 180]; ela = [0, 0]
        ex,ey = m(elo, ela) # project lat/lon to x/y         
        m.plot(ex, ey, 'k', linewidth=2)  


        # draw data (2-D array) as colored cells
        import matplotlib.colors as colors
        gamma = 0.35
        pcm = m.pcolormesh(lon,lat,
                            mean, 
                            latlon=True, 
                            cmap=cmap, 
                            #norm=colors.LogNorm(vmin=10, vmax=mean.max()),
                            norm=colors.PowerNorm(gamma=gamma),
                            alpha = 1.0,
                            shading="gouraud" 
                            )
        cb_ticks = np.arange(0, 10, 0.5) # range(-500,500,20) # tick marks on color bar
        intervals = (10, 25, 50, 100, 200, 300, 400, 500, 750, 1000, 1500, 2500, 5000, 7500, 10000)
        cbar = m.colorbar(pcm, location='bottom', 
                         pad="10%", 
                         #ticks=cb_ticks, 
                         ticks=intervals,
                         #extend='max'
                         format="%d",
                         )
        #cbar.set_label('W / m**2')  
        cbar.set_label('mm / year')   
        
        # contour plot
        #intv = 40 # stepsize for intervals
        #if varname == "ttr" or varname == "str": intv = 20
       
        
        #mpl.rcParams['contour.negative_linestyle'] = 'solid' # makes negtive contours solid, not dashed
        cs = m.contour(x,y,mean,
                    intervals,
                    latlon=True,
                    linewidths=0.7,
                    colors="black",  # dark grey
                    alpha=0.6,              
                    )
        
        plt.clabel(cs,  # the basemap (m) can't do contour labels (?) so I use the global plt object instead
                    inline=True, 
                    inline_spacing=-5, # sometime needed as gaps around numbers get too large
                    fontsize=8, 
                    fmt='%d') # format of numbers shows as labels 
        
        #cbar.add_lines(cs) # makes vertival lines at the contour interval numbers inside the color bar
        
        title_str = "Annual (1979 - 2019) mean of " + long_name + " max=" + str(max_color) +  " pow=" + str(gamma)
        plt.title(title_str, fontsize=18)    
        
        # put plot on a specific page size
        fig = mpl.pyplot.gcf() # get current figure
        fig.set_size_inches(17, 11)
        #fig.set_size_inches(11, 8.5)
        
        # make plot "bigger", less margins, pad (multiple of fontsize) has to be large or the lat labels get cut off!
        fig.tight_layout(pad=3.0, w_pad=0.4, h_pad=0.4)
                    
        # assemble file name for plt
        outfilename  = "total_precip" + "_mrd=" + str(med) +  "max=" + str(max_color) + ".pdf"
        plt.savefig(outfilename, dpi=75)
        
        plt.show() # show plot in matplotlib viewer
        print("wrote", outfilename)
        plt.clf() # clear canvas for next plot
    

print("done")




