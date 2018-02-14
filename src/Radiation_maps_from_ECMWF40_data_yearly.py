# Global Radiation maps for Cinzia Cervato and Diane Thatcher (Geol 106)
# Chris Harding (charding@iastate.edu) Sep 2015

# given a netCDF file of 1 or more variables, calculate the mean and
# plot it with basemap (matplotlib) and save it as a geotif file

import netCDF4
import os, sys
import datetime
import numpy as np
from osgeo import gdal
from osgeo import gdal_array
from osgeo import osr
import matplotlib.pyplot as plt # matplotlib state machine interface
import matplotlib as mpl 
from mpl_toolkits.basemap import Basemap, addcyclic


# helper function: get date based on number of hours since Jan.1, 1900 (the 0 for this netcdf file)
def get_datetime(nhours): 
    delta = datetime.timedelta(hours=nhours)  
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
	print "save_as_geotiff(): data array must be float!"
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
    output_raster = None # close 
    
    print "Created geotiff", filename
    

# ---------------------- MAIN -----------------------------
start_folder = os.getcwd()
print "in folder", start_folder

# netcdf file, dl from http://apps.ecmwf.int/datasets/data/era40-moda
# Monthly means of Daily means, moda, 6, 2.5°, 1957-09-01...2002-08-01, Forecast, 40 years reanalysis
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
		
# name of netcdf file on disk
fn = "ECMWF40_moda_Sep1957_Aug2002_SSR_STR_TSR_TTR.nc"  
		
ncattr = {} # dict for storing attributes

# variables to plot
#varnames = [ "ttr", "tsr", "str", "ssr"] 
varnames = [ "ttr"] 

print "reading netcdf file", fn
nc = netCDF4.Dataset(fn,diskless=False) 

# from the nc file, grab some attribute values we will need later (as a dict)
# the time array, the lat/lon arrays and the value array (as numpy arrays)
for varname in varnames:
    v  = nc.variables 
    curr_var = v[varname] # the nc variable
    
    print "working on", varname
    for i in zip(curr_var.dimensions, curr_var.shape): # print out number of cells for each dimentsion
        print "\t",i[0], i[1]
        
    # get attributes of current var as a dict (see ncheader), see http://home.strw.leidenuniv.nl/~sfinx/netcdf4.html
    for a in ("scale_factor", "add_offset", "missing_value", "units", "long_name"):
        ncattr[a] = curr_var.getncattr(a)
        print "\t", a, ncattr[a]
    
    # get the time array (1D), each with the time in secs from 1/1/1900 for each of the 540 time slices
    time_ = np.array(v["time"]) 
    total_hours = time_[-1] - time_[0]
    print "\ttime: from", get_datetime(time_[0]), "to", get_datetime(time_[-1]), "total hours:",total_hours
    print

    # get lat/lon arrays (1D)
    lat = np.array(nc.variables["latitude"])
    #print lat, lat.shape, len(lat)
    lon = np.array(nc.variables["longitude"])    

    # get data value array (3D)
    vl = np.array(v[varname])
    print "data values", vl.shape, vl.dtype    # we will later scale  and offset them  
    
    """
    # Sanity check: sample a spread of 10 values (vals) via histogram, scale and offset them to
    # first get J/m2, then convert to W/m2 by dividing by the number of seconds per month.
    bins,vals = np.histogram(vl) # 10 bin histogram with bin values in vals
    sc  = ncattr["scale_factor"]
    ofs = ncattr["add_offset"]
    time_in_secs = 730.46 * 3600  # seconds per timestep using an average number of hours per month 
    print "\nSanity check with histrogram-samples values:"
    print "raw value, *scale+offset (J/m2), /time (W/m2)"
    for v in vals:
        J_val = v * sc + ofs                       # http://home.strw.leidenuniv.nl/~sfinx/netcdf4.html: If both scale_factor and add_offset attributes are present, the data are first scaled before the offset is added.
        print int(v), J_val, J_val / time_in_secs  # raw, convert to J, convert J to W
    print
    """
    # For the sanity check we just approximated the number of seconds per month but to properly convert the data to W/m2
    # we really need to know the exact number of secs for each timeslice. 
    # _time contains all month as timestamps in hours(!)
    # We get the duration of each month in hours by subtracting athe timestamp value of the next month 
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
    #print i,duration[i]  
       
    
    # convert each slice from J/m**2 to Watts/m**2 by div. it by its duration (in secs)
    sc  = ncattr["scale_factor"]
    ofs = ncattr["add_offset"]    
    for i,curr_slice in enumerate(vl):  # go through all slices, for loop will iterate by time dimension (axis=0 [time, lat, lon])
        #print curr_slice.shape # check that we have a 2D slice 
        curr_slice = curr_slice * sc + ofs #... data are first scaled before the offset is added.  http://home.strw.leidenuniv.nl/~sfinx/netcdf4.html:
        curr_slice /= duration[i] # convert J to W, duration is number of secs per month 
        #print i, duration[i], np.mean(curr_slice)
        vl[i] = curr_slice # overwrite with Watt/s array
   
    # For ttr, negate values (as per Diana's request so the 106 students don't have to deal with negative energy)
    if varname == "ttr" or varname == "str": vl = -vl
     
    # get the mean along then time axis
    print "global min %f  mean %f  max %f" % (vl.min(), vl.mean(), vl.max()) # mean of converted 3D array
    mean = vl.mean(axis=0)
    #print "mean"; b,v = np.histogram(mean);print b;print v
    #std = vl.std(axis=0)
    #print "std"; b,v = np.histogram(std);print b;print v

    #
    # Export as geotiff using GDAL
    #
    
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
    fn = "geotiffs" + os.sep + fn # put all in geotiffs folder
    save_as_geotiff(mean, fn, geotransform, None)
    #continue # in case we don't want pdfs to be made
    
    
    #
    # Plot mean via matplotlib/Basemap and save as pdf
    #

    """ # as we will use a world wide projection (robinson), we don't need to set a bounding box
	# However, other projection may require it, so I left this in as an example
    # extent of mp to be plotted, around a precific meridian
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
    medians = range(-90,181,90) # must be -180 to +180
    #medians = [0]
    for med in medians:        
        
        # colormap(s) used for pcolormesh 
        # https://jakevdp.github.io/blog/2014/10/16/how-bad-is-your-colormap/
        #cmaps  = [plt.cm.nipy_spectral, plt.cm.CMRmap]
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
        for cmap in cmaps:
            print "colormap", cmap.name        
        
	    # draw a basemap and keep handle 
            m = Basemap(
                #llcrnrlon, llcrnrlat, ucrnrlon, urcrnrlat, 
                resolution='c', # c,l,i,h,f 
	        #area_thresh = 1000,
                projection='robin', # robinson projection
                #projection='moll',
                #projection='kav7',
                #projection='cyl',
                #projection='hammer',
                lon_0=med,    # sets the median for some projections         
                ) 

            m.drawcoastlines(linewidth=1.5, color='Black')
	    m.drawcountries(linewidth=0.01, color='#0B0B0B')
            #m.drawrivers(linewidth=1, color='Blue')
	    #m.drawstates(linewidth=1.5) # US states
 	    
            m.drawparallels(np.arange(-90.,90.,30), labels=[1,1,0,0], linewidth=1) # labels left, right, top or bottom of the plot
            m.drawmeridians(np.arange(0.,360.,30), labels=[0,0,0,1], linewidth=1)
            
            # draw equator
            elo = [-180, 180]; ela = [0, 0]
            ex,ey = m(elo, ela) # project lat/lon to x/y         
            m.plot(ex, ey, 'k', linewidth=3)            
    

            # draw data (2-D array) as colored cells
            pcm = m.pcolormesh(x,y,
                               mean, 
                               latlon=True, 
                               cmap=cmap, 
                               alpha = 1.0,
                               shading="gouraud" 
                               )
            cb_ticks = range(-500,500,20) # tick marks on color bar
            cbar = m.colorbar(pcm,location='bottom',pad="10%", ticks=cb_ticks)
            cbar.set_label('W / m**2')    
            
            # contour plot
            intv = 40 # stepsize for intervals
            if varname == "ttr" or varname == "str": intv = 20
            intervals = range(-500,500,intv)
            
            #mpl.rcParams['contour.negative_linestyle'] = 'solid' # makes negtive contours solid, not dashed
            cs = m.contour(x,y,mean,
                       intervals,
                       latlon=True,
                       linewidth=0.75,
                       colors='#080808', 
                       alpha=0.7,              
                       )
            plt.clabel(cs,  # the basemap (m) can't do contour labels (?) so I use the global plt object instead
	               inline=True, 
                       #inline_spacing=-3, # sometime needed as gaps around numbers get too large
                       fontsize=9, 
	               fmt='%d') # format of numbers shows as labels 
            cbar.add_lines(cs) # makes vertival lines at the contour interval numbers inside the color bar
            
            title_str = "Annual (1957 - 2002) mean of " + ncattr["long_name"]
            plt.title(title_str, fontsize=18)    
            
            # put plot on a specific page size
            fig = mpl.pyplot.gcf() # get's current figure
            fig.set_size_inches(17, 11)
            #fig.set_size_inches(11, 8.5)
            
            # make plot "bigger", less margins, pad (multiple of fontsize) has to be large or the lat labels get cut off!
            fig.tight_layout(pad=3.0, w_pad=0.4, h_pad=0.4)
                        
            # assemble file name for plt
            outfilename  = varname + "_mrd=" + str(med) + "_cmap="+cmap.name  + ".pdf"
            plt.savefig(outfilename, dpi=600)
            
            #m.show() # show plot in matplotlib viewer
            print "wrote", outfilename
            plt.clf() # clear canvas for next plot
	    

print "done"
sys.exit()




""" Header info
netcdf ECMWF40_moda_Sep1957_Aug2002_SSR_STR {
dimensions:
	longitude = 144 ;
	latitude = 73 ;
	time = UNLIMITED ; // (540 currently)
variables:
	float longitude(longitude) ;
		longitude:units = "degrees_east" ;
		longitude:long_name = "longitude" ;
	float latitude(latitude) ;
		latitude:units = "degrees_north" ;
		latitude:long_name = "latitude" ;
	int time(time) ;
		time:units = "hours since 1900-01-01 00:00:0.0" ;
		time:long_name = "time" ;
		time:calendar = "gregorian" ;
	short ssr(time, latitude, longitude) ;
		ssr:scale_factor = 125.433476263867 ;
		ssr:add_offset = 4109953.28326187 ;
		ssr:_FillValue = -32767s ;
		ssr:missing_value = -32767s ;
		ssr:units = "J m**-2" ;
		ssr:long_name = "Surface net solar radiation" ;
		ssr:standard_name = "surface_net_downward_shortwave_flux" ;
	short str(time, latitude, longitude) ;
		str:scale_factor = 67.7184777135184 ;
		str:add_offset = -1765199.35923886 ;
		str:_FillValue = -32767s ;
		str:missing_value = -32767s ;
		str:units = "J m**-2" ;
		str:long_name = "Surface net thermal radiation" ;
		str:standard_name = "surface_net_upward_longwave_flux" ;

// global attributes:
		:Conventions = "CF-1.0" ;
		:history = "2015-09-03 20:49:00 GMT by grib_to_netcdf-1.13.1: grib_to_netcdf /data/data04/scratch/netcdf-atls13-a562cefde8a29a7288fa0b8b7f9413f7-aVlFGs.target -o /data/data04/scratch/netcdf-atls13-a562cefde8a29a7288fa0b8b7f9413f7-L2Furv.nc -utime" ;
}



"""