# retrieving monthly averaged ERA40 data from ecmwf
# you need to have a valid account with them to download
# for the parameters go to https://apps.ecmwf.int/datasets/data/interim-mdfa/levtype=sfc/
# and select them there (and a year), the go to View retrieval request - Python script
# to see the code for the variables you want. You'll also have to figure out
# the grid size, for ERA interim it should be 0.75/0.75

# make s string with all month for all years
#s = "19570901/19571001/19571101/19571201/"
s = ""
for y in range(1979, 2019): 
	for m in range(1,13):
		r = "%d%02d01/" % (y, m)
		s += r
#s += "20020101/20020201/20020301/20020401/20020501/20020601/20020701/20020801"
s = s[:-1]
print(s)
#exit()
#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": s, #"19790101/19790201/19790301/19790401/19790501/19790601/19790701/19790801/19790901/19791001/19791101/19791201",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "169.128/228.128", # LP and CP
    "step": "0-12/12-24",
    "stream": "mdfa",
    "type": "fc",
    "target": "ERAint_89_91_totprecip_SurfSolRadDwd.nc", # name of netcdf file
    #"resol": "AV",  # ??????????? 
    "format" : "netcdf",
})

'''
Alex:
Use this link here instead: https://apps.ecmwf.int/datasets/data/interim-mdfa/levtype=sfc/. I would download steps 0-12 and 12-24. You should be able to sum the amount of precipitation from both steps and that should give you the monthly mean precip in units of m/day. 
'''