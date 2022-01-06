# process GIMMS LAI to monthly values
import multiprocessing
import os, sys, time
import numpy as np
from netCDF4 import Dataset
import math
from skimage.measure import block_reduce

root = '/workdir/yzhang/Data/'
indir = root+'/BULAI_v1'

def getfile (directory):
    fileList=[]
    for path, subdirs, files in os.walk(directory):
        for name in files:
            if ".nc4" == name[-4:] or '.hdf' == name [-4:] or '.nc'==name[-3:]:
                fileList.append(os.path.join(path,name))
    #print fileList
    fileList.sort()
    return fileList

dlat = np.arange(-90, 90 + 1e-8, 0.5)
dlon = np.arange(-180, 180 + 1e-8, 0.5)

lat = np.arange(-90 + 0.5 / 2., (90 + 1e-8) - 0.5 / 2., 0.5)
lon = np.arange(-180 + 0.5 / 2., (180 + 1e-8) - 0.5 / 2., 0.5)

def export_nc(filename,lai_dat,time):
    n_obs = lai_dat.shape[0]
    f = Dataset(filename, 'w', format='NETCDF4')
    lati = f.createDimension('lat',len(lat))
    loni = f.createDimension('lon',len(lon))
    latitudes = f.createVariable('lat','f8','lat')
    longitudes = f.createVariable('lon','f8','lon')
    obs_num = f.createDimension('n_obs',n_obs)
    observation_number = f.createVariable('n_obs','f8','n_obs')

    latitudes.units = 'degrees_north'
    longitudes.units = 'degrees_east'
    latitudes.standard_name = 'latitude'
    longitudes.standard_name = 'longitude'
    latitudes.axis = 'Y'
    longitudes.axis = 'X'
    latitudes.long_name = 'latitude'
    longitudes.long_name = 'longitude'
    latitudes[:] = lat
    longitudes[:] = lon
    observation_number = time

    lai = f.createVariable('lai','f8',('n_obs','lat','lon'),zlib=T)
    #print(vardat)
    lai[:] = lai_dat
    f.close()


# get the files for the climate variables
files = getfile(indir)
print(files)
lai_dat = np.zeros((len(files),360,720))
for i in range(0,len(files)):
    with Dataset(files[i],"r") as fin:
        print(fin)
        dat = (fin.variables["lai"])
        #exec(cmd)
        print(dat.shape)
        # set bad values to NA
        temp_dat = np.where(dat<0,dat,np.nan)
        temp_reduce = block_reduce(temp_dat,block_size=(6,6),func = np.nanmean)
        # aggregate to 0.5
        lai_dat[i,:,:] = np.transpose(temp_reduce)

# save all data
outfile1 = '/workdir/yzhang/Data/GIMMS_3g_LAI/gimms3g.lai.15day.hd.nc'
time1 = np.repeat(np.arange(1981,2016)*100,24)+np.tile(np.arange(1,25),35)
time1 = time1[12:840]
export_nc(outfile1,lai_dat,time1)

outfile2 = '/workdir/yzhang/Data/GIMMS_3g_LAI/gimms3g.lai.mon.hd.nc'
time2 = np.repeat(np.arange(1981,2016)*100,12)+np.tile(np.arange(1,13),35)
time2 = time2[6:420]
lai_mon = block_reduce(lai_dat,block_size=(2,1,1),func = np.nanmean)
export_nc(outfile2,lai_mon,time2)



