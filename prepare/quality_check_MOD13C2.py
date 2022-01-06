import multiprocessing
import os, sys, time
import numpy as np
import gdal
import time
import datetime
from netCDF4 import Dataset
import fnmatch
from skimage.measure import block_reduce

global root, outdir

var_ID = int(sys.argv[1])-1

root = '/global/scratch/yaozhang/Data/MODIS_MOD13C2/'

outdir = '/global/scratch/yaozhang/Project/water_availability/data/'

def getfile(directory):
    fileList=[]
    for path, subdirs, files in os.walk(directory):
        for name in files:
            if ".nc4" == name[-4:] or '.hdf' == name [-4:]: 
                fileList.append(os.path.join(path,name))
    #print fileList
    fileList.sort()
    return fileList

def getyeardoy (file):
    yday = os.path.basename(file)[9:16]
    return yday

def quality_check(qa_vec):
    str_qa = "{0:b}".format(qa_vec).zfill(16)
    #print str_qa
    #print vi_small[i]
    if str_qa[1]=='0' and (str_qa[8:10]=='01' or str_qa[8:10]=='10') \
        and str_qa[5]=='0' and str_qa[0]=='0':
    #if str_qa[14:15]!='11':
        return True
    else:
        return False


def process_file(hdffile,varname):
    dict_name = {'NDVI': 'CMG 0.05 Deg Monthly NDVI',
            'EVI': 'CMG 0.05 Deg Monthly EVI',
            'QA': 'CMG 0.05 Deg Monthly VI Quality'}
    bands = {}
    ds_mod09 = gdal.Open(hdffile)
    subdata = ds_mod09.GetSubDatasets()
    subdataname = [item[0] for item in subdata]
    #print subdataname[1:18]
    xvar = [varname, "QA"]
    for i in xvar:
        dat_name = fnmatch.filter(subdataname,"*"+dict_name[i]+"*")[0]
        #print dat_name
        dataset = gdal.Open(dat_name)
        bands[i]=dataset.GetRasterBand(1).ReadAsArray()
    
    VI_qa = np.vectorize(quality_check)(bands["QA"])
    print(np.sum(VI_qa))
    qachecked_ndvi= np.where(VI_qa,bands[varname],np.nan)
    ndvi_hd = np.flipud(block_reduce(qachecked_ndvi/10000,block_size=(10,10),func=np.nanmean))
    valid_num_hd = np.flipud(block_reduce(np.isnan(qachecked_ndvi),block_size=(10,10),func=np.nansum))
    qa_ndvi_hd = np.where(valid_num_hd>90,np.nan,ndvi_hd)
    return qa_ndvi_hd
    
def export_nc(outfile,ndvi_mat,varname):
    ndvi_out = np.where(np.isnan(ndvi_mat),-999,ndvi_mat)
    dlat = np.arange(-90, 90 + 1e-8, 0.5)
    dlon = np.arange(-180, 180 + 1e-8, 0.5)
    lat = np.arange(-90 + 0.5 / 2., (90 + 1e-8) - 0.5 / 2., 0.5)
    lon = np.arange(-180 + 0.5 / 2., (180 + 1e-8) - 0.5 / 2., 0.5)
    tim = np.arange(1,ndvi_mat.shape[2]+1)
    #if file exist, remove first
    try:
        os.remove(outfile)
    except OSError:
        pass
    #create the new files

    f = Dataset(outfile, 'w', format='NETCDF4')
    lati = f.createDimension('lat', len(lat))
    loni = f.createDimension('lon', len(lon))
    time = f.createDimension("time",len(tim))
    latitudes = f.createVariable('lat', 'f8', 'lat')
    longitudes = f.createVariable('lon', 'f8', 'lon')
    timev = f.createVariable('time', 'i4', 'time')
    timev.units = 'months since 2000-01-01'
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
    timev[:] = tim
    
    ndvi=f.createVariable(varname, "f8", ("lat","lon","time"),fill_value = -999,zlib=True)
    ndvi[:,:,:] = ndvi_out
    f.close()

hdffiles = getfile(root)
hdffiles.sort()
# select file only before 2018
sel_files = hdffiles[:-3]
#process_file(hdffiles[2])


varnames = ["NDVI", "EVI"]
ndvi_mat = np.zeros((360,720,len(sel_files)))
ndvi_mat[:] = np.nan
for f in range(len(sel_files)):
    ndvi_mat[:,:,f] = process_file(sel_files[f],varnames[var_ID])

outfile1 = outdir+"MOD13C2_"+varnames[var_ID]+"_2000_2020_snow_cloud.nc"

export_nc(outfile1,ndvi_mat,varnames[var_ID])


