### dlm_ndvi_analysis
import numpy as np
from numpy import matlib
from netCDF4 import Dataset
from datetime import datetime
from dlm_functions import forwardFilteringM, Model,perpixel_vi,de_seasonalize,deseasonalize_dlm,perpixel_vi_multi
import pandas as pd
import sys
import random
import glob,os
from multiprocessing import Pool
#import dask.dataframe as dd
#from dask.multiprocessing import get
#%%
spin_up_years = int(sys.argv[1])
spin_up_rep = int(sys.argv[2])
delta_v = int(sys.argv[3])


deltas = [0.95,0.98,0.99,0.995,0.999]
delta = deltas[delta_v]


print("spin_up_years: ",spin_up_years)
print("spin_up_repeat: ",spin_up_rep)
print("delta: ",delta)
#print("dataset: ",ind)
#print("random_seed: ",rand_seed)
rseas = [1,2,3]

ind_random_spei=True
ind_process_dlm=True

st = "_2001"       #""
yr_ind = 1176      #942
#%%
def parallelize_dataframe_multi(df, n_cores=4):
    df_split = np.array_split(df, n_cores,axis=1)
    pool = Pool(n_cores)
    df = np.concatenate(pool.map(unpacking_apply_along_axis_vi_multi, df_split),axis=1)
    pool.close()
    pool.join()
    return df


def unpacking_apply_along_axis_vi_multi(arr):
    #print(lag_mon,spin_up_years)
    return np.apply_along_axis(perpixel_vi_multi,0, arr,spin_up_years=spin_up_years,spin_up_rep=spin_up_rep,delta=delta,rseas=rseas)


def exportnc(results, outfile,n_regvar):
    n_var=2+n_regvar+len(rseas)*2
    out_dim = int(results.shape[0]/(n_var*2+2+n_regvar+1))

    slik = results[0:out_dim,:]
    sm = results[out_dim:((n_var+1)*out_dim),:]
    sC = results[((n_var+1)*out_dim):((n_var*2+1)*out_dim),:]
    snu = results[((n_var*2+1)*out_dim):((n_var*2+2)*out_dim),:]
    Xpre = results[((n_var*2+2)*out_dim):((n_var*2+2+n_regvar)*out_dim),:]    
    Yout = results[((n_var*2+2+n_regvar)*out_dim):((n_var*2+3+n_regvar)*out_dim),:]
 
    slik = slik.reshape((out_dim,360,720))
    sm = sm.reshape((n_var,out_dim,360,720))
    sC = sC.reshape((n_var,out_dim,360,720))
    snu = snu.reshape((out_dim,360,720))
    Xpre = Xpre.reshape((n_regvar,out_dim,360,720))
    Yout = Yout.reshape((out_dim,360,720))

    f = Dataset(outfile,'w',format="NETCDF4")
    lon = np.arange(-180+0.25,180-0.25+0.01,0.5)
    lat = np.arange(-90+0.25,90-0.25+0.01,0.5)
    ti = np.arange(0,out_dim)
    va = np.arange(0,n_var)
    rvar = np.arange(0,n_regvar)

    f.createDimension('lon',len(lon))
    f.createDimension('lat',len(lat))
    f.createDimension('tim',len(ti))
    f.createDimension('var',n_var)
    f.createDimension('regvar',n_regvar)

    longitude = f.createVariable('Longitude','f4','lon')
    latitude  = f.createVariable('Latitude','f4','lat')
    variable = f.createVariable('Variable','i4','var')
    time = f.createVariable('Time','i4','tim')
    regvar = f.createVariable('Regvar','i4','regvar')

    likelihood = f.createVariable('Likelihood','f4',('tim','lat','lon'),zlib=True)
    predictmean = f.createVariable('Predictmean','f4',('var','tim','lat','lon'),zlib=True)
    predictvariance = f.createVariable('Predictvariance','f4',('var','tim','lat','lon'),zlib=True)
    degreefreedom = f.createVariable('Degreefreedom','f4',('tim','lat','lon'),zlib=True)
    regrevariable = f.createVariable('Regrevariable','f4',('regvar','tim','lat','lon'),zlib=True)
    yanomaly = f.createVariable('Yanomaly','f4',('tim','lat','lon'),zlib=True)
    longitude[:] = lon
    latitude[:] = lat
    variable[:] = va
    time[:] = ti
    regvar[:] = rvar
    likelihood[:,:,:] = slik
    predictmean[:,:,:,:] = sm
    predictvariance[:,:,:,:] = sC
    degreefreedom[:,:,:] = snu 
    regrevariable[:,:,:,:] = Xpre
    yanomaly[:,:,:] = Yout
    f.close()




#%%
### using previous month prec for current ndvi  ## 0 using current prec for current ndvi
# NDVI dataset 
#root = '/Users/YaoZhang/'
#root = '/rigel/glab/users/zy2309/'
#root = "/global/scratch/yaozhang/Project/water_availability/"

root = '/gpfs/share/home/2106189133'
'''
def readcruncep(var):
    for i in range(5):
        pattern='/global/scratch/yaozhang/Data/CRUNCEP_v7_Mon/CRU_NCEP_V7_'+str(197+i)+'0_*'+var+'.*'
        var_file = glob.glob(pattern)
        with Dataset(var_file,'r') as fin:   # starts from 1979/7
            tmp = fin.variables[var][:]
        if i==0:
            dat = tmp[102:,:,:]
            print(dat)
        else:
            dat = np.concatenate((dat,tmp),along=0)
    return dat
'''

if ind_process_dlm:
    ndvi_file = '/global/scratch/yaozhang/Project/overshooting/Data/GIMMS3g_v1.mon.hd.growingseason.nc'

    with Dataset(ndvi_file,'r') as fin:   # starts from 1981/7
        ndvi = fin.variables["ndvi"][:]

    prec_file = '/global/scratch/yaozhang/Data/CRU/cru_ts4.04.1901.2019.pre.dat.nc'
    clou_file = '/global/scratch/yaozhang/Data/CRU/cru_ts4.04.1901.2019.cld.dat.nc'
    temp_file = '/global/scratch/yaozhang/Data/CRU/cru_ts4.04.1901.2019.tmp.dat.nc'
    with Dataset(prec_file,'r') as fin:   # starts from 1979/7
        prec = fin.variables["pre"][:][yr_ind:1380,:,:]
    with Dataset(clou_file,'r') as fin:   # starts from 1979/7
        clou = fin.variables["cld"][:][yr_ind:1380,:,:]
    with Dataset(temp_file,'r') as fin:   # starts from 1979/7
        temp = fin.variables["tmp"][:][yr_ind:1380,:,:]
    '''
    temp = readcruncep("tair")
    prec = readcruncep("prec")
    clou = readcruncep("solr")
    '''

    print(ndvi.shape)
    print(prec.shape) 
    fill_value = -999
    ndvi = np.where(np.logical_or(ndvi== -3000, ndvi== -9999), np.nan, ndvi)
    
    ndvi = ndvi/10000.0
    
    obs_ndvi = ndvi.shape[0]
    print(obs_ndvi)
    combined = np.concatenate((ndvi[(yr_ind-942):414,:,:],clou[24:(obs_ndvi+24),:,:],temp[24:(obs_ndvi+24),:,:],prec[0:(obs_ndvi+24),:,:]),axis=0)
    df = combined.reshape([obs_ndvi*4+24,360*720])
    print(df.shape)
    

    results =  parallelize_dataframe_multi(df, n_cores=20)
    outfile = root+'/analysis/DLM_GIMMS/multi_gimms'+st+'.DLM.CRU.spinup'+str(spin_up_years)+'.rep'+str(spin_up_rep)+'.delta'+str(delta)+'.nc'

    ### export to nc
    exportnc(results, outfile,n_regvar=4)

