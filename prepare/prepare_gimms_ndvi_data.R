### prepare gimms 3g NDVI dataset.
library(raster)
library(ncdf4)
library(abind)

setwd("~/Documents/Project/overshooting/")

tmp_f<-"./Data/CRU_TS_403_tmp.nc"
ncin<-nc_open(tmp_f)
tmp<-ncvar_get(ncin,varid="tmp") ## from 1979/7/1
nc_close(ncin)

vi_f<-"./Data/GIMMS3g_v1.mon.hd.highquality.nc"
ncin<-nc_open(vi_f)
vi<-ncvar_get(ncin,varid="ndvi") ## from 1981/7/1
nc_close(ncin)

good<-tmp[,,25:438]>0
good[good==0]<-NA
vi_good<-good*vi
vi_good[vi_good<500]<-NA
vi_good[is.na(vi_good)]<- -999

xtime<-ncdim_def("time",units = "months since 1981/01/01", vals=6+0:(dim(vi_good)[3]-1))
lon<-ncdim_def("longitude",units= "deg", vals=(-360:359+0.5)/2)
lat<-ncdim_def("latitude",units= "deg", vals=(-180:179+0.5)/2)
cli_var<-ncvar_def("ndvi",units="NA",dim = list(lon,lat,xtime),missval = -999,compression = 9)
ncout<-nc_create(paste0("~/Documents/Project/overshooting/Data/GIMMS3g_v1.mon.hd.growingseason.nc"),vars = list(cli_var))
ncvar_put(ncout,cli_var,vi_good)
nc_close(ncout)


