#### get average contribution, sensitivity, trend and per-biome statistics
args = commandArgs(trailingOnly=F)
print(args)


myargs1 <-sub('-','',args[length(args)-3])
myargs2 <-sub('-','',args[length(args)-2])
myargs3 <-sub('-','',args[length(args)-1])
myargs4 <-sub('-','',args[length(args)])

spinup<-as.numeric(myargs1)
spinrep<-as.numeric(myargs2)
delta<-c(0.95,0.98,0.99,0.995,0.999)[as.numeric(myargs3)+1]

setwd('/gpfs/share/home/2106189133/Project/water_availability/')

library(raster)
library(ncdf4)
library(abind)
library(trend)
library(parallel)

cl<-makeCluster(6)
clusterEvalQ(cl, {
  library(trend)
})

nvar<-c("EVI","NDVI")[as.numeric(myargs4)]
ind<-"_2015"
## deseasonlized and detrended variables are included
scenario<-paste0("spinup",spinup,".rep",spinrep,".delta",delta)
ncin<-nc_open(paste0("./analysis/DLM_MODIS/multi_modis",ind,".",nvar,".DLM.CRU.",scenario,".nc"))
predictmean<-ncvar_get(ncin,varid="Predictmean")
Xpre<-ncvar_get(ncin,varid="Regrevariable")
nc_close(ncin)
xdim = dim(Xpre)[4]

y_pre_ndvi<-predictmean[,,1:(180+spinup*12*spinrep),3]
y_prec_12<-predictmean[,,1:(180+spinup*12*spinrep),4]
rm(predictmean)

get_mean_trend<-function(vec){
  if (sum(is.na(vec))<100){
    sen<-sens.slope(vec)
    return(c(mean(vec,na.rm=T),sen$estimates,sen$p.value))
  }else{
    return(c(NA,NA,NA))
  }
}


# get the trend and mean --------------------------------------------------
ndvi_sens_mean_trend<-parApply(cl,y_pre_ndvi[,,(1+spinup*12*spinrep):(180+spinup*12*spinrep)],c(1,2),get_mean_trend)
prec_sens_mean_trend<-parApply(cl,y_prec_12[,,(1+spinup*12*spinrep):(180+spinup*12*spinrep)],c(1,2),get_mean_trend)

save(ndvi_sens_mean_trend,prec_sens_mean_trend,file=paste0("./analysis/DLM_MODIS/multi_MODIS",ind,"_",nvar,"_CRU_coef_sensi_mean_trend_",scenario,".RData"))


