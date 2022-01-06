
# calculate the precipitation sensitivity for all tropical regions --------

#multi linear regression 3models

args = commandArgs(trailingOnly=F)
print(args)

myargs3 <-sub('-','',args[length(args)])
args<-as.numeric(myargs3)

setwd('/global/scratch/yaozhang/')
#setwd("~/Documents/Project/legacy_water_effect/")

library(raster)
library(ncdf4)
library(abind)
library(parallel)
#library(aperm)
### get the detrended and deseasonalized 
deseason<-function(vec){
  detrended<-function(vec){
    if (sum(is.na(vec))>10){
      return(rep(NA,length(vec)))
    }
    x<-1:length(vec)
    reg<-lm(vec~x)
    return(vec-x*reg$coefficients[2]-reg$coefficients[1])
  }
  
  n_years<-floor(length(vec)/12)
  vec_year<-vec[1:(n_years*12)]
  dim(vec_year)<-c(12,n_years)
  msc<-apply(vec_year,1,mean,na.rm=T)
  deseasoned<-vec-rep(msc,n_years)
  return(detrended(deseasoned))
}



MLR_model4<-function(vec){
  n_obs<-(length(vec))/4
  if (sum(!is.na(vec[1:n_obs]))<50|sum(is.na(vec[(1+n_obs):(length(vec))]))>50){
    return(rep(NA,6))
  }
  dndvi<-vec[1:n_obs]
  Xprec<-vec[(n_obs+1):(n_obs*2)]
  Xclou<-vec[(n_obs*2+1):(n_obs*3)]
  Xtemp<-vec[(n_obs*3+1):(n_obs*4)]
  ### conduct the multi linear regression and get the coeficient
  reg<-lm(dndvi[2:n_obs]~dndvi[1:(n_obs-1)]+Xprec[1:(n_obs-1)]+Xclou[2:n_obs]+Xtemp[2:n_obs])
  res<-c(reg$coefficients,summary(reg)$r.square)
  return(res)
}



# MsTMIP GPP --------------------------------------------------------------


cl<-makeCluster(10)

GPP_file<-list.files("./Data/MsTMIP_v2/GPP/",full.names = T)[args]
model_name<-strsplit(basename(GPP_file),"_Monthly_GPP.nc4")[[1]]
ncin<-nc_open(GPP_file)
GPP<-ncvar_get(ncin,varid="GPP")*86400000
GPP[GPP<0]<-NA
nc_close(ncin)

ncin<-nc_open("./Data/MsTMIP_climate/CRU_NCEP_monthly_rain_1901_2010.nc")
prec<-ncvar_get(ncin,varid="rain")[,360:1,]*10   # convert cm/mon to mm/mon
nc_close(ncin)

ncin<-nc_open("./Data/MsTMIP_climate/CRU_NCEP_monthly_swdown_1901_2010.nc")
clou<-ncvar_get(ncin,varid="swdown")[,360:1,]/1e6   # convert j/m2/mon to MJ/m2/mon
nc_close(ncin)

ncin<-nc_open("./Data/MsTMIP_climate/CRU_NCEP_monthly_tair_1901_2010.nc")
tair<-ncvar_get(ncin,varid="tair")[,360:1,]   # degree K
nc_close(ncin)

## every 10 year combine prec and GPP and calculate the precipitation
all_periods_stat<-array(NA,dim=c(11,6,720,360))

for (i in 1:11){
  ### GPP is detrended deseasonalized. prec is deseasonalized
  dseason_GPP<-parApply(cl,GPP[,,(i*120-119):(i*120)],c(1,2),deseason)
  prec_ano<-parApply(cl,prec[,,(i*120-119):(i*120)],c(1,2),deseason)
  clou_ano<-parApply(cl,clou[,,(i*120-119):(i*120)],c(1,2),deseason)
  tair_ano<-parApply(cl,tair[,,(i*120-119):(i*120)],c(1,2),deseason)
  #dim(prec_ano)<-c(114,6,720,360)
  cm<-abind(prec_ano,clou_ano,along=1)
  cm<-abind(cm,tair_ano,along=1)
  comb1 <- abind(dseason_GPP,cm,along=1)
  all_periods_stat[i,,,]<-parApply(cl,comb1,c(2,3),MLR_model4)
}

save(all_periods_stat,file=paste0('./Project/water_availability/analysis/MsTMIP_gpp_v2/partial_mlr_MsTMIP_',model_name,'_best_model.RData'))

# 
# 
# if (F){
#   lai_file<-list.files("./Data/MsTMIP_v1/LAI/",full.names = T)[args]
#   model_name<-strsplit(basename(lai_file),"_Monthly_LAI.nc4")[[1]]
#   ncin<-nc_open(lai_file)
#   lai<-ncvar_get(ncin,varid="LAI")
#   lai[lai<0]<-NA
#   nc_close(ncin)
#   
#   ncin<-nc_open("./Data/MsTMIP_climate/CRU_NCEP_monthly_rain_1901_2010.nc")
#   prec<-ncvar_get(ncin,varid="rain")[,360:1,]*10   # convert cm/mon to mm/mon
#   nc_close(ncin)
#   
#   ## every 10 year combine prec and lai and calculate the precipitation
#   all_periods_stat<-array(NA,dim=c(11,5,720,360))
#   
#   for (i in 1:11){
#     lai_period<-lai[,,(i*120-119):(i*120)]
#     pre_period<-prec[,,(i*120-119):(i*120)]
#     ### lai is detrended deseasonalized. prec is deseasonalized
#     dseason_lai<-apply(lai_period,c(1,2),deseason)
#     prec_ano<-apply(pre_period,c(1,2),get_prec)
#     #dim(prec_ano)<-c(114,6,720,360)
#     
#     comb1 <- abind(dseason_lai[7:120,,],prec_ano,along=1)
#     all_periods_stat[i,,,]<-apply(comb1,c(2,3),select_best_model)
#   }
#   
#   save(all_periods_stat,file=paste0('./Project/water_availability/analysis/MsTMIP_lai/mlr_MsTMIP_',model_name,'_best_model.RData'))
#   
# }

