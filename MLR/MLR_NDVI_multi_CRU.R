#multi linear regression 3models

args = commandArgs(trailingOnly=F)
print(args)
#myargs1 <-sub('-','',args[length(args)-2])
#myargs2 <-sub('-','',args[length(args)-1])
myargs3 <-sub('-','',args[length(args)])

setwd('/global/scratch/yaozhang/Project/water_availability/')
#setwd("~/Documents/Project/legacy_water_effect/")

prec_v<-as.numeric(myargs3)

library(raster)
library(ncdf4)
library(abind)
library(parallel)
#library(aperm)
### get the detrended and deseasonalized 
deseason<-function(vec){
  if(sum(!is.na(vec))<100){
    return(rep(NA,length(vec)))
  }
  n_years<-floor(length(vec)/12)
  extra<-length(vec)-n_years*12
  vec_year<-vec[1:(n_years*12)]
  dim(vec_year)<-c(12,n_years)
  msc<-apply(vec_year,1,mean,na.rm=T)
  deseasoned<-vec-c(rep(msc,n_years),msc[1:extra])
  x<-1:length(vec)
  detr<-lm(deseasoned~x)
  coefs<-detr$coefficients
  return(deseasoned-(coefs[1]+coefs[2]*x))
}

### get the pre-season precipitation
get_prec<-function(ts_pre){    # if ndvi version is v0 then only 390 observations. if ndvi version is v1 then 414 observation
  deseason<-function(vec){
    if(sum(!is.na(vec))<100){
      return(rep(NA,length(vec)))
    }
    n_years<-floor(length(vec)/12)
    extra<-length(vec)-n_years*12
    vec_year<-vec[1:(n_years*12)]
    dim(vec_year)<-c(12,n_years)
    msc<-apply(vec_year,1,mean,na.rm=T)
    deseasoned<-vec-c(rep(msc,n_years),msc[1:extra])
    x<-1:length(vec)
    detr<-lm(deseasoned~x)
    coefs<-detr$coefficients
    return(deseasoned-(coefs[1]+coefs[2]*x))
  }
  obs_n<-length(ts_pre)
  #dts_pre<-deseason(ts_pre)
  ts_pre0<-deseason(ts_pre[25:(414+24)-18])
  ts_pre01<-deseason((ts_pre[24:(414+23)-18]+ts_pre[25:(414+24)-18])/2)
  ts_pre02<-deseason((ts_pre[23:(414+22)-18]+ts_pre[24:(414+23)-18]+ts_pre[25:(414+24)-18])/3)
  ts_pre03<-deseason((ts_pre[22:(414+21)-18]+ts_pre[23:(414+22)-18]+ts_pre[24:(414+23)-18]+ts_pre[25:(414+24)-18])/4)
  ts_pre1<-deseason(ts_pre[24:(414+23)-18])
  ts_pre12<-deseason((ts_pre[23:(414+22)-18]+ts_pre[24:(414+23)-18])/2)
  return(c(ts_pre0,ts_pre01,ts_pre02,ts_pre03,ts_pre1,ts_pre12))
}


ncin<-nc_open(paste0("./data/GIMMS3g_v1.mon.hd.highquality.nc"))
ndvi<-ncvar_get(ncin,varid="ndvi")
ndvi[ndvi<0]<-NA
nc_close(ncin)

ncin<-nc_open("/global/scratch/yaozhang/Data/CRU/cru_ts4.04.1901.2019.pre.dat.nc")
prec<-ncvar_get(ncin,varid="pre")[,,961:1380]    # from 198101 to 201512
nc_close(ncin)

ncin<-nc_open("/global/scratch/yaozhang/Data/CRU/cru_ts4.04.1901.2019.cld.dat.nc")
clou<-ncvar_get(ncin,varid="cld")[,,961:1380]    # from 198101 to 201512
nc_close(ncin)

ncin<-nc_open("/global/scratch/yaozhang/Data/CRU/cru_ts4.04.1901.2019.tmp.dat.nc")
temp<-ncvar_get(ncin,varid="tmp")[,,961:1380]    # from 198101 to 201512
nc_close(ncin)

# # use CRUNCEP-V7
# read_cruncep_data<-function(var,varname){
#   files<-list.files("/global/scratch/yaozhang/Data/CRUNCEP_v7_mon/",pattern=var,full.names=T)[9:12]
#   dat<-array(NA,dim=c(432,720,360))
#   for (i in 1:length(file)){
#     ncin<-nc_open(files[i])
#     temp<-ncvar_get(ncin,varid = varname)
#     dat[1:dim(temp)[1]+i*120-120,,]<-temp
#   }
#   return(dat)
# }
# 
# prec<-read_cruncep_data("prec","prec")[1:420,,]
# clou<-read_cruncep_data("solr","solr")[1:420,,]
# temp<-read_cruncep_data("tair","tair")[1:420,,]

### ndvi is detrended deseasonalized. prec is deseasonalized
ndvi = ndvi/10000

cl<-makeCluster(6)

dseason_ndvi<-parApply(cl,ndvi,c(1,2),deseason)

prec_ano<-parApply(cl,prec,c(1,2),get_prec)
clou_ano<-parApply(cl,clou[,,7:420],c(1,2),deseason)
temp_ano<-parApply(cl,temp[,,7:420],c(1,2),deseason)



select_best_model<-function(vec){
  MLR_model<-function(vec){
    n_obs<-(length(vec))/4
    if (sum(!is.na(vec[1:n_obs]))<100|sum(!is.na(vec[(1+n_obs):(length(vec))]))<100){
      return(rep(NA,6))
    }
    dndvi<-vec[1:n_obs]
    Xclou<-vec[(n_obs+1):(n_obs*2)]
    Xtemp<-vec[(n_obs*2+1):(n_obs*3)]
    Xprec<-vec[(n_obs*3+1):(n_obs*4)]
    ### conduct the multi linear regression and get the coeficient
    reg<-lm(dndvi[2:n_obs]~dndvi[1:(n_obs-1)]+Xprec[2:n_obs]+Xtemp[2:n_obs]+Xclou[2:n_obs])
    res<-c(reg$coefficients,summary(reg)$r.square)
    return(res)
  }
  
  n_obs<-length(vec)/9
  ndvi<-vec[1:n_obs]
  clou<-vec[n_obs+1:n_obs]
  temp<-vec[n_obs*2+1:n_obs]
  prec<-vec[n_obs*3+1:(n_obs*6)]
  dim(prec)<-c(n_obs,6)
  res<-array(NA, dim=c(6,6)) ## 5 coef+rsq * 6 scenarios
  for (i in 1:6){
    res[,i]<-MLR_model(c(ndvi,clou,temp,prec[,i]))
  }
  max_ind<-which.max(res[6,])
  if (length(max_ind)==0){
    return(rep(NA,8))
  }else{
    return(c(res[,max_ind],max_ind,res[6,5]))
  }
}

### model 1  breakpoint 1998,1999
dim(prec_ano)<-c(414,6,720,360)
prec_ano_fh<-prec_ano[1:210,,,]
prec_ano_sh<-prec_ano[211:414,,,]
dim(prec_ano_fh)<-c(210*6,720,360)
dim(prec_ano_sh)<-c(204*6,720,360)

cli_fh<-abind(clou_ano[1:210,,],temp_ano[1:210,,],along=1)
clipre_fh<-abind(cli_fh,prec_ano_fh,along=1)
comb1 = abind(dseason_ndvi[1:210,,],clipre_fh,along=1)

cli_sh<-abind(clou_ano[211:414,,],temp_ano[211:414,,],along=1)
clipre_sh<-abind(cli_sh,prec_ano_sh,along=1)
comb2 = abind(dseason_ndvi[211:414,,],clipre_sh,along=1)


first_half<-parApply(cl,comb1,c(2,3),select_best_model)
second_half<-parApply(cl,comb2,c(2,3),select_best_model)

save(first_half,second_half,file=paste0('./analysis/MLR_GIMMS/mlr_separate_period_multi_CRU_best_model.RData'))



### model 2   breakpoint 1997,1998
dim(prec_ano)<-c(414,6,720,360)
prec_ano_fh<-prec_ano[1:198,,,]
prec_ano_sh<-prec_ano[199:414,,,]
dim(prec_ano_fh)<-c(198*6,720,360)
dim(prec_ano_sh)<-c(216*6,720,360)

cli_fh<-abind(clou_ano[1:198,,],temp_ano[1:198,,],along=1)
clipre_fh<-abind(cli_fh,prec_ano_fh,along=1)
comb1 = abind(dseason_ndvi[1:198,,],clipre_fh,along=1)

cli_sh<-abind(clou_ano[199:414,,],temp_ano[199:414,,],along=1)
clipre_sh<-abind(cli_sh,prec_ano_sh,along=1)
comb2 = abind(dseason_ndvi[199:414,,],clipre_sh,along=1)


first_half<-parApply(cl,comb1,c(2,3),select_best_model)
second_half<-parApply(cl,comb2,c(2,3),select_best_model)

save(first_half,second_half,file=paste0('./analysis/MLR_GIMMS/mlr_separate_period_multi_CRU_best_model_1997.RData'))


### model 2   breakpoint 1997,1998
dim(prec_ano)<-c(414,6,720,360)
prec_ano_fh<-prec_ano[1:222,,,]
prec_ano_sh<-prec_ano[223:414,,,]
dim(prec_ano_fh)<-c(222*6,720,360)
dim(prec_ano_sh)<-c(192*6,720,360)

cli_fh<-abind(clou_ano[1:222,,],temp_ano[1:222,,],along=1)
clipre_fh<-abind(cli_fh,prec_ano_fh,along=1)
comb1 = abind(dseason_ndvi[1:222,,],clipre_fh,along=1)

cli_sh<-abind(clou_ano[223:414,,],temp_ano[223:414,,],along=1)
clipre_sh<-abind(cli_sh,prec_ano_sh,along=1)
comb2 = abind(dseason_ndvi[223:414,,],clipre_sh,along=1)


first_half<-parApply(cl,comb1,c(2,3),select_best_model)
second_half<-parApply(cl,comb2,c(2,3),select_best_model)

save(first_half,second_half,file=paste0('./analysis/MLR_GIMMS/mlr_separate_period_multi_CRU_best_model_1999.RData'))

