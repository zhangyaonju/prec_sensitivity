args = commandArgs(trailingOnly=F)
print(args)

myargs1 <-as.numeric(sub('-','',args[length(args)]))

print(myargs1)

library(Rmpfr)
#library(parallel)
# library(doParallel)

# mean PET 1671 mm/yr
setwd("/global/scratch/yaozhang/Project/water_availability/")

#setup parallel backend to use many processors
# cores=detectCores()
# cl <- makeCluster(cores[1]-1) #not to overload your computer

gammaincm<-function (x, a) 
{
  # if (!is.numeric(a) || !is.numeric(x)) 
  #   stop("All arguments must be real numbers.")
  if (length(a) > 1 || length(x) > 1) 
    stop("Arguments must be of length 1; function is not vectorized.")
  if (x == 0 && a == 0) 
    return(1)
  if (x == 0) 
    return(gamma(a))
  if (a < 0) 
    stop("Argument 'a' must be real and nonnegative.")
  if (x > 0) 
    xam <- -x + a * log(x)
  else xam <- -x + a * log(x + (0+0i))
  if (abs(xam) > 70000 || abs(a) > 50000) {
    warning("Arguments 'x' and/or 'a' are too large.")
    return(NA)
  }
  gin <- gim <- gip <- 0
  if (x == 0) {
    ga <- gamma(a)
    gim <- ga
    gip <- 0
  }
  else if (x <= 1 + a) {
    s <- 1/a
    r <- s
    for (k in 1:60) {
      r <- r * x/(a + k)
      s <- s + r
      if (abs(r/s) < 1e-15) 
        break
    }
    gin <- exp(xam) * s
    ga <- gamma(a)
    gip <- gin/ga
    gim <- ga - gin
  }
  else if (x > 1 + a) {
    t0 <- 0
    for (k in 60:1) {
      t0 <- (k - a)/(1 + k/(x + t0))
    }
    gim <- exp(xam)/(x + t0)
    ga <- gamma(a)
    gin <- ga - gim
    gip <- 1 - gim/ga
  }
  return(c(lowinc = Re(gin), uppinc = Re(gim), reginc = Re(gip)))
}


dx<-0.002
n<-c(0.35,0.43,0.45,0.5)
sh<-c(0.08,0.14,0.19,0.47)
sf<-c(0.35,0.56,0.65,1)
sw<-c(0.11,0.18,0.24,0.52)

root<-c(672,370,3140)
alpha<-c(7.4,4,16)
Delta<-c(0.588,0.385,0.799)

# gam<-array(NA,dim=c(4,3,3))
# for (i in 1:4){
#   for (j in 1:3){
#     for (k in 1:3){
#       gam[i,j,k]<-n[i]*root[k]*(sf[i]-sh[i])/alpha[j]
#     }
#   }
# }
# omeg=(sw-sh)/(sf-sh)
# 
# delta<-array(NA,dim=c(3,3))
# for (j in 1:3){
#   for (k in 1:3){
#     delta[j,k]<-Delta[k]/alpha[j]
#   }
# }

calculate_dydx<-function(x,st,rt,rd,cs,k){
  delta = Delta[cs]/alpha[rd]
  gamprime = n[st]*root[rt]*(sf[st]-sh[st])/alpha[rd]/(1-delta)
  omeg = (sw[st]-sh[st])/(sf[st]-sh[st])
  fracp = (1-delta)*exp(-delta)
  
  stat<-array(NA,dim=c(length(x),5))
  for (i in 1:length(x)){
    xprime<-(x[i]-1+fracp)/fracp
    xprimedx<-(x[i]+dx-1+fracp)/fracp
    if ((gamprime/xprime-1<0)|(gamprime/xprimedx-1)<0){
      next
    }
    # k is the percentage of stomata conductance
    if (k!=1){
      l1 = 1+0.38/(1+mpfr(exp(4/x[i]), precBits = 1000))
      l2 = 1+0.38/(1+mpfr(exp(4/(x[i]+dx)), precBits = 1000))
    }else{
      l1 = 1
      l2 = 1
    }
    k1 <- mpfr(gamprime/xprime-1, precBits = 1000)
    k2 <- mpfr(gamprime/xprime, precBits = 1000)
    k3 <- mpfr(gamprime/xprimedx-1, precBits = 1000)
    k4 <- mpfr(gamprime/xprimedx, precBits = 1000)
     
    # eq 10 *phi prime
    stat1 = 1-(gamprime^(k1)*exp(-gamprime))*xprime/(gamma(k2)-gammaincm(gamprime,k2)[2])
    # eq 8 f
    stat2 = ((gammaincm(gamprime*omeg,k2)[2]-gammaincm(gamprime,k2)[2])-
               omeg*gamprime*(gammaincm(gamprime*omeg,k1)[2]-gammaincm(gamprime,k1)[2]))/
      (gamma(k2)-gammaincm(gamprime,k2)[2])

    stat3 = 1-(gamprime^(k3)*exp(-gamprime))*xprimedx/(gamma(k4)-gammaincm(gamprime,k4)[2])
    stat4 = ((gammaincm(gamprime*omeg,k4)[2]-gammaincm(gamprime,k4)[2])-
               omeg*gamprime*(gammaincm(gamprime*omeg,k3)[2]-gammaincm(gamprime,k3)[2]))/
      (gamma(k4)-gammaincm(gamprime,k4)[2])
    stat5<-(stat3*stat4*l2*k-stat1*stat2*l1*k)/dx
    stat[i,]<-c(as.numeric(stat1),as.numeric(stat2),as.numeric(stat3),as.numeric(stat4),as.numeric(stat5))
  }
  return(stat)
}

x1<-c(5:9/50,4:19/20,10:49/10,25:85/5)
#x2<-x1/4.171532*4.130957
x2<-x1/4.159*4.1396

if (myargs1<=36){
  #soil type
  i = ceiling(myargs1/9)
  # root depth
  j = floor((myargs1-1)%%9/3)+1
  # rainfall depth #same as canopy storage
  raini = ((myargs1-1)%%9)%%3+1
  
  Lowco2 = calculate_dydx(x1,i,j,raini,raini,k=1)
  Highco2c3_no = calculate_dydx(x1,i,j,raini,raini,k=0.969)
  Highco2c4_no = calculate_dydx(x1,i,j,raini,raini,k=0.922)
  Highco2c3_pet = calculate_dydx(x2,i,j,raini,raini,k=0.969)
  Highco2c4_pet = calculate_dydx(x2,i,j,raini,raini,k=0.922)
  
  save(Lowco2,Highco2c3_no,Highco2c4_no,Highco2c3_pet,Highco2c4_pet,
       file=paste0("./analysis/minimalistic_model/minimalistic_model_gs_lai_",i,"_",j,"_",raini,".RData"))
}

# test the interception ---------------------------------------------------
if (myargs1>36){
  #soil type
  i = ceiling((myargs1-36)/3)
  # root depth
  j = 1
  # rainfall depth
  raini = 1
  # canopt storage
  canopys<-(myargs1-37)%%3+1
  
  Lowco2 = calculate_dydx(x1,i,j,raini,canopys,k=1)
  Highco2c3_no = calculate_dydx(x1,i,j,raini,canopys,k=0.969)
  Highco2c4_no = calculate_dydx(x1,i,j,raini,canopys,k=0.922)
  Highco2c3_pet = calculate_dydx(x2,i,j,raini,canopys,k=0.969)
  Highco2c4_pet = calculate_dydx(x2,i,j,raini,canopys,k=0.922)
  
  save(Lowco2,Highco2c3_no,Highco2c4_no,Highco2c3_pet,Highco2c4_pet,
       file=paste0("./analysis/minimalistic_model/minimalistic_model_gs_lai_",i,"_",j,"_",raini,"_",canopys,".RData"))
}



if (F){
  pdf("~/Dropbox/YAOZHANG/paper/2020_water_availability/figures/T_P.pdf",width=6,height=3)
  
  plot(log(1/(0.85*x1+0.15),10),log(lowco2[,1]*lowco2[,2],10),type="l",ylim=c(-2,0),axes=F,xlab="",ylab="")
  axis(1,tck=-0.04,at=-1:0,label=10^(-1:0),cex.lab=0.8)
  axis(1,tck=-0.02,at=log(c(3:10/100,2:10/10,2),10),label=c(0.03,rep(NA,16),2),cex.lab=0.8)
  axis(2,tck=-0.02,las=2,at=c(-2:1),cex.lab=0.8,lab=10^(-2:1))
  mtext(side=1,line=2.5,"P/PET")
  mtext(side=2,line=2.5,"T/P")
  dev.off()
  
  plot(log(0.85*x1+0.15,10),log(lowco2[,5]*(0.85*x1+0.15),10),ylim=c(-2,1))
  plot(log(0.85*x1+0.15,10),log(lowco2[,1]*lowco2[,2]-lowco2[,5]*(0.85*x1+0.15),10),ylim=c(-2,1))
  
  pdf("~/Dropbox/YAOZHANG/paper/2020_water_availability/figures/dT_dP.pdf",width=6,height=3)
  
  plot(log(1/(0.85*x1+0.15),10),lowco2[,1]*lowco2[,2]-lowco2[,5]*(0.85*x1+0.15),type="l",ylim=c(-1,1),axes=F,xlab="",ylab="")
  abline(h=0)
  axis(1,tck=-0.04,at=-1:0,label=10^(-1:0),cex.lab=0.8)
  axis(1,tck=-0.02,at=log(c(3:10/100,2:10/10,2),10),label=c(0.03,rep(NA,16),2),cex.lab=0.8)
  axis(2,tck=-0.02,las=2,at=c(-2:1),cex.lab=0.8,lab=10^(-2:1))
  mtext(side=1,line=2.5,"P/PET")
  mtext(side=2,line=2.5,expression(paste(partialdiff,"T/",partialdiff,"P",sep="")))
  dev.off()
  
  
  
  
  # calculate k for c3 and c4 platns ----------------------------------------
  deltaca = 30
  ca=354
  beta=0.6    # anthony walker 2020 new Phytologist
  beta=0
  k= (1+beta*deltaca/ca)/(1+deltaca/ca)
}
 



if (F){
  
  #### pdf of soil moisture
  
  x<-1:100/100
  px<-c()
  thet<-3
  gam=17.96
  
  for (i in 1:100){
    px[i]<-(x[i]^(gam/thet-1)*exp(-gam*x[i])*gam^(gam/thet))/(gamma(gam/thet)-gammainc(gam,gam/thet)[2])
  }
  plot(x,px)
  
  
  
  # optimality model -----------------------------------------------------------
  z=0.1 # km
  temp =293.15 ## deg C
  D = 1 #kP
  co2=c(340,400)
  Ko=27480*exp(36.38*(temp-298.15)/(298.15*8.314*temp))
  Kc =39.97*exp(79.43*(temp-298.15)/(298.15*8.314*temp))
  itastar = 1.0016/0.89
  Po = 21000*exp(-0.114*z)
  K = Kc*(1+Po/Ko)
  
  psi=sqrt(356.51*K/1.6/itastar)
  kai = psi/(psi+sqrt(D))
  wue =co2*(1-kai)/(1.6*D)
  
  
  plot(NA,xlim=c(-1.3,0.3),ylim=c(0,3),axes=F,xlab="",ylab="")
  axis(1,tck=-0.04,at=-1:0,label=10^(-1:0),cex.lab=0.8)
  axis(1,tck=-0.02,at=log(c(3:10/100,2:10/10,2),10),label=c(0.03,rep(NA,16),2),cex.lab=0.8)
  axis(2,tck=-0.02,las=2,at=c(0:3),cex.lab=0.8)
  box()
  mtext(side=2,line=2.5,"WUE")
  mtext(side=1,line=2.5,"P/PET")
  abline(h=wue[1],col="blue")
  abline(h=wue[2],col="red")
  
  
  plot(NA,xlim=c(-1.3,0.3),ylim=c(-0.5,0.5),axes=F,xlab="",ylab="")
  axis(1,tck=-0.04,at=-1:0,label=10^(-1:0),cex.lab=0.8)
  axis(1,tck=-0.02,at=log(c(3:10/100,2:10/10,2),10),label=c(0.03,rep(NA,16),2),cex.lab=0.8)
  axis(2,tck=-0.02,las=2,at=c(-3:3*0.1),cex.lab=0.8)
  box()
  lines(log(1/(0.85*x1+0.15),10),(highco2[,1]*highco2[,2]-highco2[,5]*(0.85*x2+0.15))*wue[2]-
          ((lowco2[,1]*lowco2[,2]-lowco2[,5]*(0.85*x1+0.15))*wue[1]),col="black")
  mtext(side=2,line=3.5,expression(paste(Delta,"(",partialdiff,"GPP/",partialdiff,"P)",sep="")))
  abline(h=0)
  
}

### test the relationship between aridity and c4 fraction
# dat<-as.matrix(read.table("./Downloads/ISLSCP_C4_1DEG_932/data/c4_percent_1d.asc"))
# dim(dat)
# image(dat)
# dat[dat<0]<-NA
# 
# lati<--90:89+0.5
# long<- -180:179+0.5
# library(ncdf4)
# lon<-ncdim_def('long',lati,units = "deg",longname = "longitude")
# lat<-ncdim_def("lati",long,units = "deg",longname = "latitude")
# c4pct<-ncvar_def("c4_pct",units = "",dim=list(lon,lat),longname = "percentage of C4")
# ncout<-nc_create("~/Documents/Data/c4pct.nc",vars = c4pct)
# ncvar_put(ncout,c4pct,dat[180:1,])
# nc_close(ncout)
# 
# library(raster)
# load("~/Documents/Project/water_availability/data/ai_HD.RData")
# ai_ras<-raster(ai_dat)
# ai_1d<-aggregate(ai_ras,c(2,2))
# image(ai_1d)
# ai<-as.matrix(ai_1d)
# image(ai)
# 
# dat[dat==0]<-NA
# com_dat<-cbind(as.vector(ai[,46:135]),as.vector(t(dat)[,135:46]))
# 
# plot(com_dat[,1],com_dat[,2])
# cor.test(com_dat[,1],com_dat[,2])
# 

if(F){
  # calculate the PET responses ---------------------------------------------
  D = 1   # kPa
  temp=20 # deg C
  u=2.7   # m/s
  Rn=10   # MJ m-2 day-1
  P=101   # kPa
  rhocp = 1.013*10^-3   # MJ/kg K
  lambda=2.45    # MJ kg-1
  mwratio=0.622   # NA
  gammapc = (rhocp*P)/(lambda*mwratio)
  delta<- 4098*(0.6108*exp((17.27*temp)/(temp+237.3)))/(temp+237.3)^2
  ep <-(0.408*delta*Rn+gammapc*(900/(temp+273))*u*D)/(delta+gammapc*(1+u*(0.34+2.4*10^-4*(co2-300))))
}
