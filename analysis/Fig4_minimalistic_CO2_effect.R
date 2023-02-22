##### Figure 4
# get all files together
if (F){
  setwd("/global/scratch/yaozhang/Project/water_availability/")
  
  lowco2<-array(NA,dim=c(122,5,4,9))
  highco2c3_no<-array(NA,dim=c(122,5,4,9))
  highco2c4_no<-array(NA,dim=c(122,5,4,9))
  highco2c3_pet<-array(NA,dim=c(122,5,4,9))
  highco2c4_pet<-array(NA,dim=c(122,5,4,9))
  
  for (i in 1:4){
    for (j in 1:3){
      for (k in 1:3){
        fi<-paste0("./analysis/minimalistic_model/minimalistic_model_gs_lai_",i,"_",j,"_",k,".RData")
        load(fi)
        lowco2[,,i,j*3-3+k]<-Lowco2
        highco2c3_no[,,i,j*3-3+k]<-Highco2c3_no
        highco2c4_no[,,i,j*3-3+k]<-Highco2c4_no
        highco2c3_pet[,,i,j*3-3+k]<-Highco2c3_pet
        highco2c4_pet[,,i,j*3-3+k]<-Highco2c4_pet
      }
    }
  }
  
  save(lowco2,highco2c3_no,highco2c4_no,highco2c3_pet,highco2c4_pet,
       file="./analysis/minimalistic_model_gs_lai.RData")
}

mean_std<-function(vec){
  return(c(mean(vec,na.rm=T),sd(vec,na.rm=T)))
}
setwd("~/Documents/Project/water_availability/")
load("./analysis/minimalistic_model_gs_lai.RData")
### process the data locally
# get the mean sensitivity
x1<-c(5:9/50,4:19/20,10:49/10,25:85/5)
x2<-x1/4.159*4.1396
#l1 = 0.38/(1+exp(1/(0.85*0.25*x1+0.15)))+1
fracp = 0.85
xprime<-(x1-1+fracp)/fracp
l1 = 1+0.38/(1+exp(4/x1))

lowco2_mean_std<-apply(lowco2[,1,,]*lowco2[,2,,]-lowco2[,5,,]*x1,c(1,2),mean_std)

# change of conductance due to CO2 for C3 and C4
only_co2_on_gs_c3_all<-(lowco2[,1,,]*lowco2[,2,,]-lowco2[,5,,]*x1)*0.969*fracp-
  ((lowco2[,1,,]*lowco2[,2,,]-lowco2[,5,,]*x1))*fracp
only_co2_on_gs_c3<-apply(only_co2_on_gs_c3_all,c(1),mean_std)

only_co2_on_gs_c4_all<-(lowco2[,1,,]*lowco2[,2,,]-lowco2[,5,,]*x1)*0.922*fracp-
  ((lowco2[,1,,]*lowco2[,2,,]-lowco2[,5,,]*x1))*fracp
only_co2_on_gs_c4<-apply(only_co2_on_gs_c4_all,c(1),mean_std)
## condutance with lai for c3 and c4
only_co2_on_gslai_c3_all<-(highco2c3_no[,1,,]*highco2c3_no[,2,,]*l1*0.969-highco2c3_no[,5,,]*x1)*fracp-
  ((lowco2[,1,,]*lowco2[,2,,]-lowco2[,5,,]*x1))*fracp
only_co2_on_gslai_c3<-apply(only_co2_on_gslai_c3_all,c(1),mean_std)

only_co2_on_gslai_c4_all<-(highco2c4_no[,1,,]*highco2c4_no[,2,,]*l1*0.922-highco2c4_no[,5,,]*x1)*fracp-
  ((lowco2[,1,,]*lowco2[,2,,]-lowco2[,5,,]*x1))*fracp
only_co2_on_gslai_c4<-apply(only_co2_on_gslai_c4_all,c(1),mean_std)

## condutance with lai and lai/T for c3 and c4
only_co2_on_gslaigs_c3_all<-(highco2c3_no[,1,,]*highco2c3_no[,2,,]*l1*0.969-highco2c3_no[,5,,]*x1)/0.969*fracp-
  ((lowco2[,1,,]*lowco2[,2,,]-lowco2[,5,,]*x1))*fracp
only_co2_on_gslaigs_c3<-apply(only_co2_on_gslaigs_c3_all,c(1),mean_std)

only_co2_on_gslaigs_c4_all<-(highco2c4_no[,1,,]*highco2c4_no[,2,,]*l1*0.922-highco2c4_no[,5,,]*x1)/0.922*fracp-
  ((lowco2[,1,,]*lowco2[,2,,]-lowco2[,5,,]*x1))*fracp
only_co2_on_gslaigs_c4<-apply(only_co2_on_gslaigs_c4_all,c(1),mean_std)

## condutance with lai and lai/T for c3 and c4 pet
only_co2_on_gslaigspet_c3_all<-(highco2c3_pet[,1,,]*highco2c3_pet[,2,,]*l1*0.969-highco2c3_pet[,5,,]*x1)/0.969*fracp-
  ((lowco2[,1,,]*lowco2[,2,,]-lowco2[,5,,]*x1))*fracp
only_co2_on_gslaigspet_c3<-apply(only_co2_on_gslaigspet_c3_all,c(1),mean_std)

only_co2_on_gslaigspet_c4_all<-(highco2c4_pet[,1,,]*highco2c4_pet[,2,,]*l1*0.922-highco2c4_pet[,5,,]*x1)/0.922*fracp-
  ((lowco2[,1,,]*lowco2[,2,,]-lowco2[,5,,]*x1))*fracp
only_co2_on_gslaigspet_c4<-apply(only_co2_on_gslaigspet_c4_all,c(1),mean_std)

# how does CO2 affect 
#load("~/Documents/Project/water_availability/analysis/minimalistic_model_gs_lai.RData")
#load("~/Documents/Project/water_availability/analysis/minimalistic_model.RData")
library(rcolors)
soilcol<-rcolors$MPL_copper[c(120,90,60,30)]
scenecol<-rcolors$GMT_paired[c(2,6,10,12)]

  
pdf("~/Documents/paper/=2020_water_availability/figures_P/Figure4_dT_dP_co2_rr.pdf",width=6.5,height=6)

par(mfrow=c(2,1),mar=c(2.5,4.7,1,1),mgp=c(3,0.4,0))


plot(NA,xlim=c(-1.3,0.3),ylim=c(-0.2,1),axes=F,xlab="",ylab="")
axis(1,tck=-0.04,at=log(c(0.05,0.1,0.5,1,2),10),label=c(0.05,0.1,0.5,1,2),cex.lab=1)
axis(1,tck=-0.02,at=log(c(3:10/100,2:10/10,2),10),label=c(0.03,rep(NA,16),2),cex.lab=1)
axis(2,tck=-0.02,las=2,at=c(-1:5*0.2),cex.lab=1)
abline(h=0,lty=2)
box()
for (i in 1:4){
  lines(log(1/x1,10),(lowco2[,1,i,1]*lowco2[,2,i,1]-lowco2[,5,i,1]*x1)*fracp,lwd=1.5,col=soilcol[i])
  for (j in 2:9){
    lines(log(1/x1,10),(lowco2[,1,i,j]*lowco2[,2,i,j]-lowco2[,5,i,j]*x1)*fracp,lwd=0.7,col=adjustcolor(soilcol[i],alpha.f = 0.3))
  }
}
mtext(side=1,line=1.3,"Aridity index")
mtext(side=2,line=3.1,"LAI sensi. to Prec.")
mtext(side=2,line=2.,expression(paste("(",partialdiff,"LAI/",partialdiff,"P)",sep="")))
mtext(side=2,line=3.7,las=2,padj = -9.5,letters[1],font=2)
legend("topleft",c("sand","sandy loam","loam","clay"),col=soilcol,lwd=rep(1.5,4),bty='n',x.intersp = 0.7,y.intersp = 0.8,cex=1)

plot(NA,xlim=c(-1.3,0.3),ylim=c(-0.1,0.1),axes=F,xlab="",ylab="")

lines(log(1/x1,10),only_co2_on_gs_c3[1,],lwd=1.5,col=scenecol[1])
polygon(c(log(1/x1,10),rev(log(1/x1,10))),
        c(only_co2_on_gs_c3[1,]+only_co2_on_gs_c3[2,],rev(only_co2_on_gs_c3[1,]-only_co2_on_gs_c3[2,])),col=adjustcolor(scenecol[1],alpha.f = 0.5),border = F)
#lines(log(1/(0.85*x1+0.15),10),only_co2_on_gs_c4[1,]*0.85,lwd=1,lty=3,col=scenecol[1])
lines(log(1/x1,10),only_co2_on_gslai_c3[1,],lwd=1.5,col=scenecol[2])
polygon(c(log(1/x1,10),rev(log(1/x1,10))),
        c(only_co2_on_gslai_c3[1,]+only_co2_on_gslai_c3[2,],rev(only_co2_on_gslai_c3[1,]-only_co2_on_gslai_c3[2,])),col=adjustcolor(scenecol[2],alpha.f = 0.5),border = F)

lines(log(1/x1,10),only_co2_on_gslaigs_c3[1,],lwd=1.5,col=scenecol[3])
polygon(c(log(1/x1,10),rev(log(1/x1,10))),
        c(only_co2_on_gslaigs_c3[1,]+only_co2_on_gslaigs_c3[2,],rev(only_co2_on_gslaigs_c3[1,]-only_co2_on_gslaigs_c3[2,])),col=adjustcolor(scenecol[3],alpha.f = 0.5),border = F)

# lines(log(1/(0.85*x1+0.15),10),only_co2_on_gslaigspet_c3[1,]*0.85,lwd=1.5,col=scenecol[4])
# polygon(c(log(1/(0.85*x1+0.15),10),rev(log(1/(0.85*x1+0.15),10))),
#         c(only_co2_on_gslaigspet_c3[1,]+only_co2_on_gslaigspet_c3[2,],rev(only_co2_on_gslaigspet_c3[1,]-only_co2_on_gslaigspet_c3[2,]))*0.85,col=adjustcolor(scenecol[4],alpha.f = 0.5),border = F)

axis(1,tck=-0.04,at=log(c(0.05,0.1,0.5,1,2),10),label=c(0.05,0.1,0.5,1,2),cex.lab=0.8)
axis(1,tck=-0.02,at=log(c(3:10/100,2:10/10,2),10),label=c(0.03,rep(NA,16),2),cex.lab=0.8)
#axis(2,tck=-0.02,las=2,at=c(-5:5*0.005),labels = rep("",11),cex.lab=0.8)
axis(2,tck=-0.02,las=2,at=c(-2:2*0.05),cex.lab=0.8)
box()
abline(h=0)
mtext(side=1,line=1.3,"Aridity index")
mtext(side=2,line=3.1,"Change in LAI sensi. to Prec.")
mtext(side=2,line=2,expression(paste("(",Delta,"(",partialdiff,"LAI/",partialdiff,"P))",sep="")))
mtext(side=2,line=3.7,las=2,padj = -9.5,letters[2],font=2)
#mtext(side=3,line=0,at=log(0.05,10),expression(paste(""%*%"1e-3",sep="")),cex=0.7)
legend("bottomleft",legend=c(expression(paste("g"[s]," on ",partialdiff,"E"[T],"/",partialdiff,"P",sep="")),
                      expression(paste("g"[s],"+LAI on ",partialdiff,"E"[T],"/",partialdiff,"P",sep="")),
                      expression(paste("g"[s],"+LAI on ",partialdiff,"E"[T],"/",partialdiff,"P+g"[s]," on ",partialdiff,"LAI/",partialdiff,"E"[T],sep=""))),
                      #expression(paste("g"[s],"+LAI on ",partialdiff,"T/",partialdiff,"P+g"[s]," on ",partialdiff,"LAI/",partialdiff,"T+E"[P],sep=""))),
                   col=scenecol[1:3],lwd=rep(1.5,3),bty='n',x.intersp = 0.7,y.intersp = 0.8,cex=0.8)
#legend("bottomright",legend=c("10%","15%","20%"),lty=c(1,2,3),lwd=rep(1.5,3),x.intersp = 0.7,y.intersp = 0.8,cex=0.8,bty='n')

dev.off()



# get the 
# y<-highco2c3_no[,3,1,1]*highco2c3_no[,4,1,1]*l2*0.95-(highco2c3_no[,3,1,1]*highco2c3_no[,4,1,1]*l2*0.95-highco2c3_no[,1,1,1]*highco2c3_no[,2,1,1]*l1*0.95)/dx*(0.85*x2+0.15)#-highco2c3_no[,5,1,1]
# plot(1/x1,y)
















###get the crosssing aridity
orix<-1/(0.85*x1+0.15)
crossing<-array(NA,dim=c(4,9))
for (i in 1:4){
  for (j in 1:9){
    diff<-(highco2[,1,i,j]*highco2[,2,i,j]-highco2[,5,i,j]*(0.85*x2+0.15))-
  ((lowco2[,1,i,j]*lowco2[,2,i,j]-lowco2[,5,i,j]*(0.85*x1+0.15)))
    crossing[i,j]<-orix[which(diff>0)[1]]
  }
}


