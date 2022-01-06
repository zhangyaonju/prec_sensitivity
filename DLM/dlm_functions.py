import numpy as np
import scipy.linalg
from scipy.stats import t as tdstr
from scipy.stats import norm
from scipy import signal
from scipy.interpolate import interp1d
#from random import randint
#import pandas as pd
import random

## note for the variables
# Prior:   [class]  the initial state (t=0) for all variance.    
# m:       [vector] mean for the theta estiate, time varying
# C:       [matrix] variance matrix for theta, time varying
# S:       [number] observation noise for Y 
# nu:      [number] degree of freedom.
# X:       [matrix] regressor variables, including 1: y[t-1], prec[t], prec[t-1:t-2], prec[t-3:t-5],prec[t-6:t-11] for each time step. 
# Y:       [vector] observational time series of Y.
# delta:   [vector] 0.98, inflation rate
# rseas:   [vector] 
# pseas:   [number] number of seasonal harmonic function (frequency)
# nseas:   [number] dimension of the seasonal component
# ntrend:  [number] dimension of trend component
# nregn:   [number] number of regression variables.
# period:  [number] number of obs per year.
# itrend:  [vector] index for the trend component in the full model
# ireg:    [vector] index for the regression component in the full model  
# iseas:   [vector] index for the seasonal component in the full model  

class Prior:
    def __init__(self, m, C, S, nu):
        self.m = m
        self.C = C
        self.S = S
        self.nu = nu

class Model:
    def __init__(self,Y,X,rseas,delta):
        self.Y = Y
        self.X = X
        self.rseas = rseas
        dd = np.ones(4)*delta
        dd[0] = 0.999
        dd[2] = 0.999
        self.delta = dd
        ntrend = 2;nregn = X.shape[1]; pseas = len(rseas);nseas = pseas*2;
        m = np.zeros([ntrend+nregn+nseas,1])
        C = scipy.linalg.block_diag(1*np.eye(ntrend),1*np.eye(nregn),1*np.eye(nseas))
        S = np.power(0.2,2); nu = ntrend+nregn+pseas;
        pr = Prior(m,C,S,nu)
        self.prior = pr
        
        
def forwardFilteringM(Model):
    Y = Model.Y
    X = Model.X
    rseas = Model.rseas
    delta = Model.delta
    Prior = Model.prior
    period = 12   #number of obs per year.
    deltrend = delta[0];delregn = delta[1];delseas = delta[2];delvar = delta[3]
    Ftrend = np.array([[1],[0]]);ntrend = len(Ftrend); Gtrend = np.array([[1,1],[0,1]]);itrend = np.arange(0,ntrend)
    nregn = X.shape[1];Fregn = np.zeros([nregn,1]);Gregn=np.eye(nregn);iregn = np.arange(ntrend,ntrend+nregn)
    pseas = len(rseas);nseas = pseas*2;iseas = np.arange(ntrend+nregn,ntrend+nregn+nseas)
    Fseas = np.matlib.repmat([[1],[0]],pseas,1);Gseas = np.zeros([nseas,nseas]);

    #print (ntrend, nregn, nseas)
    # build the harmonic compenent of the Gseasonal
    for j in range(pseas):
        c = np.cos(2*np.pi*np.power(2,rseas[j]-1)/period);  #np.power(2,rseas[j]-1)
        s = np.sin(2*np.pi*np.power(2,rseas[j]-1)/period);
        i = np.arange(2*j,2*(j+1))
        Gseas[np.reshape(i,[2,1]),i] = [[c,s],[-s,c]]
    F = np.concatenate((Ftrend,Fregn,Fseas),axis=0)
    G = scipy.linalg.block_diag(Gtrend,Gregn,Gseas) 
    m = Prior.m; C = Prior.C; S = Prior.S; nu = Prior.nu
    
    #  
    T = len(Y)
    sm = np.zeros(m.shape)
    sC = np.zeros([C.shape[0],C.shape[1],1])
    sS = np.zeros(1)
    snu = np.zeros(1)
    slik = np.zeros(1)
    for t in range(T):
        a = np.dot(G,m)
        R = np.dot(np.dot(G,C),np.transpose(G))
        R[np.reshape(itrend,[-1,1]),itrend] = R[np.reshape(itrend,[-1,1]),itrend]/deltrend
        R[np.reshape(iregn,[-1,1]),iregn] = R[np.reshape(iregn,[-1,1]),iregn]/delregn
        R[np.reshape(iseas,[-1,1]),iseas] = R[np.reshape(iseas,[-1,1]),iseas]/delseas
        nu = delvar*nu
        F[iregn,0] = X[t,]

        A = np.dot(R,F);Q = np.squeeze(np.dot(np.transpose(F),A)+S); A = A/Q; f = np.squeeze(np.dot(np.transpose(F),a))
        y = Y[t]
        
        if ~np.isnan(y):
            e = y-f; ac = (nu+np.power(e,2)/Q)/(nu+1)
            rQ = np.sqrt(Q)
            mlik = tdstr.pdf(e/rQ,nu)/rQ
            m = a+A*e; C = ac*(R-np.dot(A,np.transpose(A))*Q); nu = nu+1; S = ac*S;
        else:
            m = a; C = R;
            mlik = np.nan
        sm = np.concatenate((sm,m),axis=1)
        sC = np.concatenate((sC,np.reshape(C,[C.shape[0],C.shape[1],1])),axis=2)
        snu = np.concatenate((snu,[nu]),axis=0)
        sS = np.concatenate((sS,[S]),axis=0)
        slik = np.concatenate((slik,[mlik]),axis=0)            
    return {'sm':sm, 'sC':sC ,'snu':snu,'slik':slik}

def de_seasonalize(vi):
    n_years = int(len(vi)/12)
    if n_years == len(vi)/12.0:
        monthly_vi = vi[0:(n_years*12)].reshape([n_years,12])
        mean_seasonal_pattern = np.nanmean(monthly_vi,axis=0)
        deseasonalized = vi - np.tile(mean_seasonal_pattern, n_years)
    else:
        monthly_vi = vi[6:(n_years*12+6)].reshape([n_years,12])
        mean_seasonal_pattern = np.nanmean(monthly_vi,axis=0)
        #print(mean_seasonal_pattern)
        deseasonalized = vi - np.concatenate((mean_seasonal_pattern[6:12],np.tile(mean_seasonal_pattern, n_years)),axis=None)
    #print(len(deseasonalized))
    return deseasonalized #signal.detrend(deseasonalized)

def deseasonalize_dlm(ndvi,spin_up_years,spin_up_rep,delta,rseas):
    # use dlm to deseasonalize
    ndvi_obs = len(ndvi)
    if (sum(np.isnan(ndvi))>300):
        res = np.ones((ndvi_obs,))
        res[:] = np.nan
        return res

    Y_raw = ndvi-np.nanmean(ndvi) 

    if spin_up_years>0:
        Y = np.concatenate((np.tile(Y_raw[0:(spin_up_years*12)],spin_up_rep),Y_raw),axis=0)
    else:
        Y = Y_raw
    X = np.zeros((len(Y),1))
    M = Model(Y,X,rseas,delta)
    FF = forwardFilteringM(M)
    # trend component
    sm = (FF.get('sm'))[:,(spin_up_rep*spin_up_years*12+1):] # 1d
    # seasonal component
    seas = np.nansum(sm[np.array(rseas)*2+1,:],axis=0)
    deseasonalized = Y[(spin_up_rep*spin_up_years*12):] - sm[0,:]-seas
    deseasonalized = np.where(np.isnan(deseasonalized),0,deseasonalized)
    return deseasonalized


def calprecip(pre):
    ## get the precipitation for current to previous 2 
    n_num = len(pre)    # num of observations, should be from 1979/07/01 to 2015/12/1  total of 438
    ## precipitation should be deseasonalized
    pre12 = de_seasonalize(np.convolve(pre[23:-1], np.ones((2,)), mode='valid'))
    #pre1 = pre[24:]
    X = pre12
    #print(X.shape)
    return X


def perpixel_vi(vec,spin_up_years,spin_up_rep,delta,rseas):
    # the vec is a combination of ndvi and precipitation
    # length of ndvi is 414, length of prec is 438
    nreg=2
    ndvi_obs = int((len(vec)-24)/nreg)
    ndvi = vec[0:ndvi_obs]
    ndvi[ndvi<0] = np.nan
    prec = vec[ndvi_obs:]
    #previou one and current months
    #pre12 = de_seasonalize(np.convolve(prec[23:], np.ones((2,)), mode='valid'))   ##start from 0 
    ## previous one months
    pre12 = de_seasonalize(prec[23:-1])
    # use three seasonal harmonic components
    n_var=2+nreg+len(rseas)*2
    ###### use the deseasonalized and detrended VI data,
    if (sum(np.isnan(ndvi))>300):
        res = np.ones(((ndvi_obs+spin_up_years*spin_up_rep*12)*(n_var*2+2+nreg+1),)) #np.ones((ndvi_obs+n_var*ndvi_obs+n_var*ndvi_obs+ndvi_obs+spin_up_years*spin_up_rep*12*(n_var*2+2+6+1)+ndvi_obs*7,))
        res[:] = np.nan
        return res

    Y_raw = ndvi-np.nanmean(ndvi)     # ndvi[1:]
    deseason_vi = deseasonalize_dlm(ndvi,spin_up_years,spin_up_rep,delta,rseas)
    
    #pre_vi = calvi(deseason_vi)
    pre_vi = np.concatenate((deseason_vi[:12],deseason_vi),axis=0)
    #print(pre_vi.shape)
    X_raw = np.column_stack((pre_vi[11:-1],pre12))    #pre_vi[:-1,],pre13[:-1]
    #print(X_raw)
    #print(X_raw[:,0]-deseason_vi[:-1])
    if spin_up_years>0:
        X = np.concatenate((np.tile(X_raw[0:(spin_up_years*12),:],(spin_up_rep,1)),X_raw),axis=0)
        Y = np.concatenate((np.tile(Y_raw[0:(spin_up_years*12)],spin_up_rep),Y_raw),axis=0)
    else:
        X = X_raw
        Y = Y_raw
    #print (Y.shape,X.shape)
    M = Model(Y,X,rseas,delta)
    FF = forwardFilteringM(M)

    Xout = np.transpose(X)
    Yout = Y
    slik = FF.get('slik')[1:] # 1d
    # extract estimates on the coefficient corresponding to lag-1 NDVI
    #vid = 2 # index of autocorrelation
    sm = FF.get('sm')[:,1:] # mean of autocorrelation # 2d  
    sC = FF.get('sC')[:,:,1:] # variance of autocorrelation  # 3d
    sCdiag = np.transpose(np.diagonal(sC, axis1=0,axis2=1).reshape(ndvi_obs+spin_up_years*spin_up_rep*12,n_var))
    snu = FF.get('snu')[1:] # degree of freedom   #1d
    #print("snu",len(snu))
    # vectorize
    res = np.concatenate((slik,sm,sCdiag,snu,Xout,Yout), axis=None)
    #print(res.shape)
    return res


def perpixel_vi_multi(vec,spin_up_years,spin_up_rep,delta,rseas):
    # the vec is a combination of ndvi and precipitation
    # length of ndvi is 414, length of prec is 438
    nreg=4
    ndvi_obs = int((len(vec)-24)/nreg)
    ndvi = vec[0:ndvi_obs]
    ndvi[ndvi<0] = np.nan
    clou = vec[ndvi_obs:(ndvi_obs*2)]
    temp = vec[(ndvi_obs*2):(ndvi_obs*3)]
    prec = vec[(ndvi_obs*3):(ndvi_obs*4+24)]

    #previou one and current months
    #pre12 = de_seasonalize(np.convolve(prec[23:], np.ones((2,)), mode='valid'))   ##start from 0 
    ## previous one months
    pre12 = de_seasonalize(prec[23:-1])
    clouano = de_seasonalize(clou)
    tempano = de_seasonalize(temp)
    # use three seasonal harmonic components
    n_var=2+nreg+len(rseas)*2
    ###### use the deseasonalized VI data,
    if (sum(np.isnan(ndvi))>300):
        res = np.ones(((ndvi_obs+spin_up_years*spin_up_rep*12)*(n_var*2+2+nreg+1),)) #np.ones((ndvi_obs+n_var*ndvi_obs+n_var*ndvi_obs+ndvi_obs+spin_up_years*spin_up_rep*12*(n_var*2+2+6+1)+ndvi_obs*7,))
        res[:] = np.nan
        return res

    Y_raw = ndvi-np.nanmean(ndvi)     # ndvi[1:]
    deseason_vi = deseasonalize_dlm(ndvi,spin_up_years,spin_up_rep,delta,rseas)
    
    #pre_vi = calvi(deseason_vi)
    pre_vi = np.concatenate((deseason_vi[:12],deseason_vi),axis=0)
    #print(pre_vi.shape)
    X_raw = np.column_stack((pre_vi[11:-1],pre12,clouano,tempano))    #pre_vi[:-1,],pre13[:-1]
    #print(X_raw)
    #print(X_raw[:,0]-deseason_vi[:-1])
    if spin_up_years>0:
        X = np.concatenate((np.tile(X_raw[0:(spin_up_years*12),:],(spin_up_rep,1)),X_raw),axis=0)
        Y = np.concatenate((np.tile(Y_raw[0:(spin_up_years*12)],spin_up_rep),Y_raw),axis=0)
    else:
        X = X_raw
        Y = Y_raw
    #print (Y.shape,X.shape)
    M = Model(Y,X,rseas,delta)
    FF = forwardFilteringM(M)

    Xout = np.transpose(X)
    Yout = Y
    slik = FF.get('slik')[1:] # 1d
    # extract estimates on the coefficient corresponding to lag-1 NDVI
    #vid = 2 # index of autocorrelation
    sm = FF.get('sm')[:,1:] # mean of autocorrelation # 2d  
    sC = FF.get('sC')[:,:,1:] # variance of autocorrelation  # 3d
    sCdiag = np.transpose(np.diagonal(sC, axis1=0,axis2=1).reshape(ndvi_obs+spin_up_years*spin_up_rep*12,n_var))
    snu = FF.get('snu')[1:] # degree of freedom   #1d
    #print("snu",len(snu))
    # vectorize
    res = np.concatenate((slik,sm,sCdiag,snu,Xout,Yout), axis=None)
    #print(res.shape)
    return res


def de_seasonalize_wtrend(vi):
    n_years = int(len(vi)/12)
    if n_years == len(vi)/12.0:
        monthly_vi = vi[0:(n_years*12)].reshape([n_years,12])
        mean_seasonal_pattern = np.nanmean(monthly_vi,axis=0)
        deseasonalized = vi - np.tile(mean_seasonal_pattern, n_years)
    else:
        monthly_vi = vi[6:(n_years*12+6)].reshape([n_years,12])
        mean_seasonal_pattern = np.nanmean(monthly_vi,axis=0)
        #print(mean_seasonal_pattern)
        deseasonalized = vi - np.concatenate((mean_seasonal_pattern[6:12],np.tile(mean_seasonal_pattern, n_years)),axis=None)
    #print(len(deseasonalized))
    return deseasonalized

def perpixel_vi_multi_wtrend(vec,spin_up_years,spin_up_rep,delta,rseas):
    # the vec is a combination of ndvi and precipitation
    # length of ndvi is 414, length of prec is 438
    nreg=4
    ndvi_obs = int((len(vec)-24)/nreg)
    ndvi = vec[0:ndvi_obs]
    ndvi[ndvi<0] = np.nan
    clou = vec[ndvi_obs:(ndvi_obs*2)]
    temp = vec[(ndvi_obs*2):(ndvi_obs*3)]
    prec = vec[(ndvi_obs*3):(ndvi_obs*4+24)]

    #previou one and current months
    #pre12 = de_seasonalize(np.convolve(prec[23:], np.ones((2,)), mode='valid'))   ##start from 0 
    ## previous one months
    pre12 = de_seasonalize_wtrend(prec[23:-1])
    clouano = de_seasonalize_wtrend(clou)
    tempano = de_seasonalize_wtrend(temp)
    # use three seasonal harmonic components
    n_var=2+nreg+len(rseas)*2
    ###### use the deseasonalized VI data,
    if (sum(np.isnan(ndvi))>300):
        res = np.ones(((ndvi_obs+spin_up_years*spin_up_rep*12)*(n_var*2+2+nreg+1),)) #np.ones((ndvi_obs+n_var*ndvi_obs+n_var*ndvi_obs+ndvi_obs+spin_up_years*spin_up_rep*12*(n_var*2+2+6+1)+ndvi_obs*7,))
        res[:] = np.nan
        return res

    Y_raw = ndvi-np.nanmean(ndvi)     # ndvi[1:]
    deseason_vi = de_seasonalize_wtrend(ndvi)
    
    #pre_vi = calvi(deseason_vi)
    pre_vi = np.concatenate((deseason_vi[:12],deseason_vi),axis=0)
    #print(pre_vi.shape)
    X_raw = np.column_stack((pre_vi[11:-1],pre12,clouano,tempano))    #pre_vi[:-1,],pre13[:-1]
    #print(X_raw)
    #print(X_raw[:,0]-deseason_vi[:-1])
    if spin_up_years>0:
        X = np.concatenate((np.tile(X_raw[0:(spin_up_years*12),:],(spin_up_rep,1)),X_raw),axis=0)
        Y = np.concatenate((np.tile(Y_raw[0:(spin_up_years*12)],spin_up_rep),Y_raw),axis=0)
    else:
        X = X_raw
        Y = Y_raw
    #print (Y.shape,X.shape)
    M = Model(Y,X,rseas,delta)
    FF = forwardFilteringM(M)

    Xout = np.transpose(X)
    Yout = Y
    slik = FF.get('slik')[1:] # 1d
    # extract estimates on the coefficient corresponding to lag-1 NDVI
    #vid = 2 # index of autocorrelation
    sm = FF.get('sm')[:,1:] # mean of autocorrelation # 2d  
    sC = FF.get('sC')[:,:,1:] # variance of autocorrelation  # 3d
    sCdiag = np.transpose(np.diagonal(sC, axis1=0,axis2=1).reshape(ndvi_obs+spin_up_years*spin_up_rep*12,n_var))
    snu = FF.get('snu')[1:] # degree of freedom   #1d
    #print("snu",len(snu))
    # vectorize
    res = np.concatenate((slik,sm,sCdiag,snu,Xout,Yout), axis=None)
    #print(res.shape)
    return res

