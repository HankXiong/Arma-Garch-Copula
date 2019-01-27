# Copula package
library(copula)
# Fancy 3D plain scatterplots
library(scatterplot3d)
# ggplot2
library(ggplot2)
# Useful package to set ggplot plots one next to the other
library(grid)
set.seed(235)
######
setwd('E:/Files/Bachelar paper/paper/paper_math/data and code/')
sp=read.csv('sp500.csv',header=T)
to=read.csv('toronto.csv',header=T)
sp$Date=as.Date(sp$Date);to$Date=as.Date(to$Date)
library(xts);library(PerformanceAnalytics)

sp=xts(sp[,-1],order.by = sp$Date);to=xts(to[,-1],order.by = to$Date)
sp$Ret=CalculateReturns(sp$Adj.Close,'log')
to$Ret=CalculateReturns(to$Adj.Close,'log')

mydata=merge.xts(to$Ret,sp$Ret,join='inner')[-1,]
colnames(mydata)=c('to_ret','sp_ret')

date1=index(sp)
date2=index(to)
length(intersect(date1,date2))
date=index(mydata)
l=dim(mydata)[1]

### plot returns
par(mfcol=c(2,2),mar=c(2,2,2,2))
plot(date,mydata$sp_ret,main='Log Returns of SP500',xlab='Date',ylab='Return',type='l',col=1,lwd=0.1)
plot(date,mydata$to_ret,main='Log Returns of TSE',xlab='Date',ylab='Return',type='l',col=1,lwd=0.1)
plot(date,abs(mydata$sp_ret),main='Absolute Log Returns of SP500',xlab='Date',ylab='Return',type='l',col=1,lwd=0.1)
plot(date,abs(mydata$to_ret),main='Absolute Log Returns of TSE',xlab='Date',ylab='Return',type='l',col=1,lwd=0.1)

library(urca)
## summarize some statistics
return_summary=function(x){
  mean=Return.annualized(x,scale=250)
  sd=sd.annualized(x,scale=250)
  min=min(x)
  median=median(x)
  max=max(x)
  skew=PerformanceAnalytics::skewness(x,method='sample')
  kur=PerformanceAnalytics::kurtosis(x,method='sample')
  
  print(Box.test(x,10,type='L'))
  print(Box.test(x^2,10,type='L'))
  print(ur.df(x,lags=1))
  return(c(mean,sd,min,median,max,skew,kur))
}
apply(mydata, 2, return_summary)


######################
library(fGarch)
library(forecast)
library(copula)

da1=coredata(mydata['2008-07-29/2009-08-01']$sp_ret)
da2=coredata(mydata['2008-07-29/2009-08-01']$to_ret)
summary(ur.df(da1,type='none',lags=10))
summary(ur.df(da2,type='none',lags=10))##no unitroot, it's stationary

########################### arma-garch for an interval
par(mfcol=c(1,2))
acf(da1);pacf(da1)

m1=auto.arima(da1,max.p=5,max.q=5,max.order = 10,stationary = T,seasonal = F,allowdrift = F,allowmean = F,ic='bic',stepwise = F,parallel = T)
m2=auto.arima(da2,max.order = 6,stationary = T,seasonal = F,allowdrift =F,allowmean = F,ic='bic',stepwise = F,parallel = T)

m1_garch=garchFit(~arma(0,0)+garch(1,1),da1,cond.dist = 'sged',trace=F,include.mean = F)
m2_garch=garchFit(~arma(0,0)+garch(1,1),da2,cond.dist = 'sged',trace=F,include.mean = F)
summary(m1_garch)
summary(m2_garch)

## plot sp500 return with tsx return and their standardized return
par(mfcol=c(1,2),mar=c(4,4,2,2))
plot(da1,da2,xlim=c(quantile(da1,0.00),quantile(da1,1)),
     ylim=c(quantile(da2,0),quantile(da2,1)),xlab='SP500 Return',ylab='TSE Return Series')
plot(residuals(m1_garch,T),residuals(m2_garch,T),xlab='Standardized Innovation-SP500',ylab='Standardized In?ovation-TSE')


##############can't reject that the residuals follows t distribution
ks.test(residuals(m1_garch),'psged',mean = mean(m1_garch@residuals),sd=sd(m1_garch@residuals),nu=coef(m1_garch)['shape'],xi=coef(m1_garch)['skew'])
hist(psged(residuals(m1_garch,T),nu=coef(m1_garch)['shape'],xi=coef(m1_garch)['skew']),20)  ## almost uniform distribution
ks.test(psged(residuals(m2_garch),mean = mean(m2_garch@residuals),sd=sd(m2_garch@residuals),nu=coef(m2_garch)['shape'],xi=coef(m2_garch)['skew']),punif,min=0,max=1)

##############fit copula with moving window and compute the their distance with empirical copula
upper = 0.05
U=cbind(rep(seq(0,upper,by=0.005),times=11),rep(seq(0,upper,by=0.005),each=11))##evalution area
dis=function(x,y){
  return(sqrt(mean((x-y)^2)))
}
distance=c()
win = 250
for (i in seq(win,l,5)){
  da1=mydata[(i-win):i,1];da2=mydata[(i-win):i,2]
  m1=auto.arima(da1,max.order = 6,stationary = T,seasonal = F,allowdrift = F,allowmean = F,ic='bic')
  m2=auto.arima(da2,max.order = 6,stationary = T,seasonal = F,allowdrift = F,allowmean = F,ic='bic')
  
  m1_garch=garchFit(~garch(1,1),m1$residuals,cond.dist = 'sged',trace=T,include.mean = F)
  a=psged(residuals(m1_garch,T),nu=coef(m1_garch)['shape'],xi=coef(m1_garch)['skew'])
  
  m2_garch=garchFit(~garch(1,1),m2$residuals,cond.dist = 'sged',trace=T,include.mean = F)
  b=psged(residuals(m2_garch,T),nu=coef(m2_garch)['shape'],xi=coef(m2_garch)['skew'])
  X=cbind(a,b)

  fitt=fitCopula(tCopula(dim=2,dispstr='un'),cbind(a,b),method='itau.mpl')
  fitnormal=fitCopula(normalCopula(dim=2,dispstr='un'),cbind(a,b),method='mpl')
  fitclay=fitCopula(claytonCopula(dim=2),cbind(a,b),method='mpl')
  fitgumbel_rot=fitCopula(rotCopula(gumbelCopula(dim=2)),cbind(a,b),method='mpl')
  fitgumbel=fitCopula(gumbelCopula(dim=2),cbind(a,b),method='mpl')

  ec=C.n(U,X,smoothing = 'none')
  pc_t=pCopula(U,tCopula(param =coef(fitt)['rho.1'],df=round(coef(fitt)['df']),dim = 2 ))
  pc_clayton=pCopula(U,claytonCopula(param=coef(fitclay)['alpha'],dim=2))
  pc_normal=pCopula(U,normalCopula(param = coef(fitnormal)['rho.1'],dim=2))
  pc_gum_rot=pCopula(U,rotCopula(gumbelCopula(param=coef(fitgumbel_rot)['alpha'],dim=2)))
  pc_gum=pCopula(U,gumbelCopula(param=coef(fitgumbel)['alpha'],dim=2))
  
  distance=rbind(distance,c(dis(ec,pc_normal),dis(ec,pc_t),dis(ec,pc_clayton),dis(ec,pc_gum_rot),dis(ec,pc_gum)))
}
apply(distance, 2, mean)## normal, t, clayton, gum_rot, gum

##compare on the diagnal
C_n_diag=function(u){C.n(do.call(cbind,rep(list(u),2)),X=X)}

cop_t=tCopula(param =coef(fitt)['rho.1'],df=floor(coef(fitt)['df']),dim = 2 )
C_t_diag=function(u){pCopula(do.call(cbind,rep(list(u),2)),cop_t)}

cop_normal=normalCopula(param = coef(fitnormal)['rho.1'],dim=2)
C_normal_diag=function(u){pCopula(do.call(cbind,rep(list(u),2)),cop_normal)}

cop_clayton=claytonCopula(param=coef(fitclay)['alpha'],dim=2)
C_clayton_diag=function(u){pCopula(do.call(cbind,rep(list(u),2)),cop_clayton)}

cop_gumbel_rot=rotCopula(gumbelCopula(param=coef(fitgumbel_rot)['alpha'],dim=2))
C_gumbel_diag=function(u){pCopula(do.call(cbind,rep(list(u),2)),cop_gumbel_rot)}

C_gumbel_diag_1=function(u){pCopula(do.call(cbind,rep(list(u),2)),gumbelCopula(param=coef(fitgumbel)['alpha'],dim=2))}

curve(C_n_diag,0.1,0.9,xlab = 'u',ylab = 'Probability')
curve(C_t_diag,col=2,add=T)
curve(C_gumbel_diag,col=5,add=T)
curve(C_normal_diag,col=3,add=T)
curve(C_clayton_diag,col=4,add=T)
curve(C_gumbel_diag_1,col=6,add=T)
legend('topleft',lty=1:1,col=1:6,legend = c('Empirical','Student t','Gaussian','Clayton','Rotated Gumbel','Gumbel'),ncol = 1)
######copula simulate
cal_var=function(x){
 return(cbind(quantile(x,0.1),quantile(x,0.05),quantile(x,0.01),quantile(x,0.005)))
}
win=250
VaR=var=c()
set.seed(1)
for (i in seq(win,l,5))
  {
  da1=mydata[(i-win):i,1];da2=mydata[(i-win):i,2]
  m1=auto.arima(da1,max.order = 6,stationary = T,seasonal = F,allowdrift = F,allowmean = F,ic='bic')
  m2=auto.arima(da2,max.order = 6,stationary = T,seasonal = F,allowdrift = F,allowmean = F,ic='bic')

  m1_garch=garchFit(~garch(1,1),m1$residuals,cond.dist = 'sged',trace=F,include.mean = F,control = list(ndeps=10^-4,reltol=10^-9,factr=10^-9))
  a=psged(residuals(m1_garch,T),nu=coef(m1_garch)['shape'],xi=coef(m1_garch)['skew'])
  
  m2_garch=garchFit(~garch(1,1),m2$residuals,cond.dist = 'sged',trace=F,include.mean = F,control = list(ndeps=10^-4,reltol=10^-9,factr=10^-9))
  b=psged(residuals(m2_garch,T),nu=coef(m2_garch)['shape'],xi=coef(m2_garch)['skew'])
  pair=cbind(a,b)
  
  #########################start fit
  #fit=fitCopula(tCopula(dim=2,dispstr='un'),pair,method='itau.mpl')
  #fit=fitCopula(gumbelCopula(dim=2),pair,method='mpl')
  fit=fitCopula(rotCopula(gumbelCopula(dim=2)),pair,method='mpl',optim.control=list(trace=F,ndeps=10^-4,factr=10^8))
  
  #sico=tCopula(dim=2,param = fit@estimate[1],df=fit@estimate[2])
  sico=rotCopula(gumbelCopula(dim=2,param = fit@estimate[1])) ##inverse

  ###################start simulate
  simulate=rCopula(400*1000,sico)
  #simulate=pair[sample(1:win,size=400*1000,replace=T),]   ##empirical copula
  epsi1=qsged(simulate[,1],nu=coef(m1_garch)['shape'],xi=coef(m1_garch)['skew'])  ##empirical marginal distribution
  epsi2=qsged(simulate[,2],nu=coef(m2_garch)['shape'],xi=coef(m2_garch)['skew']) 
  #epsi1=quantile(ecdf(residuals(m1_garch,T)),simulate[,1])
  #epsi2=quantile(ecdf(residuals(m2_garch,T)),simulate[,2])
  a1=as.numeric(predict(m1_garch,1))[3]*epsi1
  a2=as.numeric(predict(m2_garch,1))[3]*epsi2
  pre1=predict(m1,1,se=F)[1]+a1;pre2=predict(m2,1,se=F)[1]+a2
  
  pre=matrix((pre1+pre2)/2,400,1000)
  
  VaR=rbind(VaR,apply(apply(pre,2,cal_var),1,mean))
}
sdate=date[seq(win,l,5)]
VaR=xts(VaR,order.by = sdate)
names(VaR)=c('var0.1','var0.05','var0.01','var_0.005')
###
VaR_empi=VaR
VaR12=VaR
VaR_delta=VaR

VaR_empi_ged=VaR

#
VaR_clayton_ged=VaR
VaR_gumbel_ged=VaR
VaR_t_ged=VaR
VaR_normal_ged=VaR

VaR_gumbel_inverse_ged=VaR
### empirical VaR
VaR=NULL
for (i in seq(win,l,5)) {
  da1=mydata[(i-win):i,1]
  da2=mydata[(i-win):i,2]
  da=coredata((da1+da2)/2)
  
  da_sample=matrix(da[sample(1:win,400*1000,replace=T)],400,1000)

  VaR=rbind(VaR,apply(apply(da_sample, 2, cal_var),1,mean))
}

####var1+var2
VaR=var1=var2=NULL
for (i in seq(win,l,5)) {
  da1=coredata(mydata[(i-win):i,1])
  da2=coredata(mydata[(i-win):i,2])
  da1_sample=matrix(da1[sample(1:win,size=1000*2000,replace = T)],1000,2000)
  da2_sample=matrix(da2[sample(1:win,size=1000*2000,replace = T)],1000,2000)
  var1=apply(apply(da1_sample, 2, cal_var),1,mean)
  var2=apply(apply(da2_sample, 2, cal_var),1,mean)
  VaR=rbind(VaR,(var1+var2)/2)
}
##delta VaR
VaR=NULL
for (i in seq(win,l,5)) {
  da1=coredata(mydata[(i-win):i,1])
  da2=coredata(mydata[(i-win):i,2])
  rho=as.numeric(cor(da1,da2,method='spearman'))
  std=as.numeric(sqrt(var(da1)/4+var(da2)/4+cov(da1,da2)/2))
  u=(mean(da1)+mean(da2))/2
  VaR=rbind(VaR,c(u+qnorm(0.1)*std,u+qnorm(0.05)*std,u+qnorm(0.01)*std,u+qnorm(0.005)*std))
}
names(VaR)=c('var0.1','var0.05','var0.01','var0.005')

##check significance
n=562
mydata$port_ret=(mydata[,1]+mydata[,2])/2
port_ret_test=mydata[seq(win,l-5,5)+1,3]
prob=c(0.1,0.05,0.01,0.005)
for( i in 1:4){
  var=VaR_clayton_ged[-dim(VaR12)[1],i];p=prob[i]
  m=sum(coredata(var)>coredata(port_ret_test))
  statistic=-2*log((1-p)^(n-m)*p^m)+2*log((1-m/n)^(n-m)*(m/n)^m)
  print(m)
  print(paste(round(pchisq(statistic,1,lower.tail = F)*100,3),'%'))
}

#
pchisq(-2*log((1-p)^(n-m)*p^m)+2*log((1-m/n)^(n-m)*(m/n)^m),1,lower.tail =F)

##plot
par(lty=3,lwd=0.5)
plot(sdate,VaR_empi$var0.01,col=1,type='l',main='Estimated VaR_0.1',xlab='Year',ylab='VaR')
lines(sdate,VaR_delta$var0.01,col=2,lwd=0.5)
lines(sdate,VaR12$var0.01,col=3,lwd=0.5)
legend('bottomright',legend=c('VaR_Empi','VaR_De?ta','VaR_Empi_Sum'),col=1:3,lty=3,box.lwd=0.5)

##
par(lty=1,lwd=0.5)
se=seq(1,length(sdate),5)
plot(sdate[se],VaR_clayton_ged$var0.01[se],col=1,type='l',main='Estimated VaR_0.1',xlab='Year',ylab='VaR',cex=0.1)
lines(sdate[se],VaR_gumbel_ged$var0.01[se],col=2,lwd=0.5)
lines(sdate[se],VaR_normal_ged$var0.01[se],col=3,lwd=0.5)
lines(sdate[se],VaR_t_ged$var0.01[se],col=4,lwd=0.5)
lines(sdate[se],VaR_empi_ged$var0.01[se],col=5,lwd=0.5)
legend('bottomright',legend = c('VaR_Clayton_Ged','VaR_Gumbel_Ged','VaR_Normal_Ged','VaR_t_Ged','VaR_Empi_Ged'),
       col=c(1?2,3,4,5),lty=1,ncol=2)

tdate=date[250:length(date)]
tdate=tdate[seq(2,length(tdate),5)]
plot(tdate,VaR_clayton_ged$var0.01[-dim(VaR)[1]],col=1,type='l',xlab='Year',ylab='VaR',ylim=c(-0.12,-0.01))
lines(tdate,VaR_clayton_ged$var0.05[-dim(VaR)[1]],col='blue')
points(tdate,mydata$port_ret[seq(251,dim(mydata)[1]-1,5)],pch=1,cex=0.5)

var3=VaR_clayton_ged[-dim(VaR12)[1],3]
breakdate3=tdate[coredata(var3)>coredata(port_ret_test)]
var2=VaR_clayton_ged[-dim(VaR12)[1],2]
breakdate2=tdate[coredata(var2)>coredata(port_ret_test)]
points(breakdate2,mydata$port_ret[breakdate2],col=2,pch=1,cex=1)
points(breakdate3,mydata$port_ret[breakdate3],col=2,pch=2,cex=1)
