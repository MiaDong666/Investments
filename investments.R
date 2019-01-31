####strategy
library("lubridate")
library("magrittr")
library("tidyverse")
library("xts")
###initial data
filepath = "feds200628.csv"
Fed_ZCBs_raw=read.csv(filepath, header = TRUE, sep=",",skip=9)
row.names(Fed_ZCBs_raw)=Fed_ZCBs_raw$X
Fed_ZCBs=as.xts(Fed_ZCBs_raw)
Fed_ZCBs=Fed_ZCBs["1983-12-30/2018-06-30",c("SVENY02","SVENY10","BETA0","BETA1","BETA2","BETA3","TAU1","TAU2")]
storage.mode(Fed_ZCBs) <- "double"
weekly <- endpoints(Fed_ZCBs, on="weeks",k=1)
Fed_ZCBs <- period.apply(Fed_ZCBs, INDEX = weekly, FUN = first)
####pre-functions
#bond price fuction
Price=function(y,t){
  P=100*exp(-y/100*t)
  return(P)
}
#the NSS yield curve
NSS=function(n){
  beta1_weight=(1-exp(-n/Fed_ZCBs$TAU1))/(n/Fed_ZCBs$TAU1)
  beta2_weight=(beta1_weight-exp(-n/Fed_ZCBs$TAU1))
  beta3_weight=(1-exp(-n/Fed_ZCBs$TAU2))/(n/Fed_ZCBs$TAU2)-exp(-n/Fed_ZCBs$TAU2)
  yield=Fed_ZCBs$BETA0+beta1_weight*Fed_ZCBs$BETA1+beta2_weight*Fed_ZCBs$BETA2+beta3_weight*Fed_ZCBs$BETA3
  return(yield)
}

#bond price and DV01 of 2-year and 10-year
Fed_ZCBs$Price_2yr=Price(Fed_ZCBs$SVENY02,2)
Fed_ZCBs$DV01_2yr=Fed_ZCBs$Price_2yr*2*1/10000
Fed_ZCBs$Price_10yr=Price(Fed_ZCBs$SVENY10,10)
Fed_ZCBs$DV01_10yr=Fed_ZCBs$Price_10yr*10*1/10000

#solving for x 10-year bond of pl=DV01_2yr/DV01_10yr
Fed_ZCBs$x=Fed_ZCBs$DV01_2yr/Fed_ZCBs$DV01_10yr

#yield and price of the bonds of last week
t_2=1+51/52
Fed_ZCBs$yield_2yrlw=NSS(t_2)
Fed_ZCBs$Price_2yrlw=Price(Fed_ZCBs$yield_2yrlw,t_2)
t_10=9+51/52
Fed_ZCBs$yield_10yrlw=NSS(t_10)
Fed_ZCBs$Price_10yrlw=Price(Fed_ZCBs$yield_10yrlw,t_10)

###DATA port concerning price and balance ratio
port=Fed_ZCBs[,c("Price_2yr","Price_10yr","Price_2yrlw","Price_10yrlw","x")]
N=length(port$Price_2yr)
###Every week we short a unit of 2-year bond and long ax unit 10-year bond 
#that the margin maches our cash level C.
port$a=1:N
port$Price_sec=1:N
port$margin=1:N
port$Capos=1:N
port$Cash_closed=1:N
port$r=NSS(1/52)
port$Interest=1:N
port$timepassage=1:N


###week-start function
week_start=function(c_0,x,p2,p10,ratio){
  
  portfolio=list(a=0)
  portfolio$a=(c_0/ratio)/(p2+x*p10)
  portfolio$margin=portfolio$a*(p2+x*p10)*ratio
  portfolio$Price_sec=-portfolio$a*p2+portfolio$a*x*p10
  portfolio$Capos=-portfolio$Price_sec+c_0
  return(portfolio)

}




#start at C0=1
###solve
# (a0*pr2+a0*x*pr10)*0.1=c0
c_0=1
p_0=week_start(c_0,x=port[1,"x"], p2=port[1,"Price_2yr"],p10=port[1,"Price_10yr"],ratio=0.1)
port[1,"a"]=p_0$a
port[1,"margin"]=p_0$margin
port[1,"Price_sec"]=p_0$Price_sec
port[1,"Capos"]=p_0$Capos
port[1,"Cash_closed"]=1
port[1,"Interest"]=port[1,"Capos"]*(exp(1/52*port[1,"r"]/100)-1)


###week-end function
###cashend=xo*a_old*p10-a_old*p2+c0+I
#rebalance
for(i in 2:N){
  
  c_old=as.double(port[i-1,"Capos"])
  c_oldclos=as.double(port[i-1,"Cash_closed"])
  I=as.double(port[i-1,"Interest"])
  a_old=as.double(port[i-1,"a"])
  x_old=as.double(port[i-1,"x"])
  port[i,"Cash_closed"]=c_old+I+a_old*(x_old*port[i,"Price_10yrlw"]-port[i,"Price_2yrlw"])
  p=week_start(c_0=port[i,"Cash_closed"],
                  x=port[i,"x"],
                  p2=port[i,"Price_2yr"],
                  p10=port[i,"Price_10yr"],
                  ratio=0.1)
  port[i,"a"]=p$a
  port[i,"margin"]=p$margin
  port[i,"Price_sec"]=p$Price_sec
  port[i,"Capos"]=p$Capos
  port[i,"Interest"]=port[i,"Capos"]*(exp(1/52*port[i,"r"]/100)-1)

}



##############return part
#NEW: re-did columns for returns
port$timepassage=-port$a*(port$Price_2yrlw - port$Price_2yr) + port$a * port$x * (port$Price_10yrlw - port$Price_10yr)
port$Int_cum=cumsum(port$Interest)
port$tps_cum=cumsum(port$timepassage)
port$time_cum=port$Int_cum+port$tps_cum


###making the return sheet
###spreadreturn
return=cbind(port$time_cum,port$a,port$x,Fed_ZCBs$DV01_2yr,Fed_ZCBs$DV01_10yr,port$Price_2yr,port$Price_10yr)
return=return[-N,]
a=as.numeric(diff(Fed_ZCBs$SVENY02,differences = 1))
return$dy2=a[-1]
b=as.numeric(diff(Fed_ZCBs$SVENY10,differences = 1))
return$dy10=b[-1]
#spred=a*DV2*dy2-a*x*DV10*dy10
return$spred=(return$a*return$DV01_2yr*return$dy2-return$a*return$x*return$DV01_10yr*return$dy10)*100
return$spred_cum=cumsum(return$spred)
#convexity
return$cvx=0.5*100*(return$dy10/100)^2*return$a*return$x*return$Price_10yr-
  0.5*4*(return$dy2/100)^2*return$a*return$Price_2yr
return$cvx_cum=cumsum(return$cvx)
#gain
c=as.numeric(diff(port$Cash_closed,differences = 1))
return$totalreturn=c[-1]
return$sum_cum=cumsum(return$totalreturn)
return$residual=return$sum_cum-return$time_cum-return$spred_cum-return$cvx_cum


summary <- cbind(return$sum_cum,return$time_cum, return$spred_cum,return$cvx_cum,return$residual)
colnames(summary) <- c("Cum Return", "Time Value", "spread Return","Convexity return","Residual")
plot.xts(summary,
         col=c("black","red","blue","green","pink"),main="return", legend.loc = "bottomleft")
head(return)
summary <- cbind(return$sum_cum,return$time_cum, return$spred_cum,return$cvx_cum,return$residual)


