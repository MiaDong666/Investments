
library(quantmod)
library(xts)
library(quadprog)
library(ggplot2)
library(readxl)

fileroute="lecture6p.xlsx"
rf<-read_excel(fileroute,6,skip=4)

###################################################################
#Q1
#get data and reshape
getSymbols("INTC",src = "yahoo", from = as.Date("1989-12-29"), to = as.Date("2018-09-28"))
getSymbols("MSFT",src = "yahoo", from = as.Date("1989-12-29"), to = as.Date("2018-09-28"))
getSymbols("LUV",src = "yahoo", from = as.Date("1989-12-29"), to = as.Date("2018-09-28"))
getSymbols("MCD",src = "yahoo", from = as.Date("1989-12-29"), to = as.Date("2018-09-28"))
getSymbols("JNJ",src = "yahoo", from = as.Date("1989-12-29"), to = as.Date("2018-09-28"))

raw_mix=cbind(INTC$INTC.Close,MSFT$MSFT.Close,LUV$LUV.Close,MCD$MCD.Close,JNJ$JNJ.Close)

#get the weekly return for five stocks
INTC_w <- weeklyReturn(INTC$INTC.Adjusted)
MSFT_w <- weeklyReturn(MSFT$MSFT.Adjusted)
LUV_w <- weeklyReturn(LUV$LUV.Adjusted)
MCD_w <- weeklyReturn(MCD$MCD.Adjusted)
JNJ_w <- weeklyReturn(JNJ$JNJ.Adjusted)

mix <- cbind(INTC_w,MSFT_w,LUV_w,MCD_w,JNJ_w)
mix=mix[-1,]
colnames(mix)=c("INTC","MSFT","LUV","MCD","JNJ")

#create a function to get minimum variance
solution=function(mix,etarget){
  
  ret=sapply(mix, mean)
  ret_a=ret*52
  var=sapply(mix,var)
  sd_a=sqrt(52)*sqrt(var)
  
  #calculate correlation matrix
  k=ncol(mix)
  n=nrow(mix)
  M=matrix(mix,ncol=k)
  mean=matrix(ret,ncol=k)
  M_mean=matrix(data=1,nrow=n)%*%mean
  D=M-M_mean
  C=1/(n-1)*t(D)%*%D
  corr=cov2cor(C)
  
  #create mean-variance frontier
  err=as.vector(ret_a)
  sigmarr=as.vector(sd_a)
  vrr=diag(sigmarr)%*%corr%*%diag(sigmarr)
  vrr <- round(vrr,6)
  onen=rep(1,k)
  A=cbind(onen,err)
  
  mrf1=solve.QP(2*vrr,rep(0,k),A,c(1,etarget),meq=2)
  
  return(mrf1)
}


etarget=seq(from=0.1, to=0.3,by= 0.0001)
sd_p <- 1:length(etarget)
for (i in 1:length(etarget)){
  mrf1=solution(mix=mix,etarget = etarget[i])
  sd_p[i]=sqrt(mrf1$value)
}


######################################################add the original frontier
mix_2=mix[,c(1:2)]

etarget2=seq(from=0.1, to=0.3,by= 0.0001)
sd_p2 <- 1:length(etarget2)
for (i in 1:length(etarget2)){
  mrf1 =solution(mix_2,etarget2[i])
  sd_p2[i]=sqrt(mrf1$value)
}

#getting two frontiers
ftp=as.data.frame(cbind(etarget,sd_p))
ft2=as.data.frame(cbind(etarget2,sd_p2))

plot(x=sd_p,y=etarget,pch=20,cex=0.1,xlab = "standard deviation",ylab = "return",xlim=c(0,0.5),ylim = c(0,0.3))

points(x=sd_p2,y=etarget2,pch=20,cex=0.1)

#find the minimum variance portfolio
ftp_min=which.min(ftp$sd_p)
p_min=ftp[ftp_min,]
ft2_min=which.min(ft2$sd_p2)
o_min=ft2[ft2_min,]

#draw the efficient frontier
eftp=ftp[ftp$etarget>=p_min$etarget,]
points(x=eftp$sd_p,y=eftp$etarget,pch=20,cex=0.1,col="purple")
eft2=ft2[ft2$etarget2>=o_min$etarget2,]
points(x=eft2$sd_p2,y=eft2$etarget2,pch=20,cex=0.1,col="blue")


#################################################################################################
#Q2
#find tangent port folio
rf=rf[rf$X__1>=19891229 & rf$X__1<=20180928,]
Rf=mean(rf$RF)*252/100
t_p1=which.max((eftp$etarget-Rf)/eftp$sd_p)
tan_p=eftp[t_p1,]
sharpe01=(eftp[t_p1,1]-Rf)/eftp[t_p1,2]
abline(a=Rf,b=sharpe01,col="purple")

t_p2=which.max((eft2$etarget2-Rf)/eft2$sd_p2)
tan_p2=eft2[t_p2,]
sharpe02=(eft2[t_p2,1]-Rf)/eft2[t_p2,2]
abline(a=Rf,b=sharpe02,col="blue")

###################################################################################################
#Q3 optimal assets

port=solution(mix=mix,etarget = tan_p$etarget)
A=5
w_p <- (tan_p$etarget-Rf)/(A*(tan_p$sd_p)^2)
w1=port$solution

w_t=matrix(w1*w_p,nrow=1)
w_t=cbind(w_t,1-w_p)
colnames(w_t)=c("INTC","MSFT","LUV","MCD","JNJ","risk_free")

ret_op=w_p*tan_p$etarget+(1-w_p)*Rf
sd_op=w_p*tan_p$sd_p
points(x=sd_op,y=ret_op,pch=20,col="red",cex=1)

text(sd_op,ret_op,labels = "Optimal Portfolio",col="red",pos=3)

ret=sapply(mix, mean)
ret_a=ret*52
var=sapply(mix,var)
sd_a=sqrt(52)*sqrt(var)
stock=cbind(sd_a,ret_a)

###adding stock potints
points(stock[1,1],stock[1,2],col="brown",pch=20)
text(stock[1,1],stock[1,2],labels = "INTC",col="brown",pos=3,pch=20)

points(stock[2,1],stock[2,2],col="brown",pch=20)
text(stock[2,1],stock[2,2],labels = "MSFT",col="brown",pos=3)

points(stock[3,1],stock[3,2],col="brown",pch=20)
text(stock[3,1],stock[3,2],labels = "LUV",col="brown",pos=3)

points(stock[3,1],stock[3,2],col="brown",pch=20)
text(stock[3,1],stock[3,2],labels = "LUV",col="brown",pos=3)

points(stock[4,1],stock[4,2],col="brown",pch=20)
text(stock[4,1],stock[4,2],labels = "MCD",col="brown",pos=3)

points(stock[5,1],stock[5,2],col="brown",pch=20)
text(stock[5,1],stock[5,2],labels = "JNJ",col="brown",pos=3)


##adding minimum and tangent portfolio
points(ftp[ftp_min,2],ftp[ftp_min,1],col="purple",cex=1,pch=20)
text(ftp[ftp_min,2],ftp[ftp_min,1],col="purple",label="min portfolio",pos=3)

points(ft2[ft2_min,2],ft2[ft2_min,1],col="blue",cex=1,pch=20)
text(ft2[ft2_min,2],ft2[ft2_min,1],col="blue",label="min portfolio",pos=3)


     