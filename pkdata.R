#x0=140;thalf=115;tmax=4;Vz=319  示例
#k=log(2)/thalf
#f=function(ka){k*exp(-k*tmax)-ka*exp(-ka*tmax)}
#root<-uniroot(f,c(1,4),tol = 0.01)
#ka=root$root
x0=60
info.tmax=c(4.00,2.00,5.04)
info.thalf=c(111,19.6)
info.Vz=319
n=100
pargen<-function(info.tmax,info.thalf,info.Vz){
  tmax=runif(1,info.tmax[2],info.tmax[3])
  thalf=rnorm(1,info.thalf[1],info.thalf[2])
  Vz=info.Vz
  k=log(2)/thalf
  f=function(ka){k*exp(-k*tmax)-ka*exp(-ka*tmax)}
  root<-uniroot(f,c(0.05,10),tol = 0.01)
  ka=root$root
  par=cbind(x0,k,ka,Vz)
  par=data.frame(par)
  return(par)
}
par=pargen(info.tmax,info.thalf,info.Vz)
concentration<-function(pargen){#此为口服一房室模型
  x0=pargen$x0;k=pargen$k;ka=pargen$ka;Vz=pargen$Vz
  t<-c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,14,24,48,72,96,120,144,192,240,288,360)
  c=c(rep(0,27))
  I=ka*x0/(Vz*(ka-k))
  for (i in 1:length(t)){
    c[i]<-I*(exp(-k*t[i])-exp(-ka*t[i]))*1000
  }
  result<-cbind(t,c)
  result<-data.frame(result)
  return(result)
}
conc=concentration(pargen(info.tmax,info.thalf,info.Vz))
#plot(conc,type=)
#plot(conc$t,log(conc$c))

AUC<-function(conc){#求AUC
  AUC=c(rep(0,length(conc$t)-1))
  for (i in 1:length(conc$t)-1){
    AUC[i]=(conc$t[i+1]-conc$t[i])*(conc$c[i+1]+conc$c[i])/2
  }
  sum=sum(AUC)
  return(sum)
}
AUCsum<-AUC(conc)

simdata<-function(n,info.tmax,info.thalf,info.Vz){
  data=replicate(n,concentration(pargen(info.tmax,info.thalf,info.Vz)))
}
data<-simdata(n,info.tmax,info.thalf,info.Vz)
rm(result)
for (i in 1:12){
  ID=i
  case.ID="SX-Cabozantinib-101PK"
  info.ID=paste("K00",i,sep = "")
  Period=1
  if (i==1){
    result.R=rbind(cbind(ID,case.ID,info.ID,Period,data[[2*i-1]],data[[2*i]]))
  }else{
    result.R=rbind(result.R,cbind(ID,case.ID,info.ID,Period,data[[2*i-1]],data[[2*i]]))
  }
}
for (i in 31:42){
  ID=i-30
  case.ID="SX-Cabozantinib-101PK"
  info.ID=paste("K00",i-30,sep = "")
  Period=2
  data[[2*i]]=data[[2*i]]*rnorm(1,0.98,0.02)
  if (i==31){
    result.T=rbind(cbind(ID,case.ID,info.ID,Period,data[[2*i-1]],data[[2*i]]))
  }else{
    result.T=rbind(result.T,cbind(ID,case.ID,info.ID,Period,data[[2*i-1]],data[[2*i]]))
  }
}
resulttemp<-rbind(result.R,result.T)
colnames(resulttemp)<-c("ID","Study ID","Subject","Period","Time","Concentration")
result<-data.frame(resulttemp)
result$ID=as.numeric(result$ID)
is.numeric(result$ID)
result$Concentration[result$Concentration==0]="BQL"
result<-result[order(result$ID),]
result<-result[,-1]
write.csv(result,"F:\\result.csv")