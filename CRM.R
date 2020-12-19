#CRM
sim.crm<-function(s,t,average,sigma,phi,n){
  L<-1
  MTD<-0
  x<-matrix(rep(0,3*n,ncol=1))
  r<-matrix(rep(0,n,ncol=1))
  d<-matrix(rep(0,n,ncol=1))
  m<-matrix(rep(0,length(s),ncol=1))
  d[1]=1
  EP<-matrix(rep(0,length(s)))
  sum.dlt<-matrix(rep(0,length(s)))
  sum.trial<-matrix(rep(0,length(s)))#�����ֵ
  L<-function(a){
    l<-1
    for (j in 1:length(s)){
      l=l*(s[j]^exp(a))^(sum.dlt[j])*(1-(s[j])^exp(a))^(sum.trial[j]-sum.dlt[j])
    }
    return(l)
  }#���������������µ���Ȼ����L(alpha|D)
  LF<-function(a)L(a)*1/(sqrt(2*pi)*sigma)*exp(-(a-average)^2/(2*sigma^2))#����L(alpha|D)f(alpha)
  ep<-function(a,j)LF(a)/int*s[j]^(exp(a))#����alpha���º�Ķ��Ը���
  for (i in 1:n){
    x[3*i-2]<-rbinom(1,1,t[d[i]])>0
    x[3*i-1]<-rbinom(1,1,t[d[i]])>0
    x[3*i]<-rbinom(1,1,t[d[i]])>0#����һ�黼������״��
    r[i]<-x[3*i-2]+x[3*i-1]+x[3*i]#DLT����
    sum.dlt[d[i]]<-sum.dlt[d[i]]+r[i]
    sum.trial[d[i]]<-sum.trial[d[i]]+3#������������
    int<-integrate(LF,-Inf,Inf)$value
    for (j in 1:5){
      if (sum.trial[j]!=0){EP[j]<-integrate(ep,-Inf,Inf,j)$value}}#���������������ݵļ�����Ķ��Ը���
    for (j in 1:5){
      EP[j]<0.3
         if (EP[j]<0.3){MTD<-j}#����MTD
      else break
    }
    if (d[i]<MTD){#ѡȡ��һ����������
      d[i+1]<-d[i]+1
    }else{
      if (d[i]<MTD){
        d[i+1]=max(d[i],1)
      }else {d[i+1]<-max(d[i]-1,1)}
    }
  }
  m[MTD]=1
  sim<-cbind(sum.dlt,sum.trial,m)#��¼����������
  return(sim)
}
a<-sim.crm(s,t,average,sigma,phi,n)
result.sim<-function(s,t,average,sigma,phi,n,N){
  record<-replicate(N,sim.crm(s,t,average,sigma,phi,n))#�ظ�N������
  result.sim<-matrix(rep(0,3*length(s)),nrow=5,ncol=3)
  for (i in 1:5){
    for (j in 1:3){
      result.sim[i,j]<-mean(record[i,j,])
      }
  }#����N�������DLT��������������ƽ��ֵ
  return(result.sim)
}
s<-c(0.1,0.2,0.3,0.4,0.5)#Ԥ�����ԹǼ�
t<-c(0.1,0.2,0.3,0.4,0.5)#ʵ�ʶ���
average<-0#alpha�����ֵ
sigma<-1.34#alpha�����׼��
phi<-0.3#Ŀ�����
n=10#Ԥ������������
N=1000#�������
a<-sim.crm(s,t,average,sigma,phi,n)
b<-result.sim(s,t,average,sigma,phi,n,N)