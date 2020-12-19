#3+3Ä£Äâ
simulation.3<-function(d,t,n){
  a<-matrix(rep(0,10))
  x<-matrix(rep(0,30))
  r<-numeric(10)
  sim<-matrix((rep(0,2*n,row=2,col=3)))
  a[1]=1
  sum.trial<-matrix(rep(0,length(d)))
  sum.dlt<-matrix(rep(0,length(d)))
  for (i in 1:n){
    x[3*i-2]=rbinom(1,1,t[a[i]])>0
    x[3*i-1]=rbinom(1,1,t[a[i]])>0
    x[3*i]=rbinom(1,1,t[a[i]])>0
    r[i]=x[3*i-2]+x[3*i-1]+x[3*i]
    sum.dlt[a[i]]=sum.dlt[a[i]]+r[i]
    sum.trial[a[i]]=sum.trial[a[i]]+3
    if(r[i]==0 & sum.trial[min(a[i]+1,5)]<6 & sum.dlt[min(a[i]+1,5)]<3){
      a[i+1]=min(a[i]+1,length(d))
  }else{
    if(r[i]==1 & sum.dlt[a[i]]==1 & sum.trial[min(a[i]+1,5)]<6 & sum.dlt[min(a[i]+1,length(d))]<3){
      a[i+1]=min(a[i]+1,length(d))
    }else{
      if (r[i]==1 & sum.trial[a[i]]==6 &sum.dlt[a[i]]==2)
        break
      else{
        if (r[i]==1 & sum.dlt[a[i]]>2 & sum.trial[a[i]]==3){
          a[i+1]=a[i]-1
        }else{
          if(r[i]>1 &sum.trial[max(a[i-1],1)]==3){
            a[i+1]=max(a[i]-1,1)
          }else break
        }
      }
    }}
  }
  sim<-cbind(sum.dlt,sum.trial)
  return(sim)
}
simulation0<-function(d,t,n,N){
  record<-replicate(N,simulation.3(d,t,n))
  result<-matrix(rep(0,2*length(d)),nrow=5,ncol=2)
  for (i in 1:5){
    for (j in 1:2){
      result[i,j]<-mean(record[i,j,])
    }
  }
  return(result)
}
d<-c(1,2,3,4,5)
t<-c(0.10,0.20,0.30,0.40,0.50)
n=10
N=1000
result<-simulation0(d,t,n,N)