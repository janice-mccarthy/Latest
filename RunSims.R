source("scores_and_hessians_Epistasis.R")
source("EpistasisSampling.R")

nsim<-10000
n0<-1000
n1<-1000
beta.1<-0.0
beta.2<-0.0
beta.3<-0.0
alpha.1<-0.0
gamma.1<-.1
freq<-0.30
alfa<-0.05
null.val<-samplingScore.function(nsim,n0,n1,beta.1,beta.2,beta.3,alpha.1,gamma.1,freq,alfa)

nsim<-1000
n0<-1000
n1<-1000
alternates<-rep(0,5)
for(k in 1:5){
  beta.1<-log(1+k/20)
  beta.2<-log(1+k/30)
  beta.3<-log(1+k/30)
  alternates[k]<-samplingScore.function(nsim,n0,n1,beta.1,beta.2,beta.3,alpha.1,gamma.1,freq,alfa)
}