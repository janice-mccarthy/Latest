samplingScore.function<-function(nsim,n0,n1,beta.1,beta.2,beta.3,alpha.1,gamma.1,freq,alfa){
  
  in.region<-0
  # initialize matrices and vectors
  p.const<-rep(0,nsim)
  const<-rep(0,nsim)
  
  score.stat<-rep(0,nsim)
  score.vector<-matrix(0,3,nsim)
#  score.test<-rep(0,nsim)
  #test.vector<-matrix(0,nsim,2)
  #test.product<-rep(0,nsim)
  
  
  corr12<-rep(0,nsim)
  corr13<-rep(0,nsim)
  corr23<-rep(0,nsim)
#  var.eff.s<-array(0,dim=c(3,3,nsim))
  p.value.score<-rep(0,nsim)
  n<-n0+n1
  
for(j in 1:nsim){   # loop over simulated datasets
    
    # initialize data matrices
    y<-matrix(0,n)
    z<-rep(0,n)
    x1<-rep(0,n)
    x2<-rep(0,n)
    i<-1
    i0<-0
    i1<-0
    
    while( (i0<n0) | (i1<n1) ){ # retrospective sampling of groups
      
      z[i]<-rnorm(1)  # covariate
      x1[i]<-rbinom(1,2,freq) # genotype
      x2[i]<-rbinom(1,2,freq) # genotype
      
      # generate Case or Control group conditional on x and z
      pi.1<-exp(alpha.1+beta.1*x1[i]+beta.2*x2[i]+beta.3*x1[i]*x2[i]+gamma.1*z[i])
      
      denom<-1+pi.1
      pi.1<-pi.1/denom # prob of case
      
      y[i]<-rbinom(1,1,pi.1) # Case/Control Status
      
      if((y[i]==0) & (i0<n0)){ # keep control if still needed
        i<-i+1
        i0<-i0+1
      }
      else{
        if((y[i]==1) & (i1<n1)){ # keep case if still needed
          i<-i+1
          i1<-i1+1
        }
      }
      
    } # end sampling loop
    #
    
    # scale genotypes
    x1<-(x1-mean(x1))/sd(x1)
    x2<-(x2-mean(x2))/sd(x2)
    
    # fit null model
  
    fit.0<-glm(y~z,family="binomial")
    nuisance<-fit.0$coefficients
    
    theta.0<-c(as.numeric(nuisance[1]),as.numeric(nuisance[2])) # null theta under full model
    e.eta<-rep(0,n)
  
    # The following only makes sense for genome-wide calculations
    for (k in 1:n){
      w<-c(1,z[k])    
      e.eta[k]<-exp(w%*%theta.0) 
    }
    
    Information<-(1/n)*(hessian.n(x1,x2,z,n,e.eta)) # information matrix (sign is accounted for in hessian.i function)
  
    var<-solve(Information) # variance under null

    rho.12<-var[1,2]/sqrt(var[1,1]*var[2,2])
    rho.13<-var[1,3]/sqrt(var[1,1]*var[3,3])
    rho.23<-var[2,3]/sqrt(var[2,2]*var[3,3])
    
    weight <-1/(2) - (1/(4*pi))*(acos(rho.12)+acos(rho.13)+acos(rho.23)) # weight for chi-bar-squared statistic
   # print(weight)
    
    # information submatrices needed for efficient score
    I.ab<-Information[1:3,4:5]
    I.bb<-solve(Information[4:5,4:5])
  
    eff.s<-eff.score.n(y,x1,x2,z,n,e.eta,I.ab,I.bb) # efficient score
    
    var.eff.s<-var.eff.score.n(y,x1,x2,z,n,e.eta,I.ab,I.bb) # variance estimator for efficient score
   
    inv.var.eff.s<-solve(var.eff.s)
  
    # If the efficient score vector is in the positive octant, the statisic is its (normalized) length. Otherwise, we project onto
    # the positive octant (using the variance form) and then take the length
    
    if (eff.s[1]>0&eff.s[2]>0&eff.s[3]<0){  # project onto 1-2 plane by subtracting component in direction 3         
      l.proj.direction<-c(0,0,1)%*%inv.var.eff.s%*%c(0,0,1)
      proj.score<- eff.s - (t(eff.s)%*%inv.var.eff.s%*%c(0,0,1)/l.proj.direction)*c(0,0,1)
      score.stat[j]<-t(proj.score)%*%inv.var.eff.s%*%proj.score  
    }
    else
    {
     if(eff.s[1]>0&eff.s[2]<0&eff.s[3]>0){ # project onto 1-3 plane 
       l.proj.direction<-c(0,1,0)%*%inv.var.eff.s%*%c(0,1,0)
       proj.score<- eff.s - t(eff.s)%*%inv.var.eff.s%*%c(0,1,0)/l.proj.direction*c(0,1,0)
       score.stat[j]<-t(proj.score)%*%inv.var.eff.s%*%proj.score    
      }
    else 
    { 
     if(eff.s[1]<0&eff.s[2]>0&eff.s[3]>0){ # project onto 2-3 plane {
      l.proj.direction<-c(1,0,0)%*%inv.var.eff.s%*%c(1,0,0)
      proj.score<- eff.s - t(eff.s)%*%inv.var.eff.s%*%c(1,0,0)/l.proj.direction*c(1,0,0)
      score.stat[j]<-t(proj.score)%*%inv.var.eff.s%*%proj.score    
      
    }
    else
    { 
      if (eff.s[1]>0&eff.s[2]>0 &eff.s[3]>0) { # in region 
      score.stat[j]<-t(eff.s)%*%inv.var.eff.s%*%eff.s
      in.region<-in.region +1
    }
    else
    { 
      if (eff.s[1]>0&eff.s[2]<0&eff.s[3]<0){  # project onto 1-d subspace         
        l.proj.direction<-c(1,0,0)%*%inv.var.eff.s%*%c(1,0,0)
        proj.score<- (t(eff.s)%*%inv.var.eff.s%*%c(1,0,0)/l.proj.direction)*c(1,0,0)
        score.stat[j]<-t(proj.score)%*%inv.var.eff.s%*%proj.score  
      }
  
    else
    { 
      if (eff.s[1]<0&eff.s[2]>0&eff.s[3]<0){  # project onto 1-d subspace         
        l.proj.direction<-c(0,1,0)%*%inv.var.eff.s%*%c(0,1,0)
        proj.score<- (t(eff.s)%*%inv.var.eff.s%*%c(0,1,0)/l.proj.direction)*c(0,1,0)
        score.stat[j]<-t(proj.score)%*%inv.var.eff.s%*%proj.score  
      }
    
    else
    { 
      if (eff.s[1]<0&eff.s[2]<0&eff.s[3]>0){  # project onto 1-d subspace         
        l.proj.direction<-c(0,0,1)%*%inv.var.eff.s%*%c(0,0,1)
        proj.score<- (t(eff.s)%*%inv.var.eff.s%*%c(0,0,1)/l.proj.direction)*c(0,0,1)
        score.stat[j]<-t(proj.score)%*%inv.var.eff.s%*%proj.score  
      }
    
    else
    {
      score.stat[j]<-0
      
    }
    }}}}}}
    #  p-value for score statistic
    score.vector[1,j]<-eff.s[1]
    score.vector[2,j]<-eff.s[2]
    score.vector[3,j]<-eff.s[3]
    
    p.value.score[j]<-weight*(1-pchisq(score.stat[j],3)) + 
     0.5*(1-weight)*(1-pchisq(score.stat[j],2)) + 0.5*(1-weight)*(1-pchisq(score.stat[j],1)) + weight*(1-pchisq(score.stat[j],0))
     
    # corr12[j]<-var.eff.s[1,2,j]/sqrt(var.eff.s[1,1,j]*var.eff.s[2,2,j]) # correlation between score statistics
    if (j%%100 == 0)
        print(j)
    
  } # end simulation loop
 
  score.test<-sum(p.value.score<alfa)/nsim
  
#  return(c(in.region,score.test,score.stat,score.vector))
   return(score.test)
} #end function definition
