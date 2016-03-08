
library("fMultivar")
library("VGAM")
score.i<-function(y,x,z,l,theta){ # full data score (not used)
  w<-c(y[2]*x,y[3]*x,y[2],y[3],y[2]*z,y[3]*z)
	w1<-c(x,0,1,0,z,rep(0,l))
	w2<-c(0,x,0,1,rep(0,l),z)
	
	e.eta1<-exp(sum(w1*theta))
	e.eta2<-exp(sum(w2*theta))
  
  print(e.eta1)
  print(theta)
  print(sum(w2*theta))

	
	s<-w - (w1*e.eta1 + w2*e.eta2)/(1+e.eta1+e.eta2)	
	return(s)
}
score.n<-function(y,x,z,theta){ # score summed over n (not used)
	n<-length(x)
	p<-length(theta)
	s<-c(rep(0,p))
  if (is.matrix(z))
     l<-length(z[1,])
  else
     l<-1
	for(i in 1:n){
    if (l>1)
		  s<-s+score.i(y[i,],x[i],z[i,],l,theta)
	  else
	    s<-s+score.i(y[i,],x[i],z[i],l,theta)
    
	}
	return(s)
}
hessian.i<-function(x,z,theta){  # full model hessian for ith observation
    
    l<-length(z)
    w1<-c(x,0,1,0,z,rep(0,l))
  	w2<-c(0,x,0,1,rep(0,l),z)
  	
  	e.eta1<-exp(sum(w1*theta))
  	e.eta2<-exp(sum(w2*theta))
  		
  	temp<-(w1*e.eta1 + w2*e.eta2)/(1+e.eta1+e.eta2)	
  	h<- temp%*%t(temp)-((w1%*%t(w1))*e.eta1 + (w2%*%t(w2))*e.eta2)/(1+e.eta1+e.eta2)	
  
  	return(h)
}
hessian.n<-function(x,z,theta){ # score summed over n
  n<-length(x)
	p<-length(theta)
	h<-matrix(0,p,p)
	for(i in 1:n){
		h<-h+hessian.i(x[i],z[i,],theta)
	}
	return(h)
}
eff.score.i<-function(y,x,z,l,theta,I.ab,inv.I.bb){  # efficient score fixing beta_1=0
	s<-score.i(y,x,z,l,theta)
	eff.s<-t(c(s[1],s[2])-I.ab%*%inv.I.bb%*%s[3:length(theta)])
	return(eff.s)
}

# variance of efficient score
var.eff.score.n<-function(y,x,z,theta,I.ab,I.bb){
	n<-length(x)
 
	if (is.matrix(z))
	  l<-length(z[1,])
	else
	  l<-1
	v<-matrix(0,2,2)
	for(i in 1:n){
    if (l>1)
		    s<-eff.score.i(y[i,],x[i],z[i,],l,theta,I.ab,I.bb)
		else
		  s<-eff.score.i(y[i,],x[i],z[i],l,theta,I.ab,I.bb)
    #print(s)
		#var<-var+s%*%t(s)
		v[1,1]<-v[1,1]+s[1]^2
		v[1,2]<-v[1,2]+s[1]*s[2]
		v[2,1]<-v[2,1]+s[2]*s[1]
		v[2,2]<-v[2,2]+s[2]^2
	}
	return(v)
}

# efficient score
eff.score.n<-function(y,x,z,theta,I.ab,I.bb){
	n<-length(x)

	if (is.matrix(z))
	  l<-length(z[1,])
	else
	  l<-1
	s<-c(0,0)
	for(i in 1:n){
    if (l>1)
	    	s<-s+eff.score.i(y[i,],x[i],z[i,],l,theta,I.ab,I.bb)
	else
	      s<-s+eff.score.i(y[i,],x[i],z[i],l,theta,I.ab,I.bb)
	}
	return(s)
}
ComputeThreeGroup<-function(x,y,z){
    eff.s<-c(0,0)
    score.test<-c(0,0)
    var.eff.s<-array(0,dim=c(2,2))
    
    n<-length(x)

    fit.0<-vglm(y~z,multinomial(refLevel=1)) # null model
    ll.0<-fit.0@criterion$loglikelihood   
       
    theta.0<-c(0,0,fit.0@coefficients) # null theta under boundary models 

    Information<-(1/n)*(-hessian.n(x,z,theta.0)) # information matrix 
    
    
    # information submatrices needed for efficient score
    lt<-length(theta.0)
    I.ab<-Information[1:2,3:lt]
    I.bb<-solve(Information[3:lt,3:lt])
    
    eff.s<- eff.score.n(y,x,z,theta.0,I.ab,I.bb) # efficient score
    
    var.eff.s<-var.eff.score.n(y,x,z,theta.0,I.ab,I.bb) # variance estimator for efficient score
    
    var.inverse<-solve(var.eff.s)
    
    print(dim(var.inverse))
    print(dim(eff.s))
    
    if ( eff.s[1]*eff.s[2] >=0 ) {    
      l.e1<-c(1,0)%*%var.inverse%*%c(1,0)
      l.e2<-c(0,1)%*%var.inverse%*%c(0,1)
      proj.1<-eff.s%*%var.inverse%*%c(1,0)/l.e1*c(1,0)
      proj.2<-eff.s%*%var.inverse%*%c(0,1)/l.e2*c(0,1)
      
      length.1<-proj.1%*%var.inverse%*%proj.1
      length.2<-proj.2%*%var.inverse%*%proj.2
      
      score.stat<-max(length.1,length.2)      
  
}   else 
      score.stat<-eff.s%*%var.inverse%*%t(eff.s)

      
    corr12<-var.eff.s[1,2]/sqrt(var.eff.s[1,1]*var.eff.s[2,2]) # correlation between score statistics
    rho<-corr12
    p.const.score<-1/2 - asin(rho)/pi
    
    sqrt.score<-sqrt(score.stat)
    
    p.max.score<-pnorm2d(sqrt.score,sqrt.score,corr12)-pnorm2d(-sqrt.score,sqrt.score,corr12)-pnorm2d(sqrt.score,-sqrt.score,corr12)+pnorm2d(-sqrt.score,-sqrt.score,corr12)
    p.max.score<-p.max.score[[1]]    


    p.value.score<-1-(p.const.score*pchisq(score.stat,2)+(1-p.const.score)*p.max.score) # p-value of score statistic  
    return(c(p.value.score,score.stat,eff.s))
}

