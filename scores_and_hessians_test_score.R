score.i<-function(y,x1,x2,z,theta){ # score for ith observation
	w<-y*c(x1,x2,x1*x2,1,z)
	e.eta<-exp(w%*%theta)
  score.factor<-1 - e.eta/(1+e.eta)

	
	s<-w*score.factor	
	return(s)
}



hessian.i<-function(x1,x2,z,theta){  # hessian for ith observation

	w<-(c(x1,x2,x1*x2,1,z))
 
	e.eta<-exp(w%*%theta)
	score.factor<- e.eta/(1+e.eta)^2
  temp<-w%*%t(w)
	h<- as.numeric(score.factor)*temp

	return(h)
}


score.n<-function(y,x1,x2,z,theta){ # score summed over n (not used)
	n<-length(x)
	p<-length(theta)
	s<-c(rep(0,p))
	for(i in 1:n){
		s<-s+score.i(y[i],x1[i],x2[i],z[i],theta)
	}
	return(s)
}


hessian.n<-function(x1,x2,z,n,theta){ # score summed over n

	h<-matrix(0,5,5)
  
	for(i in 1:n){
		h<-h+hessian.i(x1[i],x2[i],z[i],theta)
	}

	return(h)
}

eff.score.i<-function(y,x1,x2,z,theta,I.ab,inv.I.bb){  # efficient score 
	s<-score.i(y,x1,x2,z,theta)
	eff.s<-s[1:3]-I.ab%*%inv.I.bb%*%s[4:5]
	return(eff.s)
}



# variance of efficient score
var.eff.score.n<-function(y,x1,x2,z,n,theta,I.ab,I.bb){
	v<-matrix(0,3,3)
	for(i in 1:n){
		s<-eff.score.i(y[i],x1[i],x2[i],z[i],theta,I.ab,I.bb)
   
		v[1,1]<-v[1,1]+s[1]^2
		v[1,2]<-v[1,2]+s[1]*s[2]
		v[1,3]<-v[1,3]+s[1]*s[3]
		v[2,3]<-v[2,3]+s[2]*s[3]
    v[2,1]<-v[1,2]
		v[3,1]<-v[1,3]
		v[3,2]<-v[2,3]
		v[2,2]<-v[2,2]+s[2]^2
		v[3,3]<-v[3,3]+s[3]^2
	}
	return(v)
}

# efficient score
eff.score.n<-function(y,x1,x2,z,n,theta,I.ab,I.bb){

	s<-c(0,0,0)
	for(i in 1:n){
		s<-s+eff.score.i(y[i],x1[i],x2[i],z[i],theta,I.ab,I.bb)
	}
	return(s)
}
