score.i<-function(y,x1,x2,z,e.eta){ # score for ith observation
	w<-c(x1,x2,x1*x2,1,z)

  score.factor<-(y - e.eta/(1+e.eta))

  s<-w*score.factor
 
	return(s)
}

hessian.i<-function(x1,x2,z,e.eta){  # hessian for ith observation

	w<-(c(x1,x2,x1*x2,1,z))
 
	score.factor<- e.eta/(1+e.eta)^2 # There is a negative sign in the hessian that is already accounted for
  temp<-w%*%t(w)
	h<- as.numeric(score.factor)*temp

	return(h)
}

score.n<-function(y,x1,x2,z,n,e.eta){ # score summed over n 
	s<-c(rep(0,5))
	for(i in 1:n){
		s<-s+score.i(y[i],x1[i],x2[i],z[i],e.eta[i])
	}
	return(s)
}


hessian.n<-function(x1,x2,z,n,e.eta){ # score summed over n

	h<-matrix(0,5,5)
  
	for(i in 1:n){
		h<-h+hessian.i(x1[i],x2[i],z[i],e.eta[i])
	}

	return(h)
}

eff.score.i<-function(y,x1,x2,z,e.eta,I.ab,inv.I.bb){  # efficient score 
	s<-score.i(y,x1,x2,z,e.eta)
	eff.s<-s[1:3]-I.ab%*%inv.I.bb%*%s[4:5]
	return(eff.s)
}



# variance of efficient score
var.eff.score.n<-function(y,x1,x2,z,n,e.eta,I.ab,I.bb){
	v<-matrix(0,3,3)
	
  for(i in 1:n){
		s<-eff.score.i(y[i],x1[i],x2[i],z[i],e.eta[i],I.ab,I.bb)
   
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
eff.score.n<-function(y,x1,x2,z,n,e.eta,I.ab,I.bb){

	s<-c(0,0,0)
	for(i in 1:n){
		s<-s+eff.score.i(y[i],x1[i],x2[i],z[i],e.eta[i],I.ab,I.bb)
	}
	return(s)
}
