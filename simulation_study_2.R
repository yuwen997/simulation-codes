library(doParallel)
library(ncvreg)
library(rstan) 
times=500
n=1000
m=1000
p=2
sim<- function(){
  
  x=matrix(runif((n+m)*p,-1.5,1.5),n+m,p)
  x1=x[,1]
  x2=x[,2]
  colnames(x)<- c('x1','x2')
  a=rbinom(n,1,0.5)
  
  tau=2
  baseline = -1.5*x1^2-1.5*x2 #setting a
  #baseline = -1.5*x1^2-1.5*exp(x2) #setting b
  y1star =  baseline[1:n] + tau + rnorm(n,0,1)
  y0star =  baseline[1:n] + rnorm(n,0,1)
  y = rep(0, n+m)
  y[which(a==1)] = y1star[which(a==1)]
  
  y[which(a==0)] = y0star[which(a==0)] 
  
  s=1-c(rep(1,n),rep(0,m))
  bias=10*(x1[(n+1):(n+m)])^2+4*x2[(n+1):(n+m)]^3
  y[(n+1):(n+m)]=baseline[(n+1):(n+m)]+bias+rnorm(m,0,1)
  
  
  
  
  x_basis_function_3<-function(x_f) { 
    r<- matrix(c(x_f$x1^3,x_f$x2^3,
                 x_f$x1*(x_f$x2)^2,(x_f$x1)^2*x_f$x2),ncol=4,byrow=F)
    colnames(r)<- c('x1_3','x2_3',
                    'x1x2_2','x1_2x2'
    )
    return(r)
    
  }
  x_basis_function_2<-function(x_f) { 
    r<- matrix(c(x_f$x1^2,x_f$x2^2,
                 x_f$x1*(x_f$x2)),ncol=3,byrow=F)
    colnames(r)<- c('x1_2','x2_2',
                    'x1x2')
    return(r)
    
  }
  x_basis_function_1<-function(x_f) { 
    r<- matrix(c(x_f$x1,x_f$x2),ncol=2,byrow=F)
    colnames(r)<- c('x1','x2')
    return(r)
    
  }
  x_basis_function<- function(x_f){
    return(cbind(x_basis_function_1(x_f),
                 x_basis_function_2(x_f),
                 x_basis_function_3(x_f)))
  }
  
  x_basis<- cbind(x_basis_function(data.frame(x)),diag(s)%*%x_basis_function(data.frame(x)))
  d<- dim(x_basis)[2]/2
  data<- data.frame(y,a=c(a,rep(0,m)),x_basis,s=1-s)
  
  rct<- data[data$s==1,]
  es3=coef(lm(y~I(x1^2)+I(x2)+a,  
                  data = rct))[4]
  se3=vcov(lm(y~I(x1^2)+I(x2)+a,  
              data = rct))[4,4]
  
  es4=coef(lm(y~I(x1^2)+I(x2)+I((1-s)*(x1)^2)+I((1-s)*(x2)^3)+I(a), 
                  data = data))[6]
  se4=vcov(lm(y~I(x1^2)+I(x2)+I((1-s)*(x1)^2)+I((1-s)*(x2)^3)+I(a), 
              data = data))[6,6]
  
  
  
  
  sc_value=exp(seq(log(0.2),log(1),length.out = 50))
  erb<- rep(0,length(sc_value))
  for(sci in 1:length(sc_value)){
    sc=sc_value[sci]
    cv<- cv.ncvreg(cbind(data$a,x_basis), y,nfolds=10,penalty="SCAD",penalty.factor=c(1,rep(sc,d),rep(1,d)),max.iter = 20000)
    erb[sci]<- min(cv$cve)
  }
  
  sc<- sc_value[which.min(erb)]
  
  
  cv<- cv.ncvreg(cbind(data$a,x_basis), y,nfolds=10,penalty="SCAD",penalty.factor=c(1,rep(sc,d),rep(1,d)),max.iter = 20000)
  best_lambda <- cv$lambda.min
  best_model <- ncvfit(cbind(data$a,x_basis), y, penalty="SCAD",lambda=best_lambda,penalty.factor=c(1,rep(sc,d),rep(1,d)),max.iter = 20000 )
  
  id<- which(best_model$beta==0)
  newx<- cbind(data$a,x_basis)
  newx[,id]=0
  dat<- data.frame(y=y,data$a,newx)
  m1<- lm(y~.,data=dat)
  betahat<- coef(m1)
  es1<- coef(m1)[2]
  se1<- vcov(m1)[2,2]
  id1=id
  
  # se11<- colMeans(x_basis[1:(n-m),1:d])%*%ve[1:d,1:d]%*%(colMeans(x_basis[1:(n-m),1:d]))
  # se12<- (sum(apply(x_basis,2,var)[setdiff(1:d,id)]*(betahat[setdiff(1:d,id)]^2)))/(n-m)
  # var_beta1<- ve[1:d,1:d]
  # 
  # w1<- (es1-1.96*sqrt(se11+se12))< mean((x1)[1:(n-m)]*beta[1]+(x2^2)[1:(n-m)]*beta[2]+(x3^2)[1:(n-m)]*beta[3])
  # w2<-  (es1+1.96*sqrt(se11+se12))> mean((x1)[1:(n-m)]*beta[1]+(x2^2)[1:(n-m)]*beta[2]+(x3^2)[1:(n-m)]*beta[3])
  # wald11<- w1*w2
  
  # se13<- truemean%*%ve[1:d,1:d]%*%truemean
  # w1<- (es11-1.96*sqrt(se13))< (1*beta[1]+(5/4)*beta[2]+(5/4)*beta[3])
  # w2<-  (es11+1.96*sqrt(se13))> (1*beta[1]+(5/4)*beta[2]+(5/4)*beta[3])
  # wald12<- w1*w2
  # 
  
  cv<- cv.ncvreg(cbind(data$a[1:n],x_basis[1:n,1:d]), y[1:n],nfolds=10,penalty="SCAD",max.iter=20000)
  best_lambda3 <- cv$lambda.min
  best_model3 <- ncvfit(cbind(data$a[1:n],x_basis[1:n,1:d]), y[1:n], penalty="SCAD",lambda=best_lambda3,max.iter=20000)
  id<- which(best_model3$beta==0)
  newx<- cbind(data$a[1:n],x_basis[1:n,1:d])
  newx[,id]=0
  dat<- data.frame(y=y[1:n],newx)
  m3<- lm(y~.,data=dat)
  es2<- coef(m3)[2]
  se2<- vcov(m3)[2,2]

  
  # ve<- vcov(m3)
  # ve[,id]<- ve[id,]<- 0
  # se21<- colMeans(x_basis[1:(n-m),1:d])%*%ve[1:d,1:d]%*%(colMeans(x_basis[1:(n-m),1:d]))
  # se22<- (sum(apply(x_basis[1:(n-m),1:d],2,var)[-id]*(betahat[-id]^2)))/(n-m)
  # var_beta2<- ve
  # beta_rct<- betahat
  # 
  # w1<- (es2-1.96*sqrt(se21+se22))< mean((x1)[1:(n-m)]*beta[1]+(x2^2)[1:(n-m)]*beta[2]+(x3^2)[1:(n-m)]*beta[3])
  # w2<-  (es2+1.96*sqrt(se21+se22))> mean((x1)[1:(n-m)]*beta[1]+(x2^2)[1:(n-m)]*beta[2]+(x3^2)[1:(n-m)]*beta[3])
  # wald21<- w1*w2
  
  
  # se23<- truemean%*%ve[1:d,1:d]%*%truemean
  # w1<- (es22-1.96*sqrt(se23))< (1*beta[1]+(5/4)*beta[2]+(5/4)*beta[3])
  # w2<-  (es22+1.96*sqrt(se23))> (1*beta[1]+(5/4)*beta[2]+(5/4)*beta[3])
  # wald22<- w1*w2
  return(list(es1=es1,es2=es2,es3=es3,es4=es4,se1=se1,se2=se2,se3=se3,se4=se4))
  
}

registerDoParallel(32)
# error<- foreach(i=1:times,.combine=rbind) %dopar% {
#   set.seed(i)
#   tryCatch({
#     result<-  sim()
#     
#   },error=function(e){print(c(i))})
# }
error<- foreach(i=1:times,.combine=rbind) %dopar% {
  set.seed(i)
  
    result<-  sim()
    
  
}
stopImplicitCluster()
# rr<- matrix(unlist(error),times,8,F)
# apply(rr,2,mean)
# apply(rr,2,var)
# save(error,file='rebuttal.RData')



