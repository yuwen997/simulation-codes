
library(doParallel)
times=500
n=1000
m=1000
p=2
M=1
match<- function(){
  x=matrix(runif((n+m)*p,-1.5,1.5),n+m,p)
  x1=x[,1]
  x2=x[,2]
  colnames(x)<- c('x1','x2')
  a=rbinom(n,1,0.5)
  
  tau=2
  #baseline = -1.5*x1^2-1.5*x2 #setting a
  baseline = -1.5*x1^2-1.5*exp(x2) #setting b
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
  data<- data.frame(y,a=c(a,rep(0,m)),x_basis[,1:d],s=s)
  ps<- predict(glm(a~.-s-y,data=data))
  data$ps<- ps
  
  cc<- data[(data$a==0) & (data$s==0),]
  trt<- data[(data$a==1) & (data$s==0),]
  hc<- data[(data$a==0) & (data$s==1),]
  caliper=0.2*sd(ps[data$a==1])
  
  match1<- rbind(trt,cc)
  m1<- Matching::Match(Y = match1$y, Tr = match1$a, X =match1$ps ,caliper=caliper,estimand = "ATT",
                  distance.tolerance = 0, Weight = 2,replace=F)
  unmatch<- 1:dim(trt)[1]
  unmatch<- unmatch[-m1$index.treated]
  

  match2<- rbind(cbind(match1[m1$index.control,],h=1),cbind(hc,h=0))
  m2<- Matching::Match(Y = match2$y, Tr = match2$h, X =match2$ps ,M=M,estimand = "ATT",
                       distance.tolerance = 0, Weight = 2,replace=F)
  
  match3<- rbind(cbind(trt[unmatch,],h=1),cbind(hc,h=0))
  m3<- Matching::Match(Y = match3$y, Tr = match3$h, X =match3$ps ,M=M,estimand = "ATT",
                       distance.tolerance = 0, Weight = 2,replace=F)
  
  biase<- match2[c(m2$index.treated,m2$index.control),]
  modelb<- summary(lm(y~.-a-h-ps,data=biase))
  
  
  

  
  repeat_times=10
  uncentainty<- rep(  repeat_times,0)
  for(l in 1:repeat_times){
  s2<- 1/rchisq(1, df=length(m2$index.treated)-(d+2),ncp=(modelb$sigma)^2)
  xfull<-biase[,dim(biase)[2]-1]
  dd<-rnorm(1,mean=(modelb$coefficients[d+2,1]),sd=sqrt(solve(t(xfull)%*%(xfull))*s2))
  
  newc<- rbind(match1[m1$index.control,],match3[m3$index.control,-dim(match3)[2]])
  ind<- c(rep(1,length(m1$index.control)),rep(2,length(m3$index.control))) 
  y0<- newc[,1]
  y0[ind==2]<- y0[ind==2]-dd
  ate<- mean(trt$y)-mean(y0)
  uncentainty[l]<- ate
  }
ate1<- mean(uncentainty)

modelb<- summary(lm(y~.-a-ps,data=rbind(cc,hc)))
repeat_times=10
uncentainty<- rep(  repeat_times,0)
for(l in 1:repeat_times){
  s2<- 1/rchisq(1, df=length(m2$index.treated)-(d+2),ncp=(modelb$sigma)^2)
  xfull<-biase[,dim(biase)[2]-1]
  dd<-rnorm(1,mean=(modelb$coefficients[d+2,1]),sd=sqrt(solve(t(xfull)%*%(xfull))*s2))
  
  newc<- rbind(match1[m1$index.control,],match3[m3$index.control,-dim(match3)[2]])
  ind<- c(rep(1,length(m1$index.control)),rep(2,length(m3$index.control))) 
  y0<- newc[,1]
  y0[ind==2]<- y0[ind==2]-dd
  ate<- mean(trt$y)-mean(y0)
  uncentainty[l]<- ate
}
ate2<- mean(uncentainty)
 
return(c(ate1,ate2))
}



registerDoParallel(16)

error<- foreach(i=1:times,.combine=rbind) %dopar% {
  set.seed(i)
  
  result<-  match()
  
  
}
stopImplicitCluster()
rr<- matrix(unlist(error),times,2,F)
apply(rr,2,mean)-2
apply(rr,2,var)
save(error,file='/Users/yuwencheng/Desktop/fourth/written/setting_b_matching.RData')



