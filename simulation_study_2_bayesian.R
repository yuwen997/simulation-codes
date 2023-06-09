library(rstan)
library(R2jags)
library(MASS)
library(doParallel)
times=500
n=1000
m=1000
p=2
bayes<- function(){
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
  
  power_prior_weighted <- "data{
        for (i in 1:m){
          zeros[i] <- 0
        }
        k <- 10000
      }
      model{
        # Likelihood for RCT data
        for(i in 1:n_d){
          Y_D[i] ~ dnorm(mu_D[i],inv.var)
          mu_D[i] <- beta0 + inprod(beta[], X_D[i,]) + beta_T*T_D[i]
        }
        # Log-Likelihood for RW data using Zeros trick
        for (i in 1:m){
          mu_K[i] <- beta0 + inprod(beta[], X_K[i,]) + beta_T*T_K[i]
          loglik_K[i] <- alpha[i]*log(dnorm(Y_K[i], mu_K[i], inv.var))
          phi[i] <- -loglik_K[i] + k
          zeros[i] ~ dpois(phi[i])
        }

        # Prior for Betas
        for(p in 1:P){
          beta[p] ~ dnorm(0,0.0001)
        }
        
        beta0 ~ dnorm(0, 0.0001)
        beta_T ~ dnorm(0, 0.0001)

        # Prior for variance
        sigma ~ dunif(0.01, 500)
        inv.var <- 1/(sigma*sigma)
      }"
  alpha<- predict(glm(s~.,data[,c(12:13,ncol(data))],family='binomial'),data[(n+1):(n+m),12:13])
  parameters <- c("beta0", "beta", "beta_T", "sigma")
  inits=function(){list(beta0=rnorm(1,0,1), beta=mvrnorm(1 ,rep(0,9),diag(9)), beta_T=rnorm(1,0,1), sigma=runif(1,1,500))}
  sim_data=list(Y_D =y[1:n], Y_K=y[(n+1):(n+m)],
                X_D=as.matrix(data[1:n,3:11]), X_K=as.matrix(data[(n+1):(n+m),3:11]),
                T_D=data$a[1:n], T_K=data$a[(n+1):(n+m)],
                n_d=1000,alpha=alpha,m=m,P=9)  #alpha=subject specific PS
  
  pow.prior_weighted=jags(sim_data, parameters, textConnection(power_prior_weighted), inits=inits, n.chains = 2, n.thin = 1, n.iter = 2000, n.burnin = 1000)
  return(pow.prior_weighted$BUGSoutput$summary[11,])
}



registerDoParallel(32)
# error<- foreach(i=1:times,.combine=rbind) %dopar% {
#   set.seed(i)
#   tryCatch({
#     result<-  sim()
#     
#   },error=function(e){print(c(i))})
# }

registerDoParallel(32)
error<- foreach(i=1:times,.combine=rbind) %dopar% {
  set.seed(i)
  
  result<-  bayes()
  
  
}
stopImplicitCluster()
rr<- matrix(unlist(error),times,9,F)
apply(rr,2,mean)
apply(rr,2,var)
save(error,file='bayeian1.RData')



#alpha.d=as.data.frame(matchdata_rwd$distance)

