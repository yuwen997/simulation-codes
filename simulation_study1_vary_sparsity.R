library(ncvreg)

n=2000
m=1000
p=50

beta<- matrix(c(1:p)/10,p,1)


si<- function(beta,delta){
  p<- length(beta)
  x=matrix(runif(p*n,1-sqrt(3),1+sqrt(3)),n,p,T)
  s=1-c(rep(1,n-m),rep(0,m))
  y=x%*%beta+diag(s)%*%x%*%delta+rnorm(n,0,1)
  true<- mean(y[1:(n-m)])
  x<- cbind(x,diag(s)%*%x)
  
  sc_value=c(1,exp(seq(log(0.2),log(2),length.out = 10)))
  erb<- rep(0,length(sc_value))
  for(sci in 1:length(sc_value)){
    sc=sc_value[sci]
    cv<- cv.ncvreg(x, y,nfolds=10,penalty="SCAD",penalty.factor=c(rep(sc,p),rep(1,p)))
    erb[sci]<- min(cv$cve)
  }
  #id1<- which(((best_model$beta !=0) ==(c(beta,delta)!=0))==F)
  sc<- sc_value[which.min(erb)]
  cv<- cv.ncvreg(x, y,nfolds=10,penalty="SCAD",penalty.factor=c(rep(sc,p),rep(1,p)))
  best_lambda <- cv$lambda.min
  best_model <- ncvfit(x, y, penalty="SCAD",lambda=best_lambda,penalty.factor=c(rep(sc,p),rep(1,p)) )
  id<- which(best_model$beta==0)
  newx<- x
  newx[,id]=0
  dat<- data.frame(y=y,newx)
  m1<- lm(y~.-1,data=dat)
  betahat<- coef(m1)
  betahat[id]=0
  
  ebeta1<- sum(( betahat[1:p]-beta)^2)
  ebeta1_penalty<- sum((best_model$beta[1:p]-beta)^2)
  cv1<- c(rep(1,p*2))
  cv1[id]=0
  
  
  
  
  
  cv<- cv.ncvreg(x, y,nfolds=10,penalty="SCAD")
  best_lambda2 <- cv$lambda.min
  best_model2 <- ncvfit(x, y, penalty="SCAD",lambda=best_lambda2)
  #id2<- which(((best_model2$beta !=0) ==(c(beta,delta)!=0))==F)
  id<- which(best_model2$beta==0)
  newx<- x
  newx[,id]=0
  dat<- data.frame(y=y,newx)
  m2<- lm(y~.-1,data=dat)
  betahat<- coef(m2)
  betahat[id]=0
  ebeta2<- sum((betahat[1:p]-beta)^2)
  ebeta2_penalty<- sum((best_model2$beta[1:p]-beta)^2)
  cv2<- c(rep(1,p*2))
  cv2[id]=0
  
  
  
  cv<- cv.ncvreg(x[(1:(n-m)),1:p], y[1:(n-m)],nfolds=10,penalty="SCAD")
  best_lambda3 <- cv$lambda.min
  best_model3 <- ncvfit(x[(1:(n-m)),1:p], y[1:(n-m)], penalty="SCAD",lambda=best_lambda3)
  id<- which(best_model3$beta==0)
  newx<- x[(1:(n-m)),1:p]
  newx[,id]=0
  dat<- data.frame(y=y[1:(n-m)],newx)
  m3<- lm(y~.-1,data=dat)
  betahat<- coef(m3)
  betahat[id]=0
  ebeta3<- sum((betahat[1:p]-beta)^2)
  ebeta3_penalty<-sum((best_model3$beta-beta)^2)
  cv3<- c(rep(1,p*2))
  cv3[id]=0
  
  
  return(list(ebeta1=ebeta1,ebeta1_penalty=ebeta1_penalty,
              ebeta2=ebeta2,ebeta2_penalty=ebeta2_penalty,
              ebeta3=ebeta3, ebeta3_penalty=ebeta3_penalty,
              cv1=cv1,cv2=cv2,cv3=cv3,sc=sc))
  
}


times=100
ans<- matrix(0,p,6)
choose2<- choose1<- matrix(0,p,2*p)
choose3<- matrix(0,p,p)
sc=rep(0,p)
var3<- var2<- var1<- list()
for(j in 1:17){
  k=seq(2,50,by=3)[j]
  if(k!=50){
    delta<- matrix(c(rep(0,k),sample(1:(p-k))),p,1)
    delta<- delta*norm(beta)/norm(delta)
  }else{
    delta<- matrix(0,p,1)
  }
  error<- matrix(0,times,6)
  cv3<- cv2<-cv1<- matrix(0,times,2*p)
  sc_values<- rep(0,times)
  
  pb <- txtProgressBar(style=3)
  for(i in 1:times){
    set.seed(i)
    tryCatch({
      result<-  si(beta,delta)
      error[i,1]<- result$ebeta1
      error[i,2]<- result$ebeta1_penalty
      error[i,3]<- result$ebeta2
      error[i,4]<- result$ebeta2_penalty
      error[i,5]<- result$ebeta3
      error[i,6]<- result$ebeta3_penalty
      cv1[i,]<- result$cv1
      cv2[i,]<- result$cv2
      cv3[i,]<- result$cv3
      sc_values[i]<-result$sc
      warnings()
      setTxtProgressBar(pb, i/times)
    },error=function(e){print(c(k,i))})
  }
  close(pb)
  var1[[j]]<- cv1
  var2[[j]]<- cv2
  var3[[j]]<- cv3
  
  ans[k,]<- colMeans(error)
  choose1[k,]<- colMeans(cv1)
  choose2[k,]<- colMeans(cv2)
  choose3[k,]<- colMeans(cv3)[1:p]
  sc[k]=mean(sc_values)
  
}

answer<- list(mse=ans,choose1=choose1,choose2=choose2,choose3=choose3,sc=sc,var1=var1,var2=var2,var3=var3)
warnings()

save(answer,file='sc_sparcity.RData')









