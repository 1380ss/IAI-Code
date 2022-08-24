
library(matlib)
library(invgamma)
library(mvtnorm)

tscon <- function (n){
  m <- as.matrix(c(0,0,0,0,0,0,0,0))
  V <- diag(8)*100
  a1p <- 2
  a2p <- 1
  i <- 1
  
  isr <- round(runif(1))
  reward <- (round(rnorm(1,3+isr*0.5,4/3))-1)/(5-1)
  beta_coe <- c(1,isr,round(runif(6)))
  X <- beta_coe
  X <- t(as.matrix(X))

  y <- reward
  
  
  res <- y-X%*%m
  #midt <- 1/(diag(i)+X%*%V%*%t(X))
  mp <- Inverse(Inverse(V)+t(X)%*%X)%*%(Inverse(V)%*%m+t(X)%*%y)
  Vp <- Inverse(Inverse(V)+t(X)%*%X)
  a1 <- a1p + i/2
  a2 <- a2p+t(y-X%*%mp)%*%1/(diag(i)+X%*%Vp%*%t(X))%*%(y-X%*%mp)/2
  
  prec_draw <- rinvgamma(1, a1, scale = a2)
  beta_draw <- rmvnorm(1,mp,prec_draw*Vp)
  

  
  for (i in 2:n){
    
    isr <- 1*(beta_draw[2]>0)
    reward <- (round(rnorm(1,3+isr*0.5,4/3))-1)/(5-1)
    beta_coe <- c(1,isr,round(runif(6)))
    
    X <- rbind(X,beta_coe)

    
    y <- c(y,reward)
    
    
    res <- y-X%*%m
    #midt <- Inverse(diag(i)+X%*%V%*%t(X))
    mp <- Inverse(Inverse(V)+t(X)%*%X)%*%(Inverse(V)%*%m+t(X)%*%y)
    Vp <- Inverse(Inverse(V)+t(X)%*%X)
    a1 <- a1p + i/2
    a2 <- a2p+t(y-X%*%mp)%*%Inverse(diag(i)+X%*%Vp%*%t(X))%*%(y-X%*%mp)/2
    
    prec_draw <- rinvgamma(1, a1, scale = a2)
    beta_draw <- rmvnorm(1,mp,(prec_draw*Vp+t(prec_draw*Vp))/2)
    tt <- c()
    for(k in 1:1000){
      prec_draw <- rinvgamma(1, a1, scale = a2)
      beta_draw <- rmvnorm(1,mp,(prec_draw*Vp+t(prec_draw*Vp))/2)
      tt[k] <-beta_draw[2]
    }
    
  }
  return(list(X=X,y=y,m=m,V=V,tt=tt))
  
}

re <- tscon(15)  

mean(re$tt>0)



V <- diag(8)*0.01
V <- diag(8)*100
Vp <- Inverse(Inverse(V)+t(X)%*%X)
Inverse(Inverse(V)+t(X)%*%X)%*%(Inverse(V)%*%m+t(X)%*%y)
  

mean(rinvgamma(1, a1, scale = 1))


prec_draw <- rinvgamma(1, a1, scale = a2)
rmvnorm(1,mp,(prec_draw*Vp+t(prec_draw*Vp))/2)



X

df <- data.frame(cbind(y,X))
summary(lm(df,formula=y~.))


var((round(rnorm(10000,3+isr*0.5,1))-1)/(5-1))

varr <- 2/3
v0 <- 2/3
st <- sqrt(4/3)/4

mean(rnorm(100000,mean=0.5,sd=st/3)>rnorm(100000,mean=0.625,sd=st/3))

var(rnorm(9,0,1))
