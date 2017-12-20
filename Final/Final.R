################### Problem A: Asian Option Pricing using Monte Carlo Control Variate ###################
############(a)#############
# Function to calculate opntion price by BS formular
BS.Asian.geo <- function(S0, K, T, r, sigma,N){
  
  dt = T/N
  adjust.sigma <- sigma * sqrt((2*N+1)/(6*(N+1)))
  rho = 0.5 * (r-1/2*sigma^2+adjust.sigma^2)
  d1 = (log(S0/K)+(rho+adjust.sigma^2/2)*T) / (adjust.sigma*sqrt(T)) 
  d2 = (log(S0/K)+(rho-adjust.sigma^2/2)*T) / (adjust.sigma*sqrt(T))

  Price = exp(-r*T)*(S0*exp(rho*T)*pnorm(d1) - K*pnorm(d2))
  return(Price)
}
BS.Asian.geo(S0=100, K=100, T=5, r=0.03, sigma=0.3,N=5*252)

############(b)#############
MC.Asian.arith <- function(S0, K, T, r, sigma,N,M){
  
  timestart<-Sys.time()
  #precompute constants
  dt = T/N
  
  sum_CT = 0
  sum_CT2 = 0
  
  St <- matrix(0,nrow = N+1,ncol = 1)
  for(j in 1:M){ #for each simulation
    St[1] = S0
    for(i in 1:N){ #for each time step
      ybi = rnorm(1)
      St[i+1] = St[i]*exp((r-0.5*sigma^2)*dt+sigma*sqrt(dt)*ybi)
    }
    ST = sum(St)/(N+1)
    CT = max(0, ST - K)

    sum_CT = sum_CT + CT
    sum_CT2 = sum_CT2 + CT^2
  }
  value = sum_CT/M*exp(-r*T) #discount option value to time 0
  SD = sqrt((sum_CT2 - sum_CT^2/M)*exp(-2*r*T)/(M-1))
  SE = SD/sqrt(M) #standard error
  CI = c(value - 1.96*SE,value + 1.96*SE)
  return(list(price = value, CI  = CI))
  timeend<-Sys.time()
  runningtime<-timeend-timestart
  print(runningtime)
}
MC.Asian.arith(S0=100, K=100, T=5, r=0.03, sig=0.3,N=5*252,M=1e+06)

#############(c)#############
MC.Asian.geo <- function(S0, K, T, r, sigma,N,M){
  
  #precompute constants
  dt = T/N
  
  sum_CT = 0
  
  St <- matrix(0,nrow = N+1,ncol = 1)
  for(j in 1:M){ #for each simulation
    St[1] = S0
    for(i in 1:N){ #for each time step
      ybi = rnorm(1)
      St[i+1] = St[i]*exp((r-0.5*sigma^2)*dt+sigma*sqrt(dt)*ybi)
    }
    ST = prod(St)^(1/(N+1))
    CT = max(0, ST - K)
    
    sum_CT = sum_CT + CT
  }
  value = sum_CT/M*exp(-r*T) #discount option value to time 0
  return(value)
}
MC.Asian.geo(S0=100, K=100, T=5, r=0.03, sig=0.3,N=5*252,M=1e+06)

#############(d)#############
b.star <- function(S0, K, T, r, sigma,M){
  
  #precompute constants
  N = 100
  dt = T/N
  
  X=Y=c();geo.ST=arith.ST=c()
  St <- matrix(0,nrow = N+1,ncol = 1)
  for(j in 1:M){ #for each simulation
    St[1] = S0
    for(i in 1:N){ #for each time step
      ybi = rnorm(1)
      St[i+1] = St[i]*exp((r-0.5*sigma^2)*dt+sigma*sqrt(dt)*ybi)
    }
    ST.arith = sum(St)/(N+1)
    ST.geo = prod(St)^(1/(N+1))
    X = c(X,max(0,ST.geo - K)) #geometric option price
    Y = c(Y,max(0, ST.arith - K)) #arithmetic option price
    geo.ST = c(geo.ST,ST.geo)
    arith.ST = c(arith.ST,ST.arith)
  }
  b = sum((X-mean(X))*(Y-mean(Y)))/sum((X-mean(X))^2)
  return(list(b=b,geo.price=geo.ST,arith.price=arith.ST))
}
b <- b.star(S0=100, K=100, T=5, r=0.03, sigma=0.3,M=1e+06)
b$geo.price;b$arith.price;b$b

#############(e)##############
error <- function(S0,K,T,r,sig,N,M){
  simulate.price <- MC.Asian.geo(S0, K, T, r, sig,N,M)
  theo.price <- BS.Asian.geo(S0, K, T, r, sig,N)
  error.geo <- theo.price - simulate.price
  return(error.geo)
}
error(S0=100,K=100,T=5,r=0.03,sig=0.3,N=100,M=1e+03)

#############(f)##############
find.M <- function(S0,K,T,r,sig,N,M){
  price.arith.star = c()
  for(i in 1:length(M)){
    price.arith.simu <- MC.Asian.arith(S0, K, T, r, sig,N,M[i])$price
    star.b <- b.star(S0, K, T, r, sig, M[i])$b
    err <- price.arith.simu - star.b * error(S0,K,T,r,sig,N,M[i])
    price.arith.star <- c(price.arith.star, err)
  }
  return(price.arith.star)
}
price.arith.star <- find.M(S0=100,K=100,T=5,r=0.03,sig=0.3,N=100,M=c(1e+03,1e+04,1e+05))
price.arith.star

###########bonus##############
AMZN <- read.csv("amzn.csv")
S0 <- 952.82
K <- AMZN$Strike
T <- AMZN$T
sigma <- AMZN$Volatility

amzn.bs.geo=amzn.arith=amzn.geo=amzn.err=amzn.b=amzn.m=c()
for(i in 1:25){
    amzn.bs.geo <- c(amzn.bs.geo,BS.Asian.geo(S0=S0,K=K[i],T=T[i],r=0.03,sigma=sigma[i],N=10))
    amzn.arith <- c(amzn.arith,MC.Asian.arith(S0=S0,K=K[i],T=T[i],r=0.03,sig=sigma[i],N=10,M=1e+03)$price)
    amzn.geo <- c(amzn.geo,MC.Asian.geo(S0=S0,K=K[i],T=T[i],r=0.03,sig=sigma[i],N=10,M=1e+03))
    amzn.b <- c(amzn.b,b.star(S0=S0,K=K[i],T=T[i],r=0.03,sigma=sigma[i],M=1e+03)$b)
    amzn.err <- c(amzn.err,error(S0=S0,K=K[i],T=T[i],r=0.03,sig=sigma[i],N=10,M=1e+03))
    amzn.m <- c(amzn.m,find.M(S0=S0,K=K[i],T=T[i],r=0.03,sig=sigma[i],N=100,M=1e+03))
}
data.frame(amzn.bs.geo,amzn.arith,amzn.geo,amzn.b,amzn.err,amzn.m)
  
################### Problem B: A portfolio construction problem ###################
############ 1 #############
library(quantmod)
stockData <- new.env()
XLF.symbols <- c('ACC','JPM','WFC','BAC','C','GS','USB','CB','MS','AXP',
                 'AADR','ACFC','AAME','AAT','AAXJ','ACNB','ABCB','ABE','ABR','BRK-B')
getSymbols(XLF.symbols, from="2012-01-01", env=stockData, src="yahoo")

ReturnMatrix=XLF.data=NULL
for(i in 1:length(XLF.symbols)){
  temp = get(XLF.symbols[i], pos=stockData)
  XLF.data = cbind(XLF.data,Cl(temp))
  ReturnMatrix = cbind(ReturnMatrix, (Cl(temp)-Op(temp)) / Op(temp)   )
  colnames(ReturnMatrix)[i] = XLF.symbols[i]
}
PCA <- princomp(ReturnMatrix,cor = TRUE) #by default R centers the variables. Scale also makes then sd=1
summary(PCA)
PCA$loadings

############ 2 #############
library(psych)
fit <- principal(cor(ReturnMatrix), nfactors=10, rotate="varimax", n.obs=dim(ReturnMatrix)[1])
fit$loadings

library(Sim.DiffProc)
pick <- c("AXP.Close","AAXJ.Close","ABCB.Close","AAT.Close")
equity.pick <- XLF.data[,pick]

model.select <- function(data){
  #model1 Black-Scholes
  fx1 <- expression(theta[1]*x)
  gx1 <- expression(theta[2]*x)
  mod1 <- fitsde(data=as.ts(data),drift=fx1,diffusion=gx1,
                 start = list(theta1=1,theta2=1),pmle="shoji")

  #model2 mean reverting CEV
  fx2 <- expression(theta[1]+theta[2]*x)
  gx2 <- expression(theta[3]*x^theta[4])
  mod2 <- fitsde(data=as.ts(data),drift=fx2,diffusion=gx2,
                 start = list(theta1=1,theta2=1,theta3=1,theta4=1),pmle="shoji")

  #model3 Strange 1
  fx3 <- expression(theta[1]*x)
  gx3 <- expression(theta[2]+theta[3]*x^theta[4])
  mod3 <- fitsde(data=as.ts(data),drift=fx3,diffusion=gx3,
                 start = list(theta1=1,theta2=1,theta3=1,theta4=1),pmle="shoji")
  
  #model4 particular CEV
  fx4 <- expression(theta[1]*x)
  gx4 <- expression(theta[2]*x^(2/3))
  mod4 <- fitsde(data=as.ts(data),drift=fx4,diffusion=gx4,
                 start = list(theta1=1, theta2=1),pmle="shoji")
  
  #model5 Strange 2
  fx5 <- expression(theta[1]+theta[2]*x)
  gx5 <- expression((theta[3]+theta[4]*log(x))*x)
  mod5 <- fitsde(data=as.ts(data),drift=fx5,diffusion=gx5,
                 start = list(theta1=1, theta2=1,theta3=1,theta4=1),pmle="shoji")
  
  ## Computes AIC
  AIC <- c(AIC(mod1),AIC(mod2),AIC(mod3),AIC(mod4),AIC(mod5))
  Test <- data.frame(AIC,row.names = c("mod1","mod2","mod3","mod4","mod5"))
  Bestmod <- rownames(Test)[which.min(Test[,1])]
  cat("best model: ",Bestmod,"\n")
  print(Test)
}
col1.AIC <- model.select(equity.pick[,1])  #mod2 AXP
col2.AIC <- model.select(equity.pick[,2])  #mod2 AAXJ
col3.AIC <- model.select(equity.pick[,3])  #mod1 ABCB
col4.AIC <- model.select(equity.pick[,4])  #mod2 AAT

############ 3 #############
mat.cor <- cor(equity.pick)
mat.cor

############ 4 #############
#model1 Black-Scholes
fx1 <- expression(theta[1]*x)
gx1 <- expression(theta[2]*x)
mod1 <- fitsde(data=as.ts(equity.pick[,3]),drift=fx1,diffusion=gx1,
               start = list(theta1=1,theta2=1),pmle="shoji")
ABCB.coef <- coef(mod1)

#model2 mean reverting CEV
fx2 <- expression(theta[1]+theta[2]*x)
gx2 <- expression(theta[3]*x^theta[4])
mod21 <- fitsde(data=as.ts(equity.pick[,1]),drift=fx2,diffusion=gx2,
               start = list(theta1=1,theta2=1,theta3=1,theta4=1),pmle="shoji") #AXP
mod22 <- fitsde(data=as.ts(equity.pick[,2]),drift=fx2,diffusion=gx2,
               start = list(theta1=1,theta2=1,theta3=1,theta4=1),pmle="shoji") #AAXJ
mod23 <- fitsde(data=as.ts(equity.pick[,1]),drift=fx2,diffusion=gx2,
               start = list(theta1=1,theta2=1,theta3=1,theta4=1),pmle="shoji") #AAT
AXP.coef <- coef(mod21)
AAXJ.coef <- coef(mod22)
AAT.coef <- coef(mod23)

MC.mod <- function(S0, T){
  #precompute constants
  N = 255
  M = 1e+03 #1000 paths
  dt = T/N
  
  set.seed(1)
  ST1=ST2=ST3=ST4=c()
  for(j in 1:M){ #for each simulation
    x <- matrix(rnorm(4*N),nrow=N,ncol=4)
    ep <- x%*%chol(mat.cor)
    e1 <- append(0,ep[,1])
    e2 <- append(0,ep[,2])
    e3 <- append(0,ep[,3])
    e4 <- append(0,ep[,4])
    S1=S0[,1] #initial value for asset 1
    S2=S0[,2] #initial value for asset 2
    S3=S0[,3] #initial value for asset 3
    S4=S0[,4] #initial value for asset 3
    for(i in 1:(N+1)){ #for each time step
      S1 <- S1 + (AXP.coef[1]-AXP.coef[2]*S2)*dt + e1[i]*(AXP.coef[3]*S2^AXP.coef[4])*sqrt(dt) #AXP
      S2 <- S2 + (AAXJ.coef[1]-AAXJ.coef[2]*S2)*dt + e2[i]*(AAXJ.coef[3]*S2^AAXJ.coef[4])*sqrt(dt) #AAXJ
      S3 <- S3*exp((ABCB.coef[1]-0.5*ABCB.coef[2]^2)*dt + e3[i]*ABCB.coef[2]*sqrt(dt)) #AAT
      S4 <- S4 + (AAT.coef[1]-AAT.coef[2]*S2)*dt + e4[i]*(AAT.coef[3]*S2^AAT.coef[4])*sqrt(dt) #AXP #ABCB
    }
    ST1 <- c(ST1,S1)
    ST2 <- c(ST2,S2)
    ST3 <- c(ST3,S3)
    ST4 <- c(ST4,S4)
  }
  ST.mean <- c(mean(ST1),mean(ST2),mean(ST3),mean(ST4))
  ST.sd <- c(sd(ST1),sd(ST2),sd(ST3),sd(ST4))
  ST.skewness <- c(skewness(ST1),skewness(ST2),skewness(ST3),skewness(ST4))
  ST.kurtosis <- c(kurtosis(ST1),kurtosis(ST2),kurtosis(ST3),kurtosis(ST4))
  show <- data.frame(ST.mean,ST.sd,ST.skewness,ST.kurtosis)
  print(show)
}
MC.mod(S0=equity.pick[1,],T=1)

############ 5 #############
XLF <- read.csv("XLF.csv")
fx <- expression(theta[1]*x)
gx <- expression(theta[2]*x)
mod <- fitsde(data=as.ts(XLF$Close),drift=fx,diffusion=gx,
               start = list(theta1=1, theta2=1),pmle="shoji")
mu <- coef(mod)[1]
sigma <- coef(mod)[2]
mu;sigma;

############ 6 #############
XLF.return <- (XLF$Close-XLF$Open) / XLF$Open
newdata <- data.frame(XLF.return,equity.pick[,1],equity.pick[,2],equity.pick[,3],equity.pick[,4])
fit <- lm(XLF.return~., data = newdata)
weight <- coef(fit)
weight

############ 7 #############
basket <- function(type,s0,ST,T,m,n){
  
  nu = mu-sigma^2/2
  nu = as.numeric(nu)
  sum_CT = 0
  dt = T/n
  for(i in 1:m){   
    x <- matrix(rnorm(400),nrow=100,ncol=4)
    ep <- x%*%chol(mat.cor)
    e1 <- append(0,ep[,1])
    e2 <- append(0,ep[,2])
    e3 <- append(0,ep[,3])
    e4 <- append(0,ep[,4])
    s1=s2=s3=s4=c()
    S1=s0[1] #initial value for asset 1
    S2=s0[2] #initial value for asset 2
    S3=s0[3] #initial value for asset 3
    S4=s0[4] #initial value for asset 3
    for(j in 1:(n+1)){
      S1 <- S1*exp(nu*dt+e1[j]*as.numeric(sigma)*sqrt(dt))
      S2 <- S2*exp(nu*dt+e2[j]*as.numeric(sigma)*sqrt(dt))
      S3 <- S3*exp(nu*dt+e3[j]*as.numeric(sigma)*sqrt(dt))
      S4 <- S4*exp(nu*dt+e4[j]*as.numeric(sigma)*sqrt(dt))
      s1 <- c(s1,S1)
      s2 <- c(s2,S2)
      s3 <- c(s3,S3)
      s4 <- c(s4,S4)
    }
    
    U = weight[2]*s1+weight[3]*s2+weight[4]*s3+weight[5]*s4 #price of basket stocks
    if(type=="ETF"){
      CT = max(0,U[101]-ST)
    }
    else if(type=="equity"){
      CT = max(0,ST - U[101])
    }
    sum_CT = sum_CT + CT
  }
  value = sum_CT/m
  return(value)
}
n=length(XLF$Close)
basket(type = "ETF",s0 = as.numeric(equity.pick[1,]),ST = XLF$Close[n],T=1,m=1000,n=100)
basket(type = "equity",s0 = as.numeric(equity.pick[1,]),ST = XLF$Close[n],T=1,m=1000,n=100)

################### Problem C. Local volatility ###################
############ (a) #############
library(readxl)
setwd("E:/621 computational method/Final")
data <- read_excel("SPX.xls",na = "")
row1 <- names(data)
S0 <- as.numeric(row1[2])
r <- as.numeric(row1[3])/100
data <- data.frame(data)
data <- data[!is.na(data[,1]),]
SPX.data <- data[-1,2:4]
SPX.data<- SPX.data[!duplicated(SPX.data[,1:2]),]
colnames(SPX.data) <- c('muturity','strike','price')

T <- as.numeric(SPX.data[,1])
K <- as.numeric(SPX.data[,2])
MPrice <- as.numeric(SPX.data[,3])

# Function to calculate opntion price by BS formular
BS <- function(S0, K, tau, r, sigma,div){
  
  d1 <- (log(S0/K)+(r-div+sigma^2/2)*tau) / (sigma*sqrt(tau)) 
  d2 <- d1 - sigma*sqrt(tau)
  Price <- S0*exp(-div*tau)*pnorm(d1) - K*exp(-r*tau)*pnorm(d2)
  return(Price)
}

#Function to calculate error between BS and Market price
err <- function(S0,K,r,tau,sig,div,MPrice){
  BS(S0,K,tau,r,sig,div) - MPrice
}

# Function to find BS Implied Vol using Secant Method
secant <- function(S0,K,r,tau,div,MPrice){
  sig <- c()
  #loop for every strike price and market price
  for(i in 1:length(K)){
    x0 <- -1
    x1 <- 1
    err0 <- err(S0,K[i],r,tau[i],x0,div,MPrice[i])
    err1 <- err(S0,K[i],r,tau[i],x1,div,MPrice[i])
    #Loop until that the value of function to sigma is less than tolerance level 1e-4 
    while(abs(x1-x0) > 1e-4){
      x <- x1 - err1*(x1 - x0) / (err1 - err0) # Calculate the new x value
      x0 <- x1
      x1 <- x
    }
    sig <- c(x,sig)
  }
  return(sig)
}
Ivol <- secant(S0=S0,K=K,r=r,tau=T,div=0,MPrice=MPrice) #one month

t <- unique(T)
t1 <- rep(t[1],20)
t2 <- rep(t[2],20)
t3 <- rep(t[3],20)
t4 <- rep(t[4],20)
k1=k2=k3=k4=c()
for(i in 1:nrow(SPX.data)){
  if(T[i] == t[1]){
    k1 <- c(k1,SPX.data[i,2])
  }
  if(T[i] == t[2]){
    k2 <- c(k2,SPX.data[i,2])
  }
  if(T[i] == t[3]){
    k3 <- c(k3,SPX.data[i,2])
  }
  if(T[i] == t[4]){
    k4 <- c(k4,SPX.data[i,2])
  }
}
k1 <- sort(as.numeric(unique(k1)))[1:20]
k2 <- sort(as.numeric(unique(k2)))[33:52]
k3 <- sort(as.numeric(unique(k3)))[1:20]
k4 <- sort(as.numeric(unique(k4)))[1:20]
newdata <- (data.frame(K,T,Ivol))
Inv1=Inv2=Inv3=Inv4=c()
mprice1=mprice2=mprice3=mprice4=c()
for(i in 1:length(Ivol)){
  if(T[i] == t[1]){
    for(j in 1:20){
      if(K[i] == k1[j]){
        Inv1 = c(Inv1,Ivol[i])
        mprice1 = c(mprice1,SPX.data[i,3])
      }
    }
  }
  if(T[i] == t[2]){
    for(j in 1:20){
      if(K[i] == k2[j]){
        Inv2 = c(Inv2,Ivol[i])
        mprice2 = c(mprice2,SPX.data[i,3])
      }
    }
  }
  if(T[i] == t[3]){
    for(j in 1:20){
      if(K[i] == k3[j]){
        Inv3 = c(Inv3,Ivol[i])
        mprice3 = c(mprice3,SPX.data[i,3])
      }
    }
  }
  if(T[i] == t[4]){
    for(j in 1:20){
      if(K[i] == k4[j]){
        Inv4 = c(Inv4,Ivol[i])
        mprice4 = c(mprice4,SPX.data[i,3])
      }
    }
  }
}

# Creat 3D plot of volatilities versus K & T
library(scatterplot3d)
COLOR=c(rep("red",20),rep("lightgrey",20),rep("yellow",20),rep("lightblue",20))
scatterplot3d(c(k1,k2,k3,k4),c(t1,t2,t3,t4),c(Inv1,Inv2,Inv3,Inv4),main="Secant 3D plot",
              xlab="Strike Price",ylab = "Maturity",zlab = "Vol",color = COLOR,pch=16,type = "p")

############ (b) #############
library(rgl)
library(akima)

k.pick <- c(k1,k2,k3,k4)
t.pick <- c(t1,t2,t3,t4)
Inv <- c(Inv1,Inv2,Inv3,Inv4)
mprice <- as.numeric(c(mprice1,mprice2,mprice3,mprice4))
call <- data.frame(K,T,Ivol)

t.spline <- seq(min(t.pick),max(t.pick),0.01)
k.spline <- seq(min(k.pick),max(k.pick),1)

# By setting to use less points, it will "stretch" the surface over those points.
xyz <- with(call, interp(x=t.pick, y=k.pick, z=Inv, 
                         xo=t.spline, yo=k.spline, 
                         extrap=TRUE,linear = FALSE,duplicate = "mean" ))

with(xyz, persp3d(x,y,z, col=heat.colors(length(z))[rank(z)],ylim = c(min(k.pick),max(k.pick)),zlim = c(-2,2),xlab='maturity', 
                  ylab='strike', zlab='Implied volatility', main='IV Surface'))

############ (d) #############
#using implied volatility to calculate local volatility  
local.vol <- function(S0,K,tau,r,q,sig){
  dt=0.01
  dk=1
  dsig.dt = (xyz$z[length(xyz$z)/2+1]-xyz$z[length(xyz$z)/2])/(xyz$x[length(xyz$x)/2+1]-xyz$x[length(xyz$x)/2])
  dsig.dk = (xyz$z[length(xyz$z)/2+1]-xyz$z[length(xyz$z)/2])/(xyz$y[length(xyz$y)/2+1]-xyz$y[length(xyz$y)/2])
  dsig.dk2 = (xyz$z[length(xyz$z)/2]-xyz$z[length(xyz$z)/2])/(xyz$y[length(xyz$y)/2+1]-xyz$y[length(xyz$y)/2])^2
  d1 = (log(S0/K)+(r-q+1/2*sig^2)*tau)/sqrt(tau)
  term1 = 2*sig*dsig.dt*tau
  term2 = 2*sig*(r-q)*tau*K*dsig.dk
  dividend = term1+sig^2+term2
  term3 = (1 + K*d1*dsig.dk*sqrt(tau))^2
  term4 = K^2*tau*sig*(dsig.dk2-d1*dsig.dk2^2*sqrt(tau))
  divider = term3+term4
  x = sqrt(abs(dividend/divider))
  return(x)
}
lv <- local.vol(S0=S0,K=k.pick,tau = t.pick,r=r,q=0,sig = Inv)

xyz1 <- with(call, interp(x=t.pick, y=k.pick, z=lv, 
                         xo=t.spline, yo=k.spline, 
                         extrap=TRUE,linear = FALSE,duplicate = "mean" ))

with(xyz1, persp3d(x,y,z, col=heat.colors(length(z))[rank(z)],ylim = c(min(k.pick),max(k.pick)),zlim = c(-2,2),xlab='maturity', 
                  ylab='strike', zlab='Iocal volatility', main='IV Surface'))

############ (e) #############
price.local <- BS(S0=S0, K=k.pick, tau=t.pick, r=r, sigma=lv,div=0)

par(mfrow=c(1,2))
plot(price.local,ylab = "option price",main = "price by local vol")
plot(mprice,ylab = "option price",main = "market price")

############ (f) #############
table1 <- data.frame(t.pick,k.pick,mprice,Inv,lv,price.local)
colnames(table1) <- c('time to maturity','Strike price','market price','implied volatility','local volatility','price by localVol')
table1

price.bs <- BS(S0=S0, K=k.pick, tau=t.pick, r=r, sigma=Inv,div=0)
table2 <- data.frame(Inv,lv,price.bs,price.local)
colnames(table1) <- c('implied volatility','local volatility','Black Scholes price','price by localVol')
library(xlsx)
write.csv(table2, file="SPXvolitality.csv") 

############ (g) #############
data <- read_excel("JPM.xls",na = "")
row1 <- names(data)
S0 <- as.numeric(row1[2])
r <- as.numeric(row1[3])/100
data <- data.frame(data)
data <- data[!is.na(data[,1]),]
SPX.data <- data[-1,2:4]
SPX.data<- SPX.data[!duplicated(SPX.data[,1:2]),]
colnames(SPX.data) <- c('muturity','strike','price')
T <- as.numeric(SPX.data[,1])
K <- as.numeric(SPX.data[,2])
MPrice <- as.numeric(SPX.data[,3])
Ivol <- secant(S0=S0,K=K,r=r,tau=T,div=0,MPrice=MPrice) #one month
t <- unique(T)
t1 <- rep(t[1],20)
t2 <- rep(t[2],20)
t3 <- rep(t[3],20)
t4 <- rep(t[4],20)
k1=k2=k3=k4=c()
for(i in 1:nrow(SPX.data)){
  if(T[i] == t[1]){
    k1 <- c(k1,SPX.data[i,2])
  }
  if(T[i] == t[2]){
    k2 <- c(k2,SPX.data[i,2])
  }
  if(T[i] == t[3]){
    k3 <- c(k3,SPX.data[i,2])
  }
  if(T[i] == t[4]){
    k4 <- c(k4,SPX.data[i,2])
  }
}
k1 <- sort(as.numeric(unique(k1)))[1:20]
k2 <- sort(as.numeric(unique(k2)))[1:20]
k3 <- sort(as.numeric(unique(k3)))[1:20]
k4 <- sort(as.numeric(unique(k4)))[1:20]
newdata <- (data.frame(K,T,Ivol))
Inv1=Inv2=Inv3=Inv4=c()
mprice1=mprice2=mprice3=mprice4=c()
for(i in 1:length(Ivol)){
  if(T[i] == t[1]){
    for(j in 1:20){
      if(K[i] == k1[j]){
        Inv1 = c(Inv1,Ivol[i])
        mprice1 = c(mprice1,SPX.data[i,3])
      }
    }
  }
  if(T[i] == t[2]){
    for(j in 1:20){
      if(K[i] == k2[j]){
        Inv2 = c(Inv2,Ivol[i])
        mprice2 = c(mprice2,SPX.data[i,3])
      }
    }
  }
  if(T[i] == t[3]){
    for(j in 1:20){
      if(K[i] == k3[j]){
        Inv3 = c(Inv3,Ivol[i])
        mprice3 = c(mprice3,SPX.data[i,3])
      }
    }
  }
  if(T[i] == t[4]){
    for(j in 1:20){
      if(K[i] == k4[j]){
        Inv4 = c(Inv4,Ivol[i])
        mprice4 = c(mprice4,SPX.data[i,3])
      }
    }
  }
}

# Creat 3D plot of volatilities versus K & T
library(scatterplot3d)
COLOR=c(rep("red",20),rep("lightgrey",20),rep("yellow",20),rep("lightblue",20))
scatterplot3d(c(k1,k2,k3,k4),c(t1,t2,t3,t4),c(Inv1,Inv2,Inv3,Inv4),main="Secant 3D plot",
              xlab="Strike Price",ylab = "Maturity",zlab = "Vol",color = COLOR,pch=16,type = "p")

##(b)
k.pick <- c(k1,k2,k3,k4)
t.pick <- c(t1,t2,t3,t4)
Inv <- c(Inv1,Inv2,Inv3,Inv4)
mprice <- as.numeric(c(mprice1,mprice2,mprice3,mprice4))
call <- data.frame(K,T,Ivol)
t.spline <- seq(min(t.pick),max(t.pick),0.01)
k.spline <- seq(min(k.pick),max(k.pick),1)
xyz <- with(call, interp(x=t.pick, y=k.pick, z=Inv, 
                         xo=t.spline, yo=k.spline, 
                         extrap=TRUE,linear = FALSE,duplicate = "mean" ))
with(xyz, persp3d(x,y,z, col=heat.colors(length(z))[rank(z)],ylim = c(min(k.pick),max(k.pick)),zlim = c(-2,2),xlab='maturity', 
                  ylab='strike', zlab='Implied volatility', main='IV Surface'))

##(d)
lv <- local.vol(S0=S0,K=k.pick,tau = t.pick,r=r,q=0,sig = Inv)
xyz1 <- with(call, interp(x=t.pick, y=k.pick, z=lv, 
                          xo=t.spline, yo=k.spline, 
                          extrap=TRUE,linear = FALSE,duplicate = "mean" ))

with(xyz1, persp3d(x,y,z, col=heat.colors(length(z))[rank(z)],ylim = c(min(k.pick),max(k.pick)),zlim = c(-2,2),xlab='maturity', 
                   ylab='strike', zlab='Iocal volatility', main='IV Surface'))
##(e)
price.local <- BS(S0=S0, K=k.pick, tau=t.pick, r=r, sigma=lv,div=0)
##(f)
table1 <- data.frame(t.pick,k.pick,mprice,Inv,lv,price.local)
colnames(table1) <- c('time to maturity','Strike price','market price','implied volatility','local volatility','price by localVol')
table1

price.bs <- BS(S0=S0, K=k.pick, tau=t.pick, r=r, sigma=Inv,div=0)
table2 <- data.frame(Inv,lv,price.bs,price.local)
colnames(table1) <- c('implied volatility','local volatility','Black Scholes price','price by localVol')
library(xlsx)
write.xlsx(table2, file="JPMvolitality.xls") 