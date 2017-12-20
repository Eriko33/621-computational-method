######################################## Problem1 ##############################################
####################(1)######################
library(Sim.DiffProc)
sample_data <- read.csv("sample_data.csv")
plot(sample_data$stock1,type = "l",col = "red",xlim = c(1,8000),ylim = c(50,300))
lines(sample_data$stock2,col = "blue")
lines(sample_data$stock3,col = "green")
lines(sample_data$stock4,col = "gray")
lines(sample_data$stock5)

model.select <- function(data){
  #model 1
  fx1 <- expression(theta[1]*x)
  gx1 <- expression(theta[2]*x^theta[3])
  mod1 <- fitsde(data=as.ts(data),drift=fx1,diffusion=gx1,
                 start = list(theta1=1,theta2=1,theta3=1),pmle="ozaki")
  AIC(mod1)
  #model2
  fx2 <- expression(theta[1]+theta[2]*x)
  gx2 <- expression(theta[3]*x^theta[4])
  mod2 <- fitsde(data=as.ts(data),drift=fx2,diffusion=gx2,
                 start = list(theta1=1,theta2=1,theta3=1,theta4=1),pmle="ozaki")
  AIC(mod2)
  #model3
  fx3 <- expression(theta[1]+theta[2]*x)
  gx3 <- expression(theta[3]*sqrt(x))
  mod3 <- fitsde(data=as.ts(data),drift=fx3,diffusion=gx3,
                 start = list(theta1=1,theta2=1,theta3=1),pmle="ozaki")
  
  #model4
  fx4 <- expression(theta[1])
  gx4 <- expression(theta[2]*x^theta[3])
  mod4 <- fitsde(data=as.ts(data),drift=fx4,diffusion=gx4,
                 start = list(theta1=1, theta2=1,theta3=1),pmle="ozaki")
  
  #model5
  fx5 <- expression(theta[1]*x)
  gx5 <- expression(theta[2]+theta[3]*x^theta[4])
  mod5 <- fitsde(data=as.ts(data),drift=fx5,diffusion=gx5,
                 start = list(theta1=1, theta2=1,theta3=1,theta4=1),pmle="ozaki")
  
  ## Computes AIC
  AIC <- c(AIC(mod1),AIC(mod2),AIC(mod3),AIC(mod4),AIC(mod5))
  Test <- data.frame(AIC,row.names = c("mod1","mod2","mod3","mod4","mod5"))
  Bestmod <- rownames(Test)[which.min(Test[,1])]
  cat("best model: ",Bestmod,"\n")
  print(Test)
}
col1.AIC <- model.select(sample_data$stock1)  #mod3
col2.AIC <- model.select(sample_data$stock2)  #mod1
col3.AIC <- model.select(sample_data$stock3)  #mod4
col4.AIC <- model.select(sample_data$stock4)  #mod3
col5.AIC <- model.select(sample_data$stock5)  #mod1


##################(2)####################
pmle <- eval(formals(fitsde.default)$pmle)
#model1
fx1 <- expression(theta[1]*x)
gx1 <- expression(theta[2]*x^theta[3])
fitres1 <- lapply(1:4, function(i) fitsde(as.ts(sample_data$stock2),drift=fx1,diffusion=gx1,pmle=pmle[i],
                                          start = list(theta1=1,theta2=1,theta3=1)))
Coef.model1 <- data.frame(do.call("cbind",lapply(1:4,function(i) coef(fitres1[[i]]))))
Info1 <- data.frame(do.call("rbind",lapply(1:4,function(i) logLik(fitres1[[i]]))),
                  do.call("rbind",lapply(1:4,function(i) AIC(fitres1[[i]]))),
                  do.call("rbind",lapply(1:4,function(i) BIC(fitres1[[i]]))),
                  row.names=pmle)
names(Coef.model1) <- c(pmle)
names(Info1) <- c("logLik","AIC","BIC")
Coef.model1
Info1

f <- expression(0.00001825*x)
g <- expression(0.007829*x^0.5142)
mod <- snssde1d(drift=f,diffusion=g,x0=sample_data$stock2[1],M=100, N=length(sample_data$stock2),t0 = 1,T = 100001)
mod

plot(time(mod)[-1],mean(mod)[-1],type = "l", col = "red",lwd = 2, xlab = "time", ylab = "stockprice")
lines(time(mod)[-1],sample_data$stock2,col = "blue",lwd = 2)
legend("topleft",c("model mean path","real data"),inset = .01,col=c("red","blue"),lwd=2,cex=0.7)

#model3
fx3 <- expression(theta[1]+theta[2]*x)
gx3 <- expression(theta[3]*sqrt(x))
fitres3 <- lapply(1:4, function(i) fitsde(as.ts(sample_data$stock1),drift=fx3,diffusion=gx3,pmle=pmle[i],
                  start = list(theta1=1,theta2=1,theta3=1)))
Coef.model3 <- data.frame(do.call("cbind",lapply(1:4,function(i) coef(fitres3[[i]]))))
Info2 <- data.frame(do.call("rbind",lapply(1:4,function(i) logLik(fitres3[[i]]))),
                    do.call("rbind",lapply(1:4,function(i) AIC(fitres3[[i]]))),
                    do.call("rbind",lapply(1:4,function(i) BIC(fitres3[[i]]))),
                    row.names=pmle)
names(Coef.model3) <- c(pmle)
names(Info2) <- c("logLik","AIC","BIC")
Coef.model3
Info2

f <- expression( (0.0002911+0.00002188*x) )
g <- expression( 0.004014*sqrt(x) )
mod <- snssde1d(drift=f,diffusion=g,x0=sample_data$stock1[1],M=100, N=length(sample_data$stock1),t0 = 1,T = 100001)
mod

plot(time(mod)[-1],mean(mod)[-1],type = "l", col = "red",lwd = 2, xlab = "time", ylab = "stockprice")
lines(time(mod)[-1],sample_data$stock1,col = "blue",lwd = 2)
legend("topleft",c("model mean path","real data"),inset = .01,col=c("red","blue"),lwd=2,cex=0.7)

#model4
fx4 <- expression(theta[1])
gx4 <- expression(theta[2]*x^theta[3])
fitres4 <- lapply(1:4, function(i) fitsde(as.ts(sample_data$stock3),drift=fx4,diffusion=gx4,pmle=pmle[i],
                                          start = list(theta1=1,theta2=1,theta3=1)))
Coef.model4 <- data.frame(do.call("cbind",lapply(1:4,function(i) coef(fitres4[[i]]))))
Info3 <- data.frame(do.call("rbind",lapply(1:4,function(i) logLik(fitres4[[i]]))),
                    do.call("rbind",lapply(1:4,function(i) AIC(fitres4[[i]]))),
                    do.call("rbind",lapply(1:4,function(i) BIC(fitres4[[i]]))),
                    row.names=pmle)
names(Coef.model4) <- c(pmle)
names(Info3) <- c("logLik","AIC","BIC")
Coef.model4
Info3

f <- expression(1)
g <- expression(1*x^1)
mod <- snssde1d(drift=f,diffusion=g,x0=sample_data$stock3[1],M=100, N=length(sample_data$stock3),t0 = 1, T = 100001)
mod

plot(time(mod)[-1],mean(mod)[-1],type = "l", col = "red",lwd = 2, xlab = "time", ylab = "stockprice")
lines(time(mod)[-1],sample_data$stock3,col = "blue",lwd = 2)
legend("topleft",c("model mean path","real data"),inset = .01,col=c("red","blue"),lwd=2,cex=0.7)


######################################## Problem2 #############################################
####################(1)#######################
setwd("E:/621 computational method/H4")
data <- read.csv("2017_2_15_mid.csv")

#choose five year col to finish the question
num<-seq(from=1,to=nrow(data),by=2)
sigma.mkt <- data$X5Yr[num]/100
strike.mkt <- data$X5Yr[-num]/100

beta <- 0.5 #beta parameter
F <- strike.mkt #spot price equals to strike price because at the money
K <- strike.mkt
T <- 5

sigma.ATM <- function(x){ 
  
  A = (1-beta)^2/24*x[1]^2/F^(2-2*beta)
  B = 1/4*x[2]*beta*x[3]*x[1]/F^(1-beta)
  C = (2-3*x[2]^2)/24*x[3]^2
  term1 = x[1]*(1 + (A + B + C)*T)
  term2 = F^(1-beta)
  
  return(term1/term2)
}

SSE <- function(x){
  return(sum((sigma.mkt-sigma.ATM(x))^2))
}

library(pracma)
x0 <- c(1,0,0)
opt1 <- optim(x0,SSE)
Alpha1 <- opt1$par[1]
Rho1 <- opt1$par[2]
nu1 <- opt1$par[3] #vol of vol
data.frame(Alpha1,Rho1,nu1)


#########################(2)##########################
beta <- 0.7
opt2 <- optim(x0,SSE)
Alpha2 <- opt2$par[1]
Rho2 <- opt2$par[2]
nu2 <- opt2$par[3] #vol of vol
data.frame(Alpha2,Rho2,nu2)

beta <- 0.4
opt3 <- optim(x0,SSE)
Alpha3 <- opt3$par[1]
Rho3 <- opt3$par[2]
nu3 <- opt3$par[3] #vol of vol
data.frame(Alpha3,Rho3,nu3)

#########################(3)##########################
vol1 <- sigma.ATM(opt1$par)
plot(vol1,type = "l",ylim = c(0,3))
vol1
vol2 <- sigma.ATM(opt2$par)
lines(vol2,col = "red")

vol3 <- sigma.ATM(opt3$par)
lines(vol3,col = "blue")
legend("topleft",c("beta=0.5","beta=0.7","beta=0.4"),inset = .01,col=c("black","red","blue"),lwd=2,cex=0.75)
#########################(4)##########################
opt1$value
opt2$value
opt3$value

#########################(5)##########################
sigma.mkt <- data$X20Yr[num]/100
strike.mkt <- data$X20Yr[-num]/100
beta <- 0.7 #beta parameter
F <- strike.mkt #spot price equals to strike price because at the money
K <- strike.mkt
T <- 20

x <- c(Alpha2,Rho2,nu2)
sigma.ATM(x)
sigma.mkt
mean(abs(sigma.ATM(x)-sigma.mkt))
sum((sigma.ATM(x)-sigma.mkt)^2)
############################################## Problem 3 #############################################
i <- sqrt(-1+0i)
S0 <- 1
r <- 0.05
T <- 30/252
sigma <- 0.2
k <- log(K)

alpha <- 3
lambda <- 0.001
N <- 1024
eta <- 2*pi/N/lambda

chara_fun <- function(u){
  exp(i*(log(S0)+(r-sigma^2/2)*T)*u - 0.5*sigma^2*T*u^2)
}

psi_fun <-function(nu){
  term1 = exp(-r*T)*chara_fun(nu-(alpha+1)*i)
  term2 = alpha^2+alpha-nu^2+i*(2*alpha+1)*nu
  term1/term2
}

Kronecker.delta <- function(x){
  if(x==0)
    return(1)
  else
    return(0)
}

call_price <- function(u){
  #substitute some parameters
  a = N*eta
  b = 1/2*N*lambda #strike level from -b to b
  ku = -b+lambda*(u-1)

  term = c()
  #Using the trapezoid rule for the integral
  for(j in 1:N){
    nu.j = eta*(j-1)
    term1 = exp(-i*2*pi/N*(j-1)*(u-1)) * exp(i*b*nu.j)
    term2 = psi_fun(nu.j) * eta/3 * (3+(-1)^j-Kronecker.delta(j-1))
    term = c(term,term1*term2)
  }
  return(Re(exp(-alpha*ku)/pi * sum(term)))
}
call_price(513)

#calculate European call option price by BS model
library(RQuantLib)
EuropeanOption("call", underlying = 1, strike = 1, dividendYield = 0, riskFreeRate = 0.05,
               maturity = 30/252, volatility = 0.2)$value

########################### Bonus ##############################
S0 <- 1
K <- c(0.5,0.75,1,1.25,1.5)
r <- 0
T <- 5
div <- 0
sigma.nu <- 0.2
k <- log(K)

alpha <- 2
lambda <- 0.01
N <- 1024
eta <- 2*pi/N/lambda

mu <- 0.1
nu0 <- 0.1
kappa <- (4+2+1)/3
rho <- -0.3

chara_fun_heston <- function(u){
  zeta <- -(u^2+i*u)/2
  gamma <- kappa - i*rho*sigma.nu*u
  theta <- sqrt(gamma^2-2*sigma.nu^2*zeta)
  A = i*u*(log(S0)+(r-div)*T)
  B = 2*zeta*(1-exp(-theta*T))/(2*theta-(theta-gamma)*(1-exp(-theta*T)))
  C = 2*log((2*theta-(theta-gamma)*(1-exp(-theta*T)))/(2*theta))
  E = (theta-gamma)*T
  H = kappa*mu/sigma.nu^2*(C + E) 
  exp(A + B*nu0 - H)
}

psi_fun_heston <-function(nu){
  term1 = exp(-r*T)*chara_fun_heston(nu-(alpha+1)*i)
  term2 = alpha^2+alpha-nu^2+i*(2*alpha+1)*nu
  term1/term2
}

call_price_heston <- function(u){
  #substitute some parameters
  a = N*eta
  b = N*lambda/2 #strike level from -b to b
  ku = -b+lambda*(u-1)
  #u =-(k+b)/lambda+1

  term <- c()
  #Using the trapezoid rule for the integral
  for(j in 1:N){
    nu.j = eta*(j-1)
    term1 = exp(-i*2*pi/N*(j-1)*(u-1)) * exp(i*b*nu.j)
    term2 = psi_fun_heston(nu.j) * eta/3 * (3+(-1)^j-Kronecker.delta(j-1))
    term = c(term,term1*term2)  
  }
  return(Re(exp(-alpha*ku)/pi * sum(term)))
}
b = N*lambda/2 #strike level from -b to b
u <- (k+b)/lambda+1
for(a in 1:5){
  heston.FFT <- call_price_heston(u[a])
}
