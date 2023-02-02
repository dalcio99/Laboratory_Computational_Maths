# Solving 2D SIR stochastic differential equation
library(tidyverse)
library(deSolve)
library(Sim.DiffProc)
set.seed(1234, kind = "L'Ecuyer-CMRG")

# EXAMPLE 1: Disease-free equilibrium
#mu = 1; sigma = 0.5; gam = 0.8; bet = 2; lambda = 1; epsilon = 0.2; t_1 = 2; t_2 = 7

# EXAMPLE 2: Persistence equilibrium
mu=1; sigma=0.5; gam=0.8; bet = 3; lambda = 1 ;epsilon = 0.2; t_1 = 2; t_2 = 7

x0 = 0.8; y0 = 0.19   # Initial values

### STOCHASTIC MODEL ###
f <- expression((lambda +(-bet)*x*y-mu*x),((bet)*x*y-(mu+gam+epsilon)*y))  
g <- expression((-sigma*x*y),(sigma*x*y))
#correlation <- c(2,0,0,2)
#dim(correlation) <- c(2,2)
mod2d <- snssde2d(N= 1000, M=100, x0=c(x0,y0), t0=0, T=10, drift=f, diffusion=g, method="euler")

mod2d


### DETERMINISTIC MODEL ###
ode.sir.model <- function (t, x, parms) {
  ## first extract the state variables
  S <- x[1]
  I <- x[2]
  R <- x[3]
  mu <- parms[1]
  gam <- parms[2]
  bet <- parms[3]
  lambda <- parms[4]
  epsilon <- parms[5]
  ## now code the model equations
  dSdt <- lambda - mu*S - bet*S*I
  dIdt <- bet*S*I - (mu+gam+epsilon)*I
  dRdt <- gam*I - mu*R
  ## combine results into a single vector
  dxdt <- c(dSdt,dIdt,dRdt)
  ## return result as a list!
  list(dxdt)
}

## ode in time
times <- seq(from=0,to=10,by=0.001)         # mesh
xstart <- c(S=x0,I=y0,R=(1-x0-y0))         # initial values
parms <- c(mu, gam, bet, lambda, epsilon)  # parameters

ode(func=ode.sir.model, y=xstart, times=times, parms = parms) %>% as.data.frame() -> sir_ode

## Parameters defined in the paper and mentioned in stability theorems
R0 = bet*lambda/(mu*(mu+gam+epsilon))
stoc_R0 = R0-(sigma*lambda)**2/(2*mu**2*(mu+gam+epsilon))
var = sigma**2
discrim1 = mu*bet/lambda
discrim2 = bet**2/(2*(mu+gam+epsilon))

summary(mod2d, at = time)

## sde in time
plot(mod2d, col = c("blue","red"))

## sde model in plane (O,X,Y)
plot2d(mod2d,type="n") 
points2d(mod2d,col=rgb(0,100,0,50,maxColorValue=255), pch=16)

out <- rsde2d(object = mod2d, at = t_1)

head(out,n=3)

## sde model marginal density
denM <- dsde2d(mod2d,pdf="M",at = t_1)
plot(denM, main="Marginal Density", col=c(rgb(0,0,255,70,maxColorValue=255), rgb(255,0,0,70,maxColorValue=255)))


## sde model Joint density
denJ <- dsde2d(mod2d, pdf="J", n=100,at = t_1)
plot(denJ,display="contour",main="Bivariate Transition Density at time t=2")

plot(denJ,display="persp",main="Bivariate Transition Density at time t=2")

### Mean vs Deterministic solutions ###

# dev.new()
plot(mod2d,type="n")
mx <- apply(mod2d$X,1,mean)
my <- apply(mod2d$Y,1,mean)
lines(sir_ode$time, sir_ode$S, col="orange", lwd=2)
lines(sir_ode$time, sir_ode$I, col="green", lwd=2)
lines(time(mod2d),mx,col="blue",lwd=2)
lines(time(mod2d),my,col="red", lwd=2)
#lines(time(mod2d),mz,col="green")
#legend("topright",c(expression(E(S[t])),expression(E(I[t])), expression(S(t)),expression(I(t))),lty=1,inset = .01,col=c("blue", "red", "orange", "green"),cex=0.95)

#############################################################

### Convergence analysis ###
N = 1000; T = 10; t0 = 0; M = 100 # <- NB: these have to match with those at line 20 to work out

# "Discrepancy" on S(t)
S_mean_time_t <- vector(mode = "numeric", length = (N+1) )
RandS <- mod2d$X

for (n in 1:N+1) {
  errors <- vector(mode = "numeric", length = M )
  x = (10000*n) %/% (N+1) # the factor 10**4 follows from the fixed mesh with dt=0.001 at line 47 
  for (m in 0:M-1) {
    errors[m+1] <- abs(RandS[m*(N+1)+n] - sir_ode$S[x]) 
  }
  S_mean_time_t[n] <- mean(errors)
} 
max(S_mean_time_t)

# "Discrepancy" on I(t)
I_mean_time_t <- vector(mode = "numeric", length = (N+1) )
RandI <- mod2d$Y

for (n in 1:N+1) {
  x = (10000*n) %/% (N+1)
  errors <- vector(mode = "numeric", length = M )
  for (m in 0:M-1) {
    errors[m+1] <- abs(RandI[m*(N+1)+n] - sir_ode$I[x]) 
  }
  I_mean_time_t[n] <- mean(errors)
} 
max(I_mean_time_t)

