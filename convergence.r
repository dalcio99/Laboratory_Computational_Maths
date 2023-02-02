# Solving 2D SIR stochastic differential equation
library(tidyverse)
library(deSolve)
library(Sim.DiffProc)

set.seed(1234, kind = "L'Ecuyer-CMRG")
x0=0.8; y0=0.19

# EXAMPLE 1: Disease-free equilibrium
#mu=1; sigma=0.5; gam=0.8; bet = 2; lambda = 1 ;epsilon = 0.2; t_1=2; t_2=7 ; T = 10; t_0=0
# EXAMPLE 2: Persistence equilibrium
 mu=1; sigma=0.5; gam=0.8; bet = 3; lambda = 1 ;epsilon = 0.2; t_1=2; t_2=7 ; T = 10; t_0=0

# sde equations
f <- expression((lambda +(-bet)*x*y-mu*x),((bet)*x*y-(mu+gam+epsilon)*y))  
g <- expression((-sigma*x*y),(sigma*x*y))

# ode equations
ode.sir.model <- function (t, x, parms) {  # vector field that defines eqns
  S <- x[1]
  I <- x[2]
  R <- x[3]
  mu <- parms[1]
  gam <- parms[2]
  bet <- parms[3]
  lambda <- parms[4]
  epsilon <- parms[5]
  ##model equations
  dSdt <- lambda - mu*S - bet*S*I
  dIdt <- bet*S*I - (mu+gam+epsilon)*I
  dRdt <- gam*I - mu*R
  ## combined into a single vector
  dxdt <- c(dSdt,dIdt,dRdt)
  ## return result as a list!
  list(dxdt)
}

## Parameters defined in the paper and mentioned in stability theorems

R0=bet*lambda/(mu*(mu+gam+epsilon))-(sigma*lambda)**2/(2*mu**2*(mu+gam+epsilon))
var = sigma**2
discrim1 = mu*bet/lambda
discrim2 = bet**2/(2*(mu+gam+epsilon))

# ode model
#Dt = (T-t_0)/n1
Dt = 0.001
times <- seq(from=0,to=10,by=Dt)         # mesh
xstart <- c(S=x0,I=y0,R=(1-x0-y0))         # initial values
parms <- c(mu, gam, bet, lambda, epsilon)  # parameters to be passed to func

ode(func=ode.sir.model, y=xstart, times=times, parms = parms) %>%
  as.data.frame() -> out  # solutions are saved in the variable "out"

#errors
errors_S <- list()
errors_I <- list()
elapsed_times <- list() # time to run snssde2d
total_times <- list()   # total time (computing the error)


for (j in 1:1){
  m1 = 100*j
  errors_S_j <- list()
  errors_I_j <- list()
  elapsed_times_j <- list()
  total_times_j <- list()
  #correlation <- c(2,1,1,2)
  #dim(correlation) <- c(2,2)
  for (n1 in seq(from=100,to=1000,by=100)) {
    start = Sys.time()
    mod2d <- snssde2d(N=n1, M=m1, x0=c(x0,y0), t0=0, T=10, drift=f, diffusion=g,method="euler")
    
    ### time complexity ###
    epsime = as.numeric(difftime(Sys.time(), start, units=c("secs")))
    elapsed_times_j = append(elapsed_times_j, epsime)
    
    ### Convergence analysis ###
    N = n1; M = m1
    # Error on S(t)
    S_mean_time_t <- vector(mode = "numeric", length = (N+1) )
    RandS <- mod2d$X
    for (n in 1:N+1) {
      x = (10000*n) %/% (N+1)
      errors <- vector(mode = "numeric", length = M )
      for (m in 0:M-1) {
        errors[m+1] <- abs(RandS[m*(N+1)+n] - out$S[x]) 
      }
      S_mean_time_t[n] <- mean(errors)
    } 
    errors_S_j = append(errors_S_j,max(S_mean_time_t))
    
    # Error on I(t)
    I_mean_time_t <- vector(mode = "numeric", length = (N+1) )
    RandI <- mod2d$Y
    for (n in 1:N+1) {
      x = (10000*n) %/% (N+1) # the factor 10**4 follows from the fixed mesh with dt=0.001 at line 47 
      errors <- vector(mode = "numeric", length = M )
      for (m in 0:M-1) {
        errors[m+1] <- abs(RandI[m*(N+1)+n] - out$I[x]) 
      }
      I_mean_time_t[n] <- mean(errors)
    } 
    tot_time = as.numeric(difftime(Sys.time(), start, units=c("secs")))
    total_times_j = append(total_times_j, tot_time)
    errors_I_j = append(errors_I_j,max(I_mean_time_t))
  }
  errors_S = append(errors_S,errors_S_j)
  errors_I = append(errors_I,errors_I_j)
  elapsed_times = append(elapsed_times,elapsed_times_j)
  total_times = append(total_times_j,total_times)
}
errors_S
errors_I
elapsed_times
total_times
numbers <- seq(from=100,to=1000,by=100)
plot(numbers,errors_S, type="b", xlab = "Mesh dimension: Value of N", ylab = "Discrepancy on S")
plot(numbers,errors_I,type="b", xlab = "Mesh dimension: Value of N", ylab = "Discrepancy on I")
plot(numbers,elapsed_times,type="b", xlab = "Mesh dimension: Value of N", ylab = "E-M time complexity (sec)")
plot(numbers,total_times,type="b", xlab = "Mesh dimension: Value of N", ylab = "Program runtime (sec)")

