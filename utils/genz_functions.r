library("gtools")

### All the following arguments are dx1 vectors, where d is the dimension
### of the data. 

### E.g.
### x <- c(0.1, 0.1, 0.1)
### a <- c(50, 50, 50)
### u <- c(0.5, 0.5, 0.5)
### Since d (and sometimes x) are not passed in as arguments, to infer d, I need 
### arg a to be dx1

### continuous function ###

genz1 <- function(a, x, u){
  return(exp(-sum(a*sqrt((x-u)**2))))
}

genz1_int <- function(a, u){
  return(prod((1/a) * (2-exp(-a*u) - exp(a*(u-1)))))
}

### corner peak function ###
genz2 <- function(a, x, u){
  d <- length(x)
  return((1 + sum(a*x))**(-(d+1)))
}

genz2_int <- function(a, u) {
  d = length(a)
  normalizing = 1/ (factorial(d) * prod(a))
  
  acc = 0 
  acc <- acc + (-1)^d / (1 +  sum(a))
  for (k in 1:d){
    I_ = combinations(n=d, r=k)
    for (j in 1:nrow(I_)){
      acc <- acc + (-1)^(d+k) /  (1 +  sum(a) -  sum(a[I_[j,]]))
    }
  }
  return(normalizing * acc) 
}


### discontinuous ###

genz3 <- function(a, x, u){
  max_index <- min(2, length(x))
  if (sum(x[c(1: max_index)] > u[c(1: max_index)]) >= 1){
    return(0)
  } 
  return(exp(sum(a*x)))
}

genz3_int <- function(a, u){
  d = length(a)
  if (d > 1){
    return(prod((1/a) * (exp(min(1, u))-1)))
  }
  if (d == 1){
    return(1/a * exp(a*u) -1 )
  }
}


### gaussian peak ###

genz4 <- function(a, x, u){
  return(exp(-sum(a**2* (x-u)**2)))
}

genz4_int <- function(a, u){
  d <- length(a)
  return(pi**(d/2) * prod((1/a) *(pnorm(sqrt(2) *a*(1-u) - pnorm(-sqrt(2)*a*u)))))
}


### oscillatory ###
genz5_int <- function(a, u) {
  d = length(a)
  if (d %% 4 == 1){h_d = function(x){sin(x)}}
  if (d %% 4 == 2){h_d = function(x){-cos(x)}}
  if (d %% 4 == 3){h_d = function(x){-sin(x)}}
  if (d %% 4 == 0){h_d = function(x){cos(x)}}
  
  normalizing = 1/prod(a) 
  acc = 0 
  acc <- acc + h_d(2*pi*u[1] + sum(a))
  for (k in 1:d){
    I_ = combinations(n=d, r=k)
    for (j in 1:nrow(I_)){
      acc <- acc + (-1)^k * h_d(2 * pi * u[1] + sum(a) - sum(a[I_[j,]]))
    }
  }
  return(normalizing * acc) 
}





### product peak ###

genz6 <- function(a, x, u){
  return(prod(a**(-2) + (x-u)**2)**(-1))
}

genz6_int <- function(a, u){
  return(prod(a * (atan(a*(1-u)) - atan(-a*u))))
}


### step function ### 

step <- function(a, x, u){
  return(as.numeric(x[1] > 0.5))  
}

step_int <- function(a,u){
  return(.5)
}
