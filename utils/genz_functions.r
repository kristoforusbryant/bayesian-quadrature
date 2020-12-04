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


genz2 <- function(a, x, u){
  d <- length(x)
  return((1 + sum(a*x))**(-(d+1)))

genz1_int <- function(a, u){
  return(prod((1/a) * 2-exp(-a*u) - exp(a*(u-1))))
}


### corner peak function ###

}

genz2_int <- function(a){
  d <- length(a)
  acc <- (-1)**d * (1 + sum(a))**(-1)
  
  for (k in c(1:d)){
    I_ <- combinations(n=d, r=k)
    
    for (j in c(1:length(I_))){
      inner_sum <- (-1)**(k+d) * (1 + sum(a) - sum(a[I_[j,]]))**(-1)
      acc <- acc + inner_sum
    }
    return((1/(prod(a) * factorial(d))) * acc)
  }
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
  return(pi**(d/2) * prod((1/a) *pnorm(sqrt(2) *a*(1-u) - pnorm(-sqrt(2)*a*u))))
}


### oscillatory ###

genz5 <- function(a, x, u){
  return(cos(2*pi*u[1] + sum(a*x)))
}

genz5_int <- function(a, u){
  d <- length(a)
  negative_sin <- function(x){-sin(x)}
  negative_cos <- function(x){-cos(x)}
  
  f_list <- list(sin, negative_cos, negative_sin, cos)
  chosen_f <- f_list[[d%%4 + 4*(d%%4==0)]]
  
  acc <- chosen_f(2*pi*u[1] + sum(a))
  for (k in c(1:d)){
    I_ <- combinations(n=d, r=k)
    for (j in c(1:length(I_))){
      inner_sum <- (-1)**(k) * chosen_f(2*pi*u[1] + sum(a) - sum(a[I_[j,]]))
      acc <- acc + inner_sum}
  }
  return((1/prod(a))* acc)
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
