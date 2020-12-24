#setwd('Desktop/bayesian-quadrature/utils/')
source("bart_int.R")


# hyperparameters to test the 

d = 1
a = rep(150/(d^3), d)

N = 50*d
x = matrix(runif(N * d), ncol = d) # N samples of d dimensions
u = rep(.5, d)


y = apply(x, 1, function(x_){step(a, x_, u)}) 
n_total = 50*d + 20*d 


test_run <- run_iterations(integral_function=step, genz_func_name="step_10", a, n_total, n_tree=200, S = 5000, d=10, m=1000, n_iter=1)

