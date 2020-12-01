library(matlib)
rho = 3
# Kernels
l2_norm <- function(x){ sum(x * x)}
kernel_Matern <- function(x,y){
  return((1 + (sqrt(3) * l2_norm(x - y) / rho)) * exp(-(sqrt(3) * l2_norm(x - y) / rho)))
}
kernel_add <- function(x,y){sum(x+y)} # trivial for testing


# Kernel multiplication
#k_mult = function(X, X_, kernel){
  # X and X_ are m x d and n x d matrices
#  if (is.null(nrow(X))){ X = t(as.matrix(X))}
#  if (is.null(nrow(X_))){ X_ = t(as.matrix(X_))}
#  K = matrix(rep(0, len=nrow(X)*nrow(X_)), nrow = nrow(X))
#  for (i in 1:nrow(X)){
#    K[i,] = apply(X_, 1, function(a){kernel(X[i,],a)}) 
#  }
#  return(K)
#}

k_mult = function(X, X_, kernel){
  # X and X_ are m x d and n x d matrices
  if (is.null(nrow(X))){ X = t(as.matrix(X))}
  if (is.null(nrow(X_))){ X_ = t(as.matrix(X_))}
  K = apply(X_, 1, function(a){apply(X, 1, function(b){kernel(a,b)})})
  return(t(as.matrix(K)))
}


# GP 
GP <- setRefClass("GP", 
  fields= list(
    d = "numeric",
    X = "matrix", 
    y = "numeric",
    kernel = "function",
    C_N = "matrix", 
    C_N_inv = "matrix", 
    sigmasq = "numeric"
  ), 
  methods =  list(
    initialize = function(d, kernel, sigmasq){
      .self$d = d
      .self$X = matrix(nrow=0, ncol=d)
      .self$kernel = kernel
      .self$sigmasq = sigmasq
    },
    mean = function(x){k_mult(x, X, kernel) %*% C_N_inv %*% y},
    cov = function(x){ 
      k = k_mult(x, X, kernel)
      return(as.numeric(kernel(x,x) + sigmasq - k %*% C_N_inv %*% t(k)))
    },
    seq_update = function(X_new,y_new){
      if (!sum(dim(C_N))){
        X <<- rbind(X, X_new) 
        y <<- c(y, y_new)
        C_N <<- k_mult(X_new, X_new, kernel) + sigmasq * diag(nrow(X_new))
        C_N_inv <<- inv(C_N)
      }
      else{
        k = k_mult(X, X_new, kernel)
        c = k_mult(X_new,X_new, kernel) + sigmasq * diag(nrow(X_new))
        C_N <<- cbind(C_N, k) 
        C_N <<- rbind(C_N, cbind(t(k), c)) 
        C_N_inv <<- inv(C_N)
        X <<- rbind(X, X_new) 
        y <<- c(y, y_new)  
      }
    }, 
    kernel_mean_MI = function(MI_iter=10000, min=0, max=1){ # Pi(k(.,X)) kernel mean 
      return(rowSums(k_mult(matrix(runif(d*MI_iter, 0,1), ncol=d), X, kernel)) / MI_iter)
    },
    full_kernel_integral_MI = function(MI_iter=10000, min=0, max=1){
      sample = matrix(runif(d*MI_iter*2, 0,1), ncol=d*2)
      return(sum(apply(sample, 1, function(a){kernel(a[1:d], a[(d+1):(d*2)])}))/MI_iter)
    },
    GP_BQ = function(MI_iter=10000, min=0, max=1){
      integral = list()
      k_mean = t(as.matrix(kernel_mean_MI(MI_iter=MI_iter, min=0, max=1)))
      full_k_int = full_kernel_integral_MI(MI_iter=MI_iter, min=0, max=1)
      integral$E = as.numeric(k_mean %*% C_N_inv %*% y)
      integral$V = full_k_int - as.numeric(k_mean %*% C_N_inv %*% t(k_mean))
      return(integral)
    }
    )
)







