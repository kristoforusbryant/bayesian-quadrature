library(matlib)
# Kernels
l2_norm <- function(x){ sqrt(sum(x * x))}
kernel_Matern <- function(x,y, lengthscale){
  return((1 + (sqrt(3) * l2_norm(x - y) / lengthscale)) * exp(-(sqrt(3) * l2_norm(x - y) / lengthscale)))
}
kernel_rbf <- function(x,y, l=1){
  return(exp(- (l2_norm(x - y))/ (2*l*l)))
}

kernel_add <- function(x,y){sum(x+y)} # trivial for testing

# Kernel multiplication
k_mult = function(X, X_, kernel, lengthscale){
  # X and X_ are m x d and n x d matrices
  if (is.null(nrow(X))){ X = t(as.matrix(X))}
  if (is.null(nrow(X_))){ X_ = t(as.matrix(X_))}
  K = apply(X_, 1, function(a){apply(X, 1, function(b){kernel(a,b, lengthscale)})})
  if (is.null(nrow(K))){K = t(as.matrix(K))}
  return(K) # returns m x n matrices
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
    sigmasq = "numeric", 
    lengthscale = "numeric"
  ), 
  methods =  list(
    initialize = function(d, kernel, sigmasq, lengthscale){
      .self$d = d
      .self$X = matrix(nrow=0, ncol=d)
      .self$kernel = kernel
      .self$sigmasq = sigmasq
      .self$lengthscale = lengthscale
    },
    mean = function(x){k_mult(x, X, kernel, lengthscale) %*% C_N_inv %*% as.matrix(y)},
    cov = function(x){ 
      k = k_mult(x, X, kernel, lengthscale) 
      return(as.numeric(kernel(x,x, lengthscale) + sigmasq - k %*% C_N_inv %*% t(k)))
    },
    seq_update = function(X_new,y_new){
      if (!sum(dim(C_N))){
        X <<- rbind(X, X_new) 
        y <<- c(y, y_new)  
        C_N <<- k_mult(X_new, X_new, kernel, lengthscale) + sigmasq * diag(nrow(X_new))
        if (all(dim(X_new) == 1)){C_N_inv <<- 1/C_N}
        else{C_N_inv <<- solve(C_N)}
      }
      else{
        k = k_mult(X, X_new, kernel, lengthscale) # N x 1
        c = k_mult(X_new, X_new, kernel, lengthscale) + sigmasq * diag(nrow(X_new))
        C_N <<- cbind(C_N, k) 
        C_N <<- rbind(C_N, cbind(t(k), c)) 
        C_N_inv <<- solve(C_N)
        X <<- rbind(X, X_new) 
        y <<- c(y, y_new)  
      }
    }, 
    kernel_mean_MI = function(MI_iter=10000, min=0, max=1){ # Pi(k(.,X)) kernel mean 
      k_mean = rowSums(k_mult(X, matrix(runif(d*MI_iter, 0,1), ncol=d), kernel, lengthscale)) / MI_iter
      if (is.null(nrow(k_mean))){k_mean = t(as.matrix(k_mean))}
      return(k_mean)
    },
    full_kernel_integral_MI = function(MI_iter=10000, min=0, max=1){
      sample = matrix(runif(d*MI_iter*2, 0,1), ncol=d*2)
      return(sum(apply(sample, 1, function(a){kernel(a[1:d], a[(d+1):(d*2)], 
                                                     lengthscale)}))/MI_iter)
    },
    log_marginal_lik = function(){
      data_fit = -1/2 * t(as.matrix(y)) %*% C_N_inv %*% as.matrix(y) 
      complexity_penalty = -1/2 * as.numeric(
        determinant.matrix(C_N, logarithm = TRUE)$modulus)
      return(data_fit + complexity_penalty - length(y)/2 * log(2 * pi))
    },
    seq_design = function(n_of_pts=20*d, S=10000, dist_f=runif, f){
      for (i in 1:n_of_pts){
        c = matrix(runif(S * d), ncol = d)
        m = apply(c, 1, function(x_){.self$cov(x_)} )
        c_ = t(as.matrix(c[which.max(m),]))
        seq_update(c_, as.matrix(f(c_) + rnorm(1,0,sqrt(sigmasq))))
      }
    },
    GP_BQ = function(MI_iter=10000, min=0, max=1){
      integral = list()
      k_mean = as.matrix(kernel_mean_MI(MI_iter=MI_iter, min=0, max=1))
      full_k_int = full_kernel_integral_MI(MI_iter=MI_iter, min=0, max=1)
      integral$E = as.numeric(k_mean %*% C_N_inv %*% as.matrix(y))
      integral$V = full_k_int - as.numeric(k_mean %*% C_N_inv %*% t(k_mean))
      return(integral)
    }
    )
)




