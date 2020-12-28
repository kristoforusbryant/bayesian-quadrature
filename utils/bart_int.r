library('dbarts')
source("utils/intreegral.R")
source("utils/genz_functions.r")

rowVars <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

one_tree = function(trees, t){
  tree = trees[trees$tree == t,]
  return(intreegral(tree$var, tree$value, d))
}

n_trees = 200L # 200
n_samples = 1000L # 1000

rescale = function(y){(y - base::mean(y)) / max(abs(y - base::mean(y)) ) / 2}
rescale = function(y){(y - min(y)) / (max(y) - min(y)) - .5}

BART_INT <- setRefClass("BART_INT", 
    fields= list(
      d = "numeric",
      X = "matrix", 
      y_ = "numeric",
      integrand = "function", 
      control = "dbartsControl", 
      sampler = "ANY", 
      y_max = "numeric", 
      y_min = "numeric"
    ), 
    methods = list(
      initialize = function(d, integrand, control){
        .self$d = d
        .self$integrand = integrand
        .self$control = control
        .self$X = matrix(nrow=0, ncol=d)
        .self$y_max = .5
        .self$y_min = -.5
      },
      run = function(X_new, y_new){
        # X is a N x d matrix, y is a vector
    
        # Updating attributes
        X <<- rbind(X, X_new)
        y_ <<- c(y_, y_new)
        y_max <<- max(y_)
        y_min <<- min(y_)        
        # Running MCMC sampler
        sampler_ <- dbarts(formula = X, data = rescale(y_), tree.prior = cgm(2.0, 0.95), 
                          node.prior = normal(0.25/sqrt(n_trees)), 
                          n.samples = n_samples, resid.prior = chisq(3, 0.90)) 
        sampler_$setControl(control)

        sampler_$run()
        sampler <<- sampler_
      },
      mean = function(x){(rowMeans(sampler$predict(x)[,,1], dim=1) + .5) * (y_max - y_min) + y_min},
      var = function(x){rowVars(sampler$predict(x)[,,1], dim=1) * (y_max - y_min)**2},
      seq_design = function(n_of_pts=20*d, S=10000){
        for (i in 1:n_of_pts){
          c = matrix(runif(S * d), ncol = d)
          vars = .self$var(c)
          c_ = c[which.max(vars),,drop=FALSE]
          run(c_, integrand(c_))
        }
      }, 
      bart_int = function(){
        int_trees = sapply(1:control@n.samples, function(x_){  
          trees = BI_obj$sampler$getTrees(1:n_trees, 1, x_) 
          sum(sapply(1:n_trees, function(t){one_tree(trees, t)}))
        })
        res <- list()
        res$E <- (base::mean(int_trees) + .5) * (y_max - y_min) + y_min
        res$V <- stats::var(int_trees)  * (y_max - y_min)**2
        return(res)
      }, 
      bart_int_unscaled = function(){
        int_trees = sapply(1:control@n.samples, function(x_){  
          trees = BI_obj$sampler$getTrees(1:n_trees, 1, x_) 
          sum(sapply(1:n_trees, function(t){one_tree(trees, t)}))
        })
        res <- list()
        res$E <- base::mean(int_trees)
        res$V <- stats::var(int_trees)  
        return(res)
      }
    )
)

