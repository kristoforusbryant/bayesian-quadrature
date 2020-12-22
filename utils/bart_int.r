library('dbarts')
source("utils/intreegral.R")
source("utils/genz_functions.r")

rowVars <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

one_tree = function(trees,t){
  tree = trees[trees$tree == t,]
  return(intreegral(tree$var, tree$value, d))
}

BART_INT <- setRefClass("BART_INT", 
    fields= list(
      d = "numeric",
      X = "matrix", 
      y = "numeric",
      integrand = "function", 
      control = "dbartsControl", 
      sampler = "ANY"
    ), 
    methods = list(
      initialize = function(d, integrand, control){
        .self$d = d
        .self$integrand = integrand
        .self$control = control
        .self$X = matrix(nrow=0, ncol=d)
      },
      run = function(X_new, y_new){
        # X is a N x d matrix, y is a vector
        ## !!! Make sure y is scaled appropriately (mean 0, between -.5 and .5) !!!
        X <<- rbind(X, X_new) 
        y <<- c(y, y_new) 
        sampler_ <- dbarts(formula = X, data = y, tree.prior = cgm(2.0, 0.95), 
                          node.prior = normal(0.25/sqrt(200)), 
                          n.samples = 1000L, resid.prior = chisq(3, 0.90)) # n_tree = 200, m = 1000
        sampler_$setControl(control)
        sampler_$run()
        sampler <<- sampler_
      }, 
      mean = function(x){rowMeans(sampler$predict(x)[,,1], dim=1)},
      var = function(x){rowVars(sampler$predict(x)[,,1], dim=1)},
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
          trees = BI_obj$sampler$getTrees(1:200, 1, x_)
          sum(sapply(1:200, function(t){one_tree(trees, t)}))
        })
        res <- list()
        res$E <- mean(int_trees)
        res$V <- var(int_trees) 
        return(res)
      }
    )
)

