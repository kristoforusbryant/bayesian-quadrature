library("dbarts")

#'@param integral_function one of the genz functions to compute the integrand
#'@param a parameter 'a' characterising the integrand function
#'@param x_train Nxd matrix specifying the x data, where N is the number of datapoints and d is the dimension of the data
#'@param y_train Nx1 vector specifying the y data
#'@param n_total integer specifying the total number of meshpoints to approximate the integral.
#'@param n_tree integer specifying the number of trees to use in the bart model
#'@param S integer specifying the number of candidate points to consider for each sequential fit
#'@param d integer specifying the dimension of the x_train data
#'@param seed integer specifying the seed to ensure reproducibility of results


bart_sequential_fitting <- function(integral_function, a, x_train, y_train, n_total, n_tree, S = 10000, d, seed){
  
  # compute number of additional meshpoints required
  nseq <- n_total- length(y_train) 
  
  add_mesh <- 0
  set.seed(seed)
  
  while (add_mesh < nseq){
    
    # dim of candidate set is S x d
    candidate_set <- t(replicate(S, runif(d, 0, 1)))
    
    #fit bart model
    fit_result <- bart(x.train = x_train, y.train = y_train, x.test = candidate_set, sigest = NA, sigdf = 3, sigquant = 0.90,
                       k = 2.0, power = 2.0, base = 0.95, ntree= n_tree, nskip = 1000, keepevery = 5, ndpost = 1000)
    
    # find meshpoint that gives greatest variance
    candidate_var <- apply(fit_result$yhat.test, 2, var)
    c_star <- candidate_set[which.max(candidate_var), ]
    
    #evaluate y_star
    epsilon_star <- rnorm(1, mean=0, sd = mean(fit_result$sigma))
    y_star <- integral_function(rep(a, d), c_star, rep(0.5, d)) + epsilon_star
    
    
    x_train <- rbind(x_train, c_star)
    y_train <- c(y_train, y_star)
    
    add_mesh <- add_mesh + 1
    
  }
  return(list(x_train, y_train))
}


