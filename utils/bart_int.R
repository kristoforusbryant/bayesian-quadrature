library("dbarts")
setwd('Desktop/bayesian-quadrature/utils/')
source("intreegral.R")



bart_sequential_fitting <- function(integral_function, a, x_train, y_train, n_total, n_tree, S = 5000, d, m){
  
  # compute number of additional meshpoints required
  nseq <- n_total- length(y_train) 
  
  add_mesh <- 0
  
  
  # prefix the settings of the dbart model
  control_ <- dbartsControl(verbose = FALSE, keepTrainingFits = TRUE, useQuantiles = FALSE,
                            keepTrees = TRUE, n.samples = m,
                            n.burn = 1000, n.trees = n_tree, #n.chains = 1,
                            n.threads = guessNumCores(), n.thin = 5, printEvery = 100,
                            rngKind = "default", rngNormalKind = "default",
                            updateState = TRUE)
  
  while (add_mesh < nseq){
    
    # dim of candidate set is S x d
    
    candidate_set <- matrix(runif(S * d), ncol = d) 
    
    
    sampler <- dbarts(formula = x_train, data = y_train, tree.prior = cgm(2.0, 0.95), node.prior = normal(0.25/sqrt(n_tree)), 
                      n.samples = m, resid.prior = chisq(3, 0.90))
    
    sampler$setControl(control_)
    
    sampler$run()
    
    # get meshpoint with greatest variance
    
    pred <- sampler$predict(candidate_set)
    pred_chain <- pred[,,1]
    
    
    expected_mean <- apply(pred_chain, 1, mean)
    candidate_var <- apply(pred_chain, 1, var)
    
    
    c_star <- candidate_set[which.max(candidate_var), ]

    
    #evaluate y_star
    y_star <- integral_function(rep(a, d), c_star, rep(0.5, d)) + rnorm(1, 0, 0.1) #epsilon_star
    
    x_train <- rbind(x_train, c_star)
    y_train <- c(y_train, y_star)
    
    add_mesh <- add_mesh + 1
    

    
  }
  # estimate the integrals
  
  # find stats of dbarts model
  trees <- sampler$getTrees(chainNums=1, treeNums=1:200, sampleNums=1:m)
  df <- split(trees, with(trees, interaction(tree, sample)), drop=TRUE)
  result <- lapply(df, function(x){intreegral_fix(x$var, x$val, d)})
  names(result) <- sub("^.*\\.", "", names(result))

  int_sum <- (1/m) * sum(tapply(unlist(result), names(unlist(result)), sum))
  
  return(list(x_train, y_train, int_sum))
}


run_iterations <- function(integral_function, genz_func_name, a, n_total, n_tree, S = 5000, d, m,  n_iter=10){
  
  output_dir <- file.path("genz_function_results")
  
  if (!dir.exists(output_dir)) {dir.create(output_dir, recursive=TRUE)}
  
  pi_f_global = list()
  x_train_global = list()
  y_train_global = list()
  time_taken_global = list()
  
  for (i in (1:n_iter)){
    
    set.seed(i*47)
    x_train <- matrix(runif(N * d), ncol = d)
    y_train <- apply(x_train, 1, function(x_){integral_function(a, x_, u)}) 
    
    set.seed(i)
    
    ptm <- proc.time()

    iter_output <- bart_sequential_fitting(integral_function, a, x_train, y_train, n_total, n_tree, S, d, m)
   
     time_taken <- getElement(proc.time() - ptm, "elapsed")
     print(paste("Time taken on iter: ", i))
    
     print(time_taken)
    
    
    x_train_global[[i]] <- iter_output[1]
    y_train_global[[i]] <- iter_output[2]
    pi_f_global[[i]] <- iter_output[3]
    time_taken_global[[i]] <- time_taken
    
  }
  
  df_path <- paste(output_dir, "/", genz_func_name, ".csv" , sep="")
  
  write.csv(data.frame(pi_f_global, time_taken_global), df_path, row.names = FALSE)
  
  return(data.frame(pi_f_global, time_taken_global))
  }

