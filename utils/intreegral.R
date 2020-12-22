library(rstack)

intreegral <- function(vars, vals, d){
  stopifnot(length(vars) == length(vals), max(vars) <= d)
  
  # initialisation 
  integral = 0
  left_empty = TRUE
  moves = stack$new() # 0 means left, 1 means right 
  moves$push(-1L) # as buffer
  history <- stack$new()
  history$push(-1L) # as buffer
  min_stacks = list()
  max_stacks = list()
  for (i in 1:d){
    min_stacks[[i]] <- stack$new()
    min_stacks[[i]]$push(-1L) #as buffer
    min_stacks[[i]]$push(0L)
    max_stacks[[i]] <- stack$new()
    max_stacks[[i]]$push(-1L) #as buffer
    max_stacks[[i]]$push(1L)
  }
  
  # main 
  for (i in 1:length(vars)){
    if (vars[i] > 0){
      max_stacks[[vars[i]]]$push(vals[i])
      history$push(vars[i])
      left_empty = TRUE
      moves$push(0)
    }
    else{
      if (left_empty){
        integral <- integral + vals[i] * 
          prod(sapply(max_stacks, function(x){x$peek()}) -
                 sapply(min_stacks, function(x){x$peek()}))
        now = history$peek()
        if (now < 0){break} # only one leaf 
        min_stacks[[now]]$push(max_stacks[[now]]$pop())
        left_empty = FALSE
        moves$pop()
        moves$push(1)
      }
      else{
        integral <- integral + vals[i] * 
          prod(sapply(max_stacks, function(x){x$peek()}) -
                 sapply(min_stacks, function(x){x$peek()}))
        if (i == length(vars)){break}
        while(moves$peek()){
          hist_pop = history$pop()
          min_stacks[[hist_pop]]$pop()
          moves$pop()
        }
        now = history$peek()
        if (now < 0){break} # only one leaf 
        min_stacks[[now]]$push(max_stacks[[now]]$pop())
        
        moves$pop() # deleting the excess left node
        moves$push(1)
      }
    }
  }
  return(integral)
}

intreegral_debug <- function(vars, vals, d){
  stopifnot(length(vars) == length(vals), max(vars) <= d)
  
  # initialisation 
  integral = 0
  left_empty = TRUE
  moves = stack$new() # 0 means left, 1 means right 
  moves$push(-1L) # as buffer
  history <- stack$new()
  history$push(-1L) # as buffer
  min_stacks = list()
  max_stacks = list()
  for (i in 1:d){
    min_stacks[[i]] <- stack$new()
    min_stacks[[i]]$push(-1L) #as buffer
    min_stacks[[i]]$push(0L)
    max_stacks[[i]] <- stack$new()
    max_stacks[[i]]$push(-1L) #as buffer
    max_stacks[[i]]$push(1L)
  }
  
  # main 
  for (i in 1:length(vars)){
    print(paste("###ITERATION: ", i))
    print("MIN_STACKS")
    print(sapply(min_stacks, function(x){x$peek()}))
    print("MAX_STACKS")
    print(sapply(max_stacks, function(x){x$peek()}))
    print(paste("HISTORY: ", history$peek()))
    print(paste("MOVES: ", moves$peek()))
    
    if (vars[i] > 0){
      print("---ADDING CONDITIONS---")
      max_stacks[[vars[i]]]$push(vals[i])
      history$push(vars[i])
      left_empty = TRUE
      moves$push(0)
    }
    else{
      if (left_empty){
        print("---CALCULATING INTEGRAL: LEFT EMPTY---")
        integral <- integral + vals[i] * 
          prod(sapply(max_stacks, function(x){x$peek()}) -
                 sapply(min_stacks, function(x){x$peek()}))
        now = history$peek()
        if (now < 0){break} # only one leaf 
        min_stacks[[now]]$push(max_stacks[[now]]$pop())
        left_empty = FALSE
        moves$pop()
        moves$push(1)
      }
      else{
        print("---CALCULATING INTEGRAL: LEFT NOT EMPTY---")
        integral <- integral + vals[i] * 
          prod(sapply(max_stacks, function(x){x$peek()}) -
                 sapply(min_stacks, function(x){x$peek()}))
        if (i == length(vars)){break}
        while(moves$peek()){
          hist_pop = history$pop()
          print(paste("POPPING HISTORY: ", hist_pop))
          min_stacks[[hist_pop]]$pop()
          moves$pop()
        }
        now = history$peek()
        if (now < 0){break}
        min_stacks[[now]]$push(max_stacks[[now]]$pop())
        
        moves$pop() # deleting the excess left node
        moves$push(1)
      }
      print(paste("INTEGRAL: ", integral))
    }
  }
  
  print("---END LOOP---")
  print(paste("HISTORY SIZE: ", history$size()))
  print("MIN_STACKS SIZE: ")
  print(sapply(min_stacks, function(x_){x_$size()}))
  print("MAX_STACKS SIZE: ")
  print(sapply(max_stacks, function(x_){x_$size()}))
  print(paste("MOVES SIZE: ", moves$size()))
  return(integral)
}


# TESTS
##d = 2
##vars <- c(2,1,-1,-1,-1)
##vals <- c(.5,.5, 1, 1, 1)
##print(intreegral(vars, vals, d)) # 1
##vars <- c(2,1,-1, 2,-1,-1,-1)
##vals <- c(.5,.5, 1, .5, 1, 1, 1) 
##print(intreegral(vars, vals, d)) # 1
##vars <- c(2,1,-1, 2,-1,-1,-1)
##vals <- c(.5,.5, 1, .3, 1, 1, 1)
##print(intreegral(vars, vals, d)) # 1
##vars <- c(2,1,-1, 2,-1,-1,-1)
##vals <- c(.5,.5, 1, .3, 2, 3, 4)
##print(intreegral(vars, vals, d)) # 2.85 

##d=3 
##vars <- c(1,2,3,-1,-1,-1,-1)
##vals <- c(.7,.2,.3, 1, 1, 1, 1)
##print(intreegral(vars, vals, d)) # 1 

##d=10
##vars <- c(1,2,3,-1,-1,-1,-1)
##vals <- c(.7,.2,.3, 1, 1, 1, 1)
##print(intreegral(vars, vals, d)) # 1 



