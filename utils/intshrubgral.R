pop <- function(x){
  assign(as.character(substitute(x)), x[-length(x)], parent.frame())
}

push <- function(x, val){
  x_new = x
  x_new[[length(x_new)+1]] = val
  assign(as.character(substitute(x)), x_new, parent.frame())
}

inshrubgral <- function(vars, vals, d){
  
  #initialise 
  current_param = list(c(rep(c(0,1), (d)), 1)) #list(c(rep(c(0,1), (d+1)), 1))
  leaf_param = list()
  leaf_visited = 0
  
  vars_len = length(vars)
  vals_len = length(vals)
  
  #sanity check
  if (vars_len != vals_len){
    stop("Error: vars and vals must be of same length")
  }
  
  for (i in 1:vars_len){
    current_var = unlist(vars[i])
    current_val = unlist(vals[i])
    
    if (current_var != -1){
      param_to_change = current_var
      prev_param = tail(current_param, n=1)
      
      #for left branch
      new_left = unlist(prev_param)
      new_left[[2*param_to_change]] = current_val
      #mark visit
      new_left[length(new_left)] = 1
      
      #for right branch
      new_right = unlist(prev_param)
      new_right[[2*param_to_change-1]] = current_val #new_right[[param_to_change]] - current_val
      #mark that visit has not taken place
      new_right[length(new_right)] = 0
      
      # append node values to list (right first, then left)
      push(current_param, new_right)
      push(current_param, new_left)
      next
    }
    if (current_var == -1){
      leaf_node = unlist(tail(current_param, n=1))
      min_vals = leaf_node[seq(1, length(leaf_node)-1, 2)]
      max_vals = leaf_node[seq(2, length(leaf_node)-1, 2)]
      pi_ =  current_val * prod(max_vals - min_vals)
      push(leaf_param, pi_)
      pop(current_param)
      leaf_visited = 1
    }
    if (leaf_visited == 1){
      # base case
      if (i == length(vars)){
        break
      }
      else {
        while ((tail(unlist(tail(current_param, n=1)), n=1) == 1) && (i != length(vars))){
          pop(current_param)
        }
        # reset visit index to 1 for "right branch" because it is visited next
        last_node = unlist(tail(current_param, n=1))
        pop(current_param)
        last_node[length(last_node)] = 1
        push(current_param, last_node)
      }
    }
  }
  return(sum(unlist(leaf_param)))
}
