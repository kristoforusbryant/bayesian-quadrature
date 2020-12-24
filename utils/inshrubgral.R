
# V2 of intreegral 

# Given two lists vars and vals, we compute the probability associated with each leaf, and then return
# it as a sum. The inputs vars and vals contain information that allows the information
# to be visualised as a tree, ordered in DFS manner.
# vars: contains integer values comprising of -1, and 1: number of variables used for each node split.
# vals: contains numeric values. When its associated var value is greater than 0, this value represents
# the constraint used for each node split. 
# 
# Expected complexity of this algorithm is O(n), where n is the length of vars.
# Sketch of algorithm:
#   
# At each row of vars and vals, if we are not at the leaf, we are at a node (indicated by value of var).
# This means the val associated with this row is a constraint, and we can compute this node's subsequent left and right
# node parameters. For each, initialise an additional parameter indicating whether its the left (1) or right (0) parameter. This will
# be useful for backtracking later. Store these branching nodes' parameters in the order right, then left, to a list current_param. 
# 
# When we encounter the first leaf, we need to start backtracking. We take note of this by using the indicator leaf_visited.
# 
# If we are at a leaf (var=-1), pop the latest entry of current_param and append it to leaf_param.
# If we are backtracking and realise the next value in our stack is a node visited before (left node), pop it from the stack. Do this
# until the subsequent value in our stack has not been visited before. Finally, change this left/right param to 1, as we will
# visit it next. 


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
  current_param = list(c(rep(c(0,1), (d+1)), 1))
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
      pi_ = current_val * prod(max_vals - min_vals)
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


# inshrubgral <- function(vars, vals, d){
#   
#   #initialise 
#   current_param = list(rep(1, (d+1)))
#   leaf_param = list()
#   leaf_visited = 0
#   
#   vars_len = length(vars)
#   vals_len = length(vals)
#   
#   #sanity check
#   if (vars_len != vals_len){
#     stop("Error: vars and vals must be of same length")
#   }
#   
#   for (i in 1:vars_len){
#     current_var = unlist(vars[i])
#     current_val = unlist(vals[i])
#     
#     if (current_var != -1){
#       param_to_change = current_var
#       prev_param = tail(current_param, n=1)
#       
#       #for left branch
#       new_left = unlist(prev_param)
#       new_left[[param_to_change]] = current_val
#       #mark visit
#       new_left[length(new_left)] = 1
#       
#       #for right branch
#       new_right = unlist(prev_param)
#       new_right[[param_to_change]] = new_right[[param_to_change]] - current_val
#       #mark that visit has not taken place
#       new_right[length(new_right)] = 0
#       
#       # append node values to list (right first, then left)
#       push(current_param, new_right)
#       push(current_param, new_left)
#       next
#     }
#     if (current_var == -1){
#       leaf_node = unlist(tail(current_param, n=1))
#       push(leaf_param, prod(current_val, head(leaf_node, n=length(leaf_node)-1)))
#       pop(current_param)
#       leaf_visited = 1
#     }
#     if (leaf_visited == 1){
#       # base case
#       if (i == length(vars)){
#         break
#       }
#       else {
#         while ((tail(unlist(tail(current_param, n=1)), n=1) == 1) && (i != length(vars))){
#           pop(current_param)
#         }
#         # reset visit index to 1 for "right branch" because it is visited next
#         last_node = unlist(tail(current_param, n=1))
#         pop(current_param)
#         last_node[length(last_node)] = 1
#         push(current_param, last_node)
#       }
#     }
#   }
#   return(sum(unlist(leaf_param)))
# }
# 
# 
# 
