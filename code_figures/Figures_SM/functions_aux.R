## Functions for the comparison between molecular and fitted branch lengths and predicted biomasses
library(ape)

x_k_r <- function(presence, A, nsp){
  Ak <- A[presence > 0, presence > 0]
  onesk <- presence[presence > 0]
  return(solve(Ak, onesk))
}

get_endpoints_r <- function(P, A, nsp){
  X <- P
  for (i in 1:ncol(P)){
    xtmp <- x_k_r(P[,i], A, nsp[i])
    X[P[,i] > 0,i] <- xtmp
  }
  return(X)
}

# Function to expand the tree structure from ape package
expand_tree <- function(original_tree){
  tree_bipart <- prop.part(original_tree)
  new_str <- tree_bipart
  # 1. Identify the cherries and work with the first one (the position of the next cherry will change after you append)
  cherries <- which(sapply(new_str, length) == 2)
  
  for(split_num in 1:length(cherries)){
    flag <- 0
    # Separate the cherry
    for(l in 1:2){
      new_str <- append(new_str, new_str[[cherries[split_num]]][l], cherries[split_num]+l-1)
      flag <- flag+1
    }
    
    # Checking if the first element of the list is a certain number
    for(linenum in cherries[split_num]:3){
      leaf_test <- end(new_str[[linenum-1]])[1]
      if(sum(sapply(new_str, function(x) all(x[[1]] == leaf_test, na.rm = TRUE))) ==0){
        # Where to append it
        pos <- tail(which(sapply(new_str, length)==1), n=1)
        new_str <- append(new_str, leaf_test, pos)
        flag <- flag+1
      }
    }
    
    cherries <- cherries+flag
  }
  return(new_str)  
}

# Getting the corresponding order between the molecular tree and our tree
get_col_redorder <- function(V, new_str){
  # Get the labels from V
  labels <- apply(t(V), 2, function(x) (1:nrow(t(V)))[x>0])
  
  #the other labels are stored in the list target
  #Find which labels are associated with the target
  
  reordering <- unlist(lapply(new_str, function(y) which(lapply(labels, function(x) all.equal(y, x)) == TRUE)))
  return(reordering)
}

# Code to plot predicted and observed biomasses

# Plotting individual abundances
plot_biomass <- function(results, predictions, label_optional = "molecular tree"){
  observed <- as.vector(results$Obs)
  predicted <- as.vector(predictions)
  predicted <- predicted[observed > 0]
  observed <- observed[observed > 0]
  main_label <- results$label
  
  correlation <- cor(log(observed), log(predicted))
  dt_ind <- tibble("predicted biomass" = predicted, "observed biomass" = observed)
  label_ind <- paste(main_label, label_optional, "correlation", round(correlation, 3))
  pl_ind <- ggplot(data = dt_ind) + aes(x = `observed biomass`, y = `predicted biomass`) + 
    geom_point() + geom_abline(slope = 1, intercept = 0, linetype = 2) + 
    scale_x_log10() + scale_y_log10() + theme_bw() + ggtitle(label_ind)
  show(pl_ind)
  return(pl_ind)
}
