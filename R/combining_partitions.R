library(plyr)

#' @title Build a hierchical partition from posterior probabilities
#' 
#' Create a partition of classes from weights or probabilities
#'
#' This function applies the methodology described in [citar article]
#' to build a hierarchy of classes using the weights or probabilities 
#' that an element belongs to each class
#' @param tau dataframe of probabilities/weights (\code{tau} must be strictly positive)
#' 
#' @param omega function with two parameters (\code{v_tau}, \code{a}). Parameter 
#' \code{v_tau} is a vector of probabilities, parameter \code{a} is the a selected class.
#' \code{omega}(\code{v_tau}, \code{a}) gives the representativeness of element with
#' probabities \code{v_tau} to class \code{a}
#' 
#' @param lambda function with three parameters (\code{v_tau}, \code{a}, \code{b}).
#' Parameter \code{v_tau} is a vector of probabilities, parameters \code{a} and \code{b}
#' are classes to be combined.
#' @export
get_hierarchical_partition = function(tau, 
                                      omega, 
                                      lambda){
  ctau = tau
  K = ncol(ctau)
  partitions = list()
  partitions[[K]] = as.list(1:K)
  names(partitions[[K]]) = plyr::laply(partitions[[K]], part_name)
  for(k in K:2){
    COMB = t(expand.grid(1:k, 1:k))
    COMB = COMB[, COMB[1,] != COMB[2,]]
    rownames(COMB) = c('a', 'b')
    colnames(COMB) = col.names = apply(COMB, 2, paste, collapse='-')
    to_merge = which.max( v <- plyr::aaply(COMB, 2, function(ind){
      a = ind[1]; b = ind[2]
      sum( apply(ctau, 1, function(v_tau) omega(v_tau, a) * lambda(v_tau, a, b) ) ) / sum( apply(ctau, 1, function(v_tau) omega(v_tau, a) ) )
    }) )
    part = COMB[,to_merge]
    partitions[[k-1]] = b_absorbes_a(partitions[[k]], part['a'], part['b'] )
    ctau[,part['b']] = ctau[,part['a']] + ctau[,part['b']]
    ctau  = ctau[,-part['a']]
  }
  class(partitions) = 'hpartition'
  partitions
}

part_name = function(part) sprintf("(%s)", paste(sort(part), collapse=','))

## PART B absorbes A. In this function part A is incorporated to B and after that partition A is eliminated
b_absorbes_a = function(partition, partA, partB){
  if(partA == partB){
    stop( "Same part A and B")
  }
  if(! (partA %in% 1:length(partition) & partB %in% 1:length(partition)) ){
    stop( "Some part out of range")
  }
  new_partition = partition
  new_partition[[partB]] = c(new_partition[[partA]], new_partition[[partB]])
  new_partition[[partA]] = NULL
  names(new_partition) = plyr::laply(new_partition, part_name)
  new_partition
}