require(plyr)
require(gridExtra)
require(ggplot2)

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
  names(new_partition) = laply(new_partition, part_name)
  new_partition
}

#' Create a partition of classes from weights or probabilities
#'
#' This function applies the methodology described in [citar article]
#' to build a hierarchy of classes using the weights or probabilities 
#' that an element belongs to each class
#' @param tau dataframe of probabilities/weights (\code{tau} must be strictly positive)
#' 
#' @param varphi function with two parameters (\code{v_tau}, \code{a}). Parameter 
#' \code{v_tau} is a vector of probabilities, parameter \code{a} is the a selected class.
#' \code{varphi}(\code{v_tau}, \code{a}) gives the representativeness of element with
#' probabities \code{v_tau} to class \code{a}
#' @export
#' @param theta function with three parameters (\code{v_tau}, \code{a}, \code{b}).
#' Parameter \code{v_tau} is a vector of probabilities, parameters \code{a} and \code{b}
#' are classes to be combined.
#' 
get_hierarchical_partition = function(tau, 
                                      varphi, 
                                      theta){
  ctau = tau
  K = ncol(ctau)
  partitions = list()
  partitions[[K]] = as.list(1:K)
  names(partitions[[K]]) = laply(partitions[[K]], part_name)
  for(k in K:2){
    COMB = t(expand.grid(1:k, 1:k))
    COMB = COMB[, COMB[1,] != COMB[2,]]
    rownames(COMB) = c('a', 'b')
    colnames(COMB) = col.names=apply(COMB, 2, paste, collapse='-')
    to_merge = which.max( v <- aaply(COMB, 2, function(ind){
      a = ind[1]; b = ind[2]
      sum( apply(ctau, 1, function(v_tau) varphi(v_tau, a) * theta(v_tau, a, b) ) ) / sum( apply(ctau, 1, function(v_tau) varphi(v_tau, a) ) )
    }) )
    part = COMB[,to_merge]
    partitions[[k-1]] = b_absorbes_a(partitions[[k]], part['a'], part['b'] )
    ctau[,part['b']] = ctau[,part['a']] + ctau[,part['b']]
    ctau  = ctau[,-part['a']]
  }
  class(partitions) = 'hpartition'
  partitions
}

get_hierarchical_partition_mult_1 = function(tau, 
                                      varphi = function( v_tau, a) if(which.max(v_tau) == a) 1 else 0, 
                                      theta = function(v_tau, a, b) log(v_tau[a] / v_tau[b])^2){
  ctau = tau
  K = ncol(ctau)
  partitions = list()
  partitions[[K]] = as.list(1:K)
  names(partitions[[K]]) = laply(partitions[[K]], part_name)
  for(k in K:2){
    COMB = t(expand.grid(1:k, 1:k))
    COMB = COMB[, COMB[1,] != COMB[2,]]
    rownames(COMB) = c('a', 'b')
    colnames(COMB) = col.names=apply(COMB, 2, paste, collapse='-')
    to_merge = which.min( v <- aaply(COMB, 2, function(ind){
      a = ind[1]; b = ind[2]
      sum( apply(ctau, 1, function(v_tau) varphi(v_tau, a) * theta(v_tau, a, b) ) ) / sum( apply(ctau, 1, function(v_tau) varphi(v_tau, a) ) )
    }) )
    part = COMB[,to_merge]
    partitions[[k-1]] = b_absorbes_a(partitions[[k]], part['a'], part['b'] )
    ctau[,part['b']] = sqrt(ctau[,part['a']] * ctau[,part['b']])
    ctau  = ctau[,-part['a']]
  }
  class(partitions) = 'hpartition'
  partitions
}

get_hierarchical_partition_mult_2 = function(tau, 
                                             varphi,# = function( v_tau, a) if(which.max(v_tau) == a) 1 else 0, 
                                             theta){# = function(v_tau, a, b) log(v_tau[a] / v_tau[b])^2){
  ctau = tau
  K = ncol(ctau)
  partitions = list()
  partitions[[K]] = as.list(1:K)
  names(partitions[[K]]) = laply(partitions[[K]], part_name)
  for(k in K:2){
    COMB = t(expand.grid(1:k, 1:k))
    COMB = COMB[, COMB[1,] != COMB[2,]]
    rownames(COMB) = c('a', 'b')
    colnames(COMB) = col.names=apply(COMB, 2, paste, collapse='-')
    to_merge = which.min( v <- aaply(COMB, 2, function(ind){
      a = ind[1]; b = ind[2]
      sum( apply(ctau, 1, function(v_tau) varphi(v_tau, a) * theta(v_tau, a, b) ) ) / sum( apply(ctau, 1, function(v_tau) varphi(v_tau, a) ) )
    }) )
    part = COMB[,to_merge]
    partitions[[k-1]] = b_absorbes_a(partitions[[k]], part['a'], part['b'] )
    ctau = prop_partition_mult(tau, partitions[[k-1]])
  }
  class(partitions) = 'hpartition'
  partitions
}

hp_entropy = function(tau){
  varphi = function(v_tau, a) 1
  theta = function(v_tau, a, b) -xlog(v_tau[a] + v_tau[b]) + xlog(v_tau[a]) + xlog(v_tau[b])
  get_hierarchical_partition(tau, varphi, theta)
}
hp_entropy.prop = function(tau){
  varphi = function(v_tau, a) v_tau[a]
  theta = function(v_tau, a, b) -xlog(v_tau[a] + v_tau[b]) + xlog(v_tau[a]) + xlog(v_tau[b])
  get_hierarchical_partition(tau, varphi, theta)
}
hp_entropy.dichotomic = function(tau){
  varphi = function(v_tau, a) if(which.max(v_tau) == a) 1 else 0
  theta = function(v_tau, a, b) -xlog(v_tau[a] + v_tau[b]) + xlog(v_tau[a]) + xlog(v_tau[b])
  get_hierarchical_partition(tau, varphi, theta)
}
hp_atchison.prop = function(tau){
  varphi = function(v_tau, a) v_tau[a]
  theta = function(v_tau, a, b) log(v_tau[a] / v_tau[b])^2
  get_hierarchical_partition(tau, varphi, theta)
}
hp_atchison.dichotomic = function(tau){
  varphi = function(v_tau, a) if(which.max(v_tau) == a) 1 else 0
  theta = function(v_tau, a, b) log(v_tau[a] / v_tau[b])^2
  get_hierarchical_partition(tau, varphi, theta)
}
hp_demp = function(tau){
  varphi = function(v_tau, a) v_tau[a]
  theta = function(v_tau, a, b) -if(which.max(v_tau) == b) 1 else 0
  get_hierarchical_partition(tau, varphi, theta)
}
hp_demp2 = function(tau){
  varphi = function(v_tau, a) v_tau[a]
  theta = function(v_tau, a, b) if(which.max(v_tau) != b) 1 else 0
  get_hierarchical_partition(tau, varphi, theta)
}
hp_demp.prop = function(tau){
  varphi = function(v_tau, a) v_tau[a]
  theta = function(v_tau, a, b) -v_tau[b]
  get_hierarchical_partition(tau, varphi, theta)
}

prop_partition = function(tau, partition)
  do.call('cbind', llply(partition, function(part){
    if(is.vector(tau[,part])) return(tau[,part])
    apply(tau[,part], 1, sum)
  }))

prop_partition_mult = function(tau, partition)
  do.call('cbind', llply(partition, function(part){
    if(is.vector(tau[,part])) return(tau[,part])
    apply(tau[,part], 1, prod)^(1/length(part))
  }))

cluster_partition = function(tau, partition)
  names(partition)[apply( do.call('cbind', llply(partition, function(part){
     if(is.vector(tau[,part])) return(tau[,part])
     apply(tau[,part], 1, sum)
     })), 1, which.max)]

plot.hpartition = function(hp, tau, data, nrow = 2){
  df = ldply(hp, function(partition){
    df = data.frame(data)
    df$cluster = cluster_partition(tau, partition)
    df$step = factor(length(partition), levels =  length(hp):1)
    df
  })
  
  plots = llply(split(df, df$step), function(d) 
    ggplot(data=d, aes(x=X1, y=X2, col=cluster)) + 
      geom_point(size=2) + xlab("") + ylab("") + 
      theme(legend.position="none"))#legend.title=element_blank()))
  
  do.call("grid.arrange", c(plots, nrow=nrow))
}

xlog <- function(x) {
  xlog1d <- function (xi) if (xi == 0) 0 else (xi*log(xi))
  
  if (is.null(dim(x))) return(sapply(x,xlog1d))
  else return(matrix(sapply(x,xlog1d),dim(x)))
}
