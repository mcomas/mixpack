<<<<<<< HEAD
#' @title Orthonormal basis for the Simplex space
#' 
=======



>>>>>>> 2880484322866831e3a7885f8e3c0ae4d17692b3
#' Basis from the simplex space with $D$ components
#'
#' @param D numbre of components
#' @export
ilr_basis = function(D){
  lapply(1:(D-1), function(i){
    I = i+1
    l = exp(1 / sqrt(i*I))
    r = 1 / exp(sqrt(i/I))
    s_i = i * l + r + D - I
    c(rep(l, i),
      r,
      rep(1, D-I) ) / s_i
  })
}
<<<<<<< HEAD

#' @title Coordinates for an orthonormal basis
#' 
=======
>>>>>>> 2880484322866831e3a7885f8e3c0ae4d17692b3
#' Coordinates respect basis \code{\link{ilr_basis}}
#'
#' @param X compositional sample
#' @export
ilr_coordinates = function(X){
  D = ncol(X)
  r = data.frame(do.call('cbind', lapply(1:(D-1), function(i){
    if(i == 1){
      XN = matrix(X[,1], ncol=1)
    }else{
      XN = X[,1:i]
    }
    as.numeric(sqrt(i/(i+1)) * log( apply(XN, 1, prod)^(1/i) / X[,i+1] ))
  })))
  names(r) = paste('coord', 1:(D-1), sep='.')
  r
}

