require(ggplot2)

ggbiplot <- function(X, alpha=0, cols=NULL, comp=1:2, size=2, transparency = 0.15){
  
  s = svd(scale(X, scale=F))
  obs = data.frame(s$u %*% diag(s$d^alpha))
  names(obs) = paste0('Comp.', 1:ncol(obs))
  
  vars = data.frame(t(diag(s$d^(1-alpha) ) %*% t(s$v)))
  names(vars) = names(obs)
  vars$varnames = colnames(X)
  
  if(!is.null(cols))
    obs$g = as.factor(cols)
  
  # PC being a prcomp object
  x = names(obs)[comp[1]]
  y = names(obs)[comp[2]]
  if(is.null(cols)){
    plot <- ggplot(obs, aes_string(x=x, y=y))
  }else{
    plot <- ggplot(obs, aes_string(x=x, y=y, col='g'))
  }
  plot = plot + geom_point(size=size, alpha=transparency)
  
  #plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
  
  mxy_obs = max( abs(max(obs[,c(x,y)])), abs(min(obs[,c(x,y)])))
  mxy_vars = max( abs(max(vars[,c(x,y)])), abs(min(vars[,c(x,y)])))
  
  vars <- transform(vars,
                    v1 = .7 *  mxy_obs / mxy_vars * (get(x)),
                    v2 = .7 *  mxy_obs / mxy_vars * (get(y))
  )
  if(is.null(cols)){
    plot <- plot + geom_segment(data=vars, aes(x=0, y=0, xend=v1, yend=v2), 
                                #arrow=arrow(length=unit(0.2,"cm")), 
                                alpha=0.75, size=3)  
  }else{
    vars$g = factor(levels(obs$g))
    plot <- plot + geom_segment(data=vars, aes(x=0, y=0, xend=v1, yend=v2, color=g), 
                                #arrow=arrow(length=unit(0.2,"cm")), 
                                alpha=0.75, size=3)  
  }
  plot
  
}