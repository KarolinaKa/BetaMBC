# /////////////////////////////////////
# MOM univariate betas(m,v) ///////////
# /////////////////////////////////////


beta_rep_MOM <- function(data, lower = 0.01, upper = 100) {
  x_bar <- mean(data)
  N     <- length(data)
  s2    <- var(data) * (N - 1)/N
  
  f <- function(v)
    v*x_bar*(1 - x_bar)/(3*v + 1) - s2
  
  u <- uniroot(f, c(lower, upper))
  c(x_bar*(2*u$root + 1) - u$root, u$root)
}

comp <- sample(1:G,                       # components
               prob = p, 
               size = N, 
               replace = TRUE)          
# simulated data
sim <- rbeta(N, m[comp] / v[comp] + 1, (1 - m[comp]) / v[comp] + 1) 


g <- sapply(1:G, function(j) subset(sim, comp == j))
sapply(g, beta_rep_MOM)

