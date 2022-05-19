basic_norm <- function(weights){
  prod(colMeans(weights))
}

harmonic_norm <- function(g,weights,varphi){
  weight_prod <- apply(weights,1,prod)
  (sum(abs(g)^(1-varphi)*weight_prod))/
    (sum(abs(g)^(-varphi)*weight_prod))
}