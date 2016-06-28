obj_L0 <-
function(Y, mu, p, L, covmat_inverse, S_bar, graph, lambda1, lambda2, tau, pie) {
  log.liklyhood <- 0.0
  liklyhood = 0.0
  penalty <- 0.0
  n=nrow(Y)
  nn = round(pie*n, 0)
  #S_bar <- matrix(0,p,L*p) 
  #  for (l in 1:L){
  #    temp_Y = Y[memship==l] #get the lth cluster data, what if no data in lth cluster?
  #    stemp = t(t(temp_Y)-mu_temp)
  #    S_bar[,((l-1)*p+1):(l*p)] <- (t(stemp) %*% stemp) / cl$size[l]
  #  }
  #covmat_inverse = matrix(0, p, L*p)  
  #n = sum(nn)
  for (j in 1:n){
    liklyhood = 0.0
    for (l in 1:L){
      tmp_covmat <- covmat_inverse[,((l-1)*p+1):(l*p)]
      liklyhood = liklyhood + pie[l] * exp(fmvnorm(p, tmp_covmat, Y[j,], mu[l,])) #equation (2) in EJS Paper
      #log.liklyhood <- log.liklyhood + nn[l] * (- log(det(tmp_covmat)) 
      #                                          + sum(tmp_covmat*tmp_S))
    }
    log.liklyhood = log.liklyhood + log(liklyhood)
  }
  for (l in 1:L) {
    tmp_covmat <- covmat_inverse[,((l-1)*p+1):(l*p)]
    tmp_S <- S_bar[,((l-1)*p+1):(l*p)]
    tmp_covmat[abs(tmp_covmat) > tau] <- tau   
    penalty <- penalty + lambda1 * (sum(apply(abs(tmp_covmat),2,sum)) - sum(abs(diag(tmp_covmat))))
  }
  graph <- graph + 1 # switch to R indexing #
  Num.of.edges <- dim(graph)[2]
  for (e in 1:Num.of.edges) {
    l1 <- graph[1,e]
    l2 <- graph[2,e]
    tmp_covmat <- covmat_inverse[,(l1*p+1):((l1+1)*p)] - covmat_inverse[,((l2-1)*p+1):(l2*p)]
    tmp_covmat[abs(tmp_covmat) > tau] <- tau
    penalty <- penalty + lambda2 * (sum(apply(abs(tmp_covmat),1,sum)) - sum(abs(diag(tmp_covmat))))   
  }
  #return (log.liklyhood + penalty / 2)
  return (log.liklyhood - penalty)
}
