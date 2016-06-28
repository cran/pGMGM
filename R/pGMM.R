pGMM <- function(Y, k=5, method="New-SP", lambda1, lambda2, ncluster, tau=0.01, threshold = 1e-5, MAX_iter=100, seed=1){
  # A complete graph
  graph_complete = matrix(0,2,ncluster*(ncluster-1)/2)
  for (ni in 1:(ncluster-1)){
    graph_complete[,(ncluster*(ni-1)-(ni-1)*ni/2+1):(ncluster*ni-ni*(ni+1)/2)] = rbind(rep(ni,ncluster-ni),(ni+1):ncluster)
  }
  graph <- graph_complete - 1
  
  if(method == "New-SP"){
     out = emnccov(Y, k, lambda1, lambda2, ncluster, tau, graph, threshold, MAX_iter,seed =seed)
  }else if(method == "New-JGL"){
     out = emjglcov(Y, k, lambda1, lambda2, ncluster, tau, graph, threshold, MAX_iter,seed =seed)
  }else {
   cat("Method must be either 'New-SP' or 'New-JGL'", "\n")
  }
  covinvmat = out$covinv
  p = nrow(covinvmat)
  covinv = list()
  for(ni in 1:ncluster){
    covinv[[ni]] = covinvmat[, ((ni-1)*p+1):(ni*p)]
  }
  out$covinv = covinv
  return(out)
}


