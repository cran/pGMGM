emjglcov <-
function(Y, k, lambda1, lambda2, ncluster, tau, graph, threshold, MAX_iter,seed){
  n=nrow(Y)
  nlam1 = length(lambda1)
  nlam2 = length(lambda2)
  gtotal = length(ncluster)
  p = ncol(Y)
  pred_loglik = array(0, dim = c(nlam1, nlam2, gtotal))
  max_predll = -Inf 

  props <- rep(1/k,k)
  # How large should each group be?
  ns <- round(n * props)
  ns[k] <- n - sum(ns[-k])
  ids <- rep(1:k, ns)
  # Shuffle ids so that the groups are randomized
  set.seed(seed)
  whichgroup <- sample(ids)
  for(n1 in 1:nlam1){
    for(n2 in 1:nlam2){
      for(g in 1:gtotal){
        cvpred_loglik = rep(0, k)
        for(ki in 1:k){
          training = Y[(whichgroup!=ki),]
          tuning = Y[(whichgroup==ki),]
          cl_par = emcljgl(training, p, ncluster[g], lambda1[n1], lambda2[n2], threshold, MAX_iter, tau,graph,seed)
          # EM algorithm with known mu, covinv, but unknown tau_ij and pie_i
          cvpred_loglik[ki] = emcv(tuning, cl_par$mu, p, ncluster[g], cl_par$covinv, cl_par$pie)$ploglik                 
        }
        #print(cvpred_loglik) #print the pred log likelihood for each cross validation subset
        pred_loglik[n1, n2, g] = mean(cvpred_loglik)
        if(pred_loglik[n1, n2 , g] > max_predll){
          max_predll = pred_loglik[n1, n2 , g]
          par_optim = c(lambda1[n1], lambda2[n2],ncluster[g])
        }
      }
    }
  }  
  result = emcljgl(Y, p, par_optim[3], par_optim[1], par_optim[2], threshold, MAX_iter, tau, graph,seed) 
  # # add cv log likelihood
  # res = c(result, max_predll)
  return(result)
}
