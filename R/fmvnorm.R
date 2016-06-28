fmvnorm <-
function(p, covmat_inv, Y, mu){ # Y n*p, Y[j,]; mu L*p, mu[l,]
  #mvnormden = -0.5 * p * log(2*pi) + 0.5 * log(det(covmat_inv )) -0.5 * t(Y-mu) %*% covmat_inv %*% (Y-mu)
  #mvnormden = -p * log(2*pi) + log(det(covmat_inv )) - drop(t(Y-mu) %*% covmat_inv %*% (Y-mu))
  testread = try(det(covmat_inv), silent=T)
  if(class(testread)!="try-error"){
    ldcovmat = log(det(covmat_inv ))
  }else{
    ldcovmat = 1e-8
  }
  mvnormden = -p * log(2*pi) + ldcovmat - drop(t(Y-mu) %*% covmat_inv %*% (Y-mu))
  return(0.5*mvnormden) # return the log likelihood
}
