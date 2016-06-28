emcv <-
function(Y, mu, p, ncluster, covinv, pie){ 
  # EM algorithm for cross-validation with known mu, covinv, pie_i
  # initialize clustering using kmeans
  n = nrow(Y)
  cvlogliksum = 0
  for(j in 1:n){
    temp = 0
    for(i in 1:ncluster){
      temp = temp + pie[i] * exp(fmvnorm(p, covinv[,((i-1)*p+1):(i*p)], Y[j,],mu[i,]))
    }
    cvlogliksum = cvlogliksum + log(temp)  
  }
  output = list(ploglik = cvlogliksum) 
}
