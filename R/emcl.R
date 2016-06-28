emcl <-
function(Y, p, ncluster, lambda1, lambda2, threshold, MAX_iter, tau, graph){
  # initialize clustering using kmeans
  L=ncluster 
  n=nrow(Y)
  # if(ncluster == 1){
  #    memship_ini = rep(1, ncluster)
  #    center_ini = apply(Y, 2, mean)
  #    pie = 1
  #  }else{
  kmeans_out = kmeans(x=Y, centers = ncluster, iter.max = 10, nstart = ncluster) 
  memship_ini = kmeans_out$cluster
  #print(memship_ini)
  center_ini = kmeans_out$centers
  pie = kmeans_out$size/n
  # }
  
  #for(l in 1:ncluster){ # initialize pie
  #  pie = sum(center_ini==l)/n
  #}
  
  mu = center_ini # initialize mu
  S_bar = matrix(0, p, L*p)
  for (l in 1:L){ # calculate S_bar
    temp_Y = Y[memship_ini==l,] #get the lth cluster data
    if(is.vector(temp_Y)) mu_temp = mean(temp_Y)
    else mu_temp = apply(temp_Y, 2, mean)
    stemp = t(t(temp_Y)-mu_temp)
    S_bar[,((l-1)*p+1):(l*p)] <- (t(stemp) %*% stemp) / kmeans_out$size[l]
  }
  nn=kmeans_out$size
  sol_path <- MGGM.path(S_bar, nn, lambda1, lambda2, graph, tau) # get the inv cov matrix
  covinv_est = sol_path$sol_nonconvex[,,1,1] # p row, p*l colu
  npie = pie
  nmu = mu
  ncovinv_est = covinv_est
  z = 1
  diff = 100
  nplog = obj_L0(Y, mu, p, L, covinv_est, S_bar, graph, lambda1, lambda2, tau, pie) 
  rho = matrix(0, L, n)
  while(diff>threshold && z<MAX_iter){
    for(j in 1:n){ # update tau_ij
      sm = rep(0, n)
      for(l in 1:L){ # compute the sum in the denominator
        sm[j] = sm[j] + pie[l] * exp(fmvnorm(p, covinv_est[,((l-1)*p+1):(l*p)], Y[j,], mu[l,]))
      }
      for(l in 1:L){ # compute tau_ij
        rho[l,j] = pie[l] * exp(fmvnorm(p, covinv_est[,((l-1)*p+1):(l*p)], Y[j,], mu[l,])) / sm[j]
      }   
    }
    # update pie_i
    pie = apply(rho, 1, mean) #(sum_j tau_ij = pie_i *n)
    nn = round(pie*n, 0)
    # update mu_i
    #for(l in 1:L){
    #  mu[l, ] = tau[l]%*%Y / (pie[l]*n)
    #}
    mu = rho %*% Y/ (pie * n) # each mu_i is divided by (sum_j tau_ij = pie_i *n)
    #calculate S_bar
    S_bar <- matrix(0,p,L*p)
    for (l in 1:L){
      S = matrix(0,p,p) 
      for(j in 1:n){
        S = S + rho[l,j] * outer((Y[j,]-mu[l,]), (Y[j,]-mu[l,]), "*")
      }
      S_bar[,((l-1)*p+1):(l*p)] <- S/(pie[l]*n)
    }
    # update inverse covariance matrix
    sol_path <- MGGM.path(S_bar, nn, lambda1, lambda2, graph, tau)
    #sol_nonconvex[,,j,i] is refering to the solution corresponding to grid point
    #(Lambda1.vec[i], Lambda2.vec[j])
    covinv_est = sol_path$sol_nonconvex[,,1,1] # p row, p*l column
    #temp_memship = apply(rho, 2, which.max)
    plog = obj_L0(Y, mu, p, L, covinv_est, S_bar, graph, lambda1, lambda2, tau, pie)  
    #diff = abs(plog - nplog)
    diff = plog-nplog
    nplog = plog
    z = z+1 
  }
  class_memship = apply(rho, 2, which.max)
  # output = list(pie = pie, mu = mu, covinv = covinv_est, ploglik = nplog, membership = class_memship, par=c(lambda1, lambda2, ncluster))
  output = list(pie = pie, mu = mu, covinv = covinv_est, membership = class_memship, par=c(lambda1, lambda2, ncluster))
}
