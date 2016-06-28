emcljgl = function(Y, p, ncluster, lambda1, lambda2, threshold, MAX_iter, tau, graph,seed){
  # initialize clustering using kmeans
  L=ncluster 
  n=nrow(Y)
#  kmeans_out = kmeans(x=Y, centers = ncluster, iter.max = 10, nstart = ncluster) 
#  memship_ini = kmeans_out$cluster
#  center_ini = kmeans_out$centers
#  pie = kmeans_out$size/n

    set.seed(seed)
    memship_ini = sample(seq(1,ncluster),n, replace =T)
    #print(memship_ini)
    clust_size = as.numeric(table(memship_ini))
    center_ini = matrix(0,nrow =ncluster, ncol=p)
    for(l in 1:L){
      center_ini[l,] = apply(Y[memship_ini==l,], 2, mean)
    }
    pie = clust_size/n

  dat1 = list()
  temp_num_cluster = length(unique(memship_ini))
  for(l in 1:temp_num_cluster){
     dat1[[l]] = Y[memship_ini == l,]
  }
  mu = center_ini # initialize mu
  sol_path <- JGL(dat1 , penalty = "fused", lambda1, lambda2, 0.1, return.whole.theta=T)
  covinv_est = matrix(0, nrow = p, ncol = temp_num_cluster*p)
  S_bar = matrix(0, p, L*p)
  for (l in 1:temp_num_cluster){ 
    covinv_est[,((l-1)*p+1):(l*p)] = sol_path$theta[[l]]
    # calculate S_bar
    temp_Y = dat1[[l]] #get the lth cluster data
    mu_temp = apply(temp_Y, 2, mean)
    stemp = t(t(temp_Y)-mu_temp)
    S_bar[,((l-1)*p+1):(l*p)] <- (t(stemp) %*% stemp) / clust_size[l]
  }
  #covinv_est = sol_path$sol_nonconvex[,,1,1] # p row, p*l colu
  npie = pie
  nmu = mu
  ncovinv_est = covinv_est
  z = 1
  diff = 100
  nplog = obj_L12(Y, mu, p, L, covinv_est, S_bar, graph, lambda1, lambda2, pie) 
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
    mu = rho %*% Y/ (pie * n) # each mu_i is divided by (sum_j tau_ij = pie_i *n)
    temp_memship = apply(rho, 2, which.max)
    #calculate S_bar
    dat1 = list()
    temp_num_cluster = length(unique(temp_memship))
    for(l in 1:temp_num_cluster){
      dat1[[l]] = Y[temp_memship==l,]
    }
    sol_path <- JGL(dat1, penalty = "fused", lambda1, lambda2, 0.1, return.whole.theta=T)
    covinv_est = matrix(0, nrow = p, ncol = temp_num_cluster*p)
    S_bar = matrix(0, p, L*p)
    for (l in 1:L){ # update S_bar
      covinv_est[,((l-1)*p+1):(l*p)] = sol_path$theta[[l]]
      S = matrix(0,p,p) 
      for(j in 1:n){
        S = S + rho[l,j] * outer((Y[j,]-mu[l,]), (Y[j,]-mu[l,]), "*")
      }
      S_bar[,((l-1)*p+1):(l*p)] <- S/(pie[l]*n)
    }
    plog = obj_L12(Y, mu, p, L, covinv_est, S_bar, graph, lambda1, lambda2, pie)  
    #diff = abs(plog - nplog)
    diff = plog-nplog
    nplog = plog
    z = z+1 
  }
  class_memship = apply(rho, 2, which.max)
  output = list(pie = pie, mu = mu, covinv = covinv_est, membership = class_memship, par=c(lambda1, lambda2, ncluster)) 
}
