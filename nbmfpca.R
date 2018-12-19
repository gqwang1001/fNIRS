nbmfpca = function(rejection.index = ind_large, dat = delbox, nsubj = nsubj, nchannels = nchannels, tlength = 1){
      
      require(Matrix)
      require(ICSNP)
      require(MASS)  
    
      #### dat -- data input.
      #### subj.selected -- Selected subjects that may be useful in the cross-valadation. 
      #### nsubj -- Number of subjects
      #### nchannels -- Number of channels per subject
      #### rejection.ind -- preselected channels that are not included in the analysis. 
  
  
      rejection_matrix = matrix(1,nrow = nsubj,ncol = nchannels)
      rejection_matrix[rejection.index] = 0 #0 means reject
      rejection_list = as.vector(t(rejection_matrix))
      
      
      dat1 = dat[rejection_list==1,]
      num_subj = colSums(rejection_matrix)  # number of subjects for each channel
      num_channel = rowSums(rejection_matrix) #number of channels in each subjects
      cum_channel_number = cumsum(num_channel)
      
      ###    Calculate the grand mu
      tn = dim(dat)[2]
      mu_subj = matrix(nrow = nsubj,ncol = tn)
      for (i in 1:nsubj) mu_subj[i,] = colMeans(dat1[(cum_channel_number[i]-num_channel[i]+1):cum_channel_number[i],])
      mu = colMeans(mu_subj)
      
      dat2 = dat
      dat2[rejection_list==0,] = rep(0,tn)
      ###   Calculate the channel specific mean
      eta = matrix(0,nchannels,tn)
      for (j in 1:nchannels) eta[j,] = colSums(dat2[(0:(nsubj-1)*nchannels) + j, ])/num_subj[j] - mu
      
      ###   Calculate the resid
      resid.temp = matrix(0,nrow = nsubj*nchannels,ncol=tn)
      for(j in 1:nchannels) resid.temp[(0:(nsubj-1)*nchannels) + j,] = t(t(dat[ (0:(nsubj-1)*nchannels)+j,])-eta[j,]-mu)
      resid = resid.temp[rejection_list==1, ] #remove rejected channels
      
      ###     Estimate the three covariance functions: overall covariance G, 
      ###     between covariance Gb and within covariance Gw
      resid.temp.subj = list()
      for(m in 1:nsubj) resid.temp.subj[[m]] = resid[(cum_channel_number[m] - num_channel[m] + 1) : cum_channel_number[m], ]
      
      G.temp = matrix(0, tn, tn)
      for (m in 1: nsubj) G.temp = G.temp + crossprod(resid.temp.subj[[m]])/num_channel[m]
      G = G.temp/nsubj
      
      Gw.temp =  matrix(0, tn, tn)
      for (m in 1: nsubj){
        if(num_channel[m] >1){                                          
          resid.pair.diff =ICSNP::pair.diff(resid.temp.subj[[m]]) # compute pairwise (j1<j2) difference 
          Gw.temp = Gw.temp + crossprod(resid.pair.diff)/nrow(resid.pair.diff)
        }
      }
      
      Gw = Gw.temp/(2*sum(num_channel > 1))
      Gb = G - Gw 
      
      e1 = eigen(Gb)
      e2 = eigen(Gw)
      
      ### Output the eigenvalues
      
      N = tn
      fpca1.value = e1$values* tlength / N
      fpca2.value = e2$values* tlength / N
      
      ###     Keep only non-negative eigenvalues
      fpca1.value = ifelse(fpca1.value>=0, fpca1.value, 0)
      fpca2.value = ifelse(fpca2.value>=0, fpca2.value, 0)
      
      ###     Calculate the percentage of variance that are explained by the components
      
      percent1 = (fpca1.value)/sum(fpca1.value)
      percent2 = (fpca2.value)/sum(fpca2.value)
      
      rho = sum(fpca1.value)/(sum(fpca1.value)+sum(fpca2.value))
      
      ###     Decide the number of components that are kept at level 1 and 2. The general
      ###     rule is to stop at the component where the cumulative percentage of variance 
      ###     explained is greater than 90% and the variance explained by any single component
      ###     after is less than 1/N. The number of components are also no less than the 
      ###     pre-determined minimum values min.K1 or min.K2.
      min.K1 = 4
      min.K2 = 4
      K1 = max( which(cumsum(percent1) < 0.9 | percent1 > 1/tn ) + 1, min.K1 )
      K2 = max( which(cumsum(percent2) < 0.9 | percent2 > 1/tn ) + 1, min.K2 )
      
      ###     estimate eigen vectors for discretized covariance matrices and
      ###     transform them into norm one eigenfunctions
      fpca1.vectors = e1$vectors[, 1:K1]*sqrt(N/tlength)
      fpca2.vectors = e2$vectors[, 1:K2]*sqrt(N/tlength)
      
      ###     The eigenfunctions are unique only up to a change of signs.
      ###     Select the signs of eigenfunctions so that the integration over the domain 
      ###     is non-negative
      
      for(i in 1:K1) {
        v2 = fpca1.vectors[,i]
        tempsign = sum(v2)
        fpca1.vectors[,i] = ifelse(tempsign<0, -1,1) * v2
      }
      for(i in 1:K2) {
        v2 = fpca2.vectors[,i]
        tempsign = sum(v2)
        fpca2.vectors[,i] = ifelse(tempsign<0, -1,1) * v2
      }
      
      ###     First, calculate the inner product (the cosine of the angles) between 
      ###     level 1 eigenfunctions and level 2 eigenfunctions
      cross.integral = t(fpca1.vectors)%*%fpca2.vectors*tlength/N
      ###     Next, calculate the inner product of each centered function with the 
      ###     level 1 or level 2 eigenfunctions
      
      resid.temp2 = resid.temp
      resid.temp2[rejection_list==0,] = rep(0,tn) # set the rejected observations to 0s
      
      ###     Next, calculate the normalized inner product of each centered function with the 
      ###     level 1 or level 2 eigenfunctions
      
      
      M= nsubj
      J = nchannels 
      int1 = matrix(0, M*J, K1)
      int2 = matrix(0, M*J, K2)
      for(i in 1:(M*J))   {
        for(j in 1:K1) int1[ i ,j] = sum( resid.temp2[i,] * fpca1.vectors[,j] ) * tlength /N
        for(j in 1:K2) int2[ i ,j] = sum( resid.temp2[i,] * fpca2.vectors[,j] ) * tlength /N  
      }
      
      
      ###     Finally, calculate the principal component scores based on the formulas
      ###     given in the paper. This is an alternative (much more efficient) way to implement the BLUP formula. 
      s1 = matrix(NA, nsubj*nchannels, K1)
      s2 = matrix(NA, nsubj*nchannels, K2)
      
     # source("fnirs.reshape.R")
      design.xi = ginv( diag(rep(1,K1)) - cross.integral %*% t(cross.integral) )
      for(m in 1:nsubj) {
        resd = rep(0, K1)
        for(j in 1:J) {
          index = (m-1) * J + j
          resd = resd + ( int1[index,] - drop(cross.integral %*% int2[index,]) )/num_channel[m]###########################
        }
        index.m = ( (m-1) * J + 1 ) : (m*J)
        xi.temp = design.xi %*% resd
        s1[index.m,] = matrix(rep(xi.temp, each=J), nrow=J)
        s2[index.m,] = t(t(int2[index.m,]) - drop( t(cross.integral) %*% xi.temp ))
      }
      
      s1.reshape = array(NA, c(nsubj,nchannels,K1))
      s2.reshape = array(NA, c(nsubj,nchannels,K2))

      s1.reshape1 = matrix(nrow = nsubj,ncol = K1)
      for(i in 1:K1){ 
        s1.reshape[,,i] = fnirs.reshape(s1[,i],nsubj = nsubj,nchannels = nchannels)
        s1.reshape1[,i] = s1.reshape[,1,i]
      }
      
      for(i in 1:K2) {
        s.temp = fnirs.reshape(s2[,i],nsubj = nsubj, nchannels = nchannels)
        s.temp[rejection.index]=NA;
        s2.reshape[,,i]=s.temp
      }
      
      # 
      # 
      # I = nrow(resid)
      # score1.blup = matrix(0, nrow = M, ncol = K1)
      # score2.blup = matrix(0, nrow = I, ncol = K2)
      # row.ind = 0
      # 
      # for(m in 1:M){
      #   Y.tilde.m = resid.temp.subj[[m]]
      #   Ji = nrow(Y.tilde.m)
      #   
      #   A11 = diag(K1)*Ji
      #   A12 = t(as.vector(rep(1,Ji))) %x% cross.integral
      #   A21 = t(A12)
      #  # A22 = diag(rep(1,Ji)) %x% diag(rep(1, K2))
      #   A22 = diag(Ji*K2)
      #     
      #     
      #   A = rbind(cbind(A11,A12),cbind(A21,A22))
      #   Aiv = chol2inv(chol(A))
      #   
      #   B1 = t(as.vector(rep(1,Ji)))%x% t(fpca1.vectors)
      #   B2 = diag(Ji) %x% t(fpca2.vectors)
      #   B = rbind(B1,B2)
      #   
      #   Zmat = Aiv %*% B%*% c(t(Y.tilde.m)) / N
      #   #Zmat = solve(A, B %*% c(t(Y.tilde.m)))
      #   
      #   sc1 = Zmat[1:K1,]
      #   sc2 = Zmat[-c(1:K1)]
      #   
      #   score1.blup[m,] = sc1
      #   subj.indices = row.ind + 1:Ji
      #   score2.blup[subj.indices,] = matrix(sc2, nrow = Ji, byrow = T)
      #   
      #   
      #   row.ind = row.ind + num_channel[m]
      # }
      
      
      
      
      dat.results=list(K1=K1,K2=K2,lambda1=fpca1.value,lambda2=fpca2.value,
                          phi1=fpca1.vectors,phi2=fpca2.vectors,
                          rho = rho,scores1=s1.reshape1,
                          scores2=s2.reshape,mu=mu,eta=t(eta),
                          reject.index=rejection.index) 
      return(dat.results)
}



fnirs.reshape = function(dat,nsubj,nchannels){
  dat_reshape = matrix(NA,nrow = nsubj,ncol = nchannels)
  for(i in 1:nsubj){
    for (j in 1:nchannels){
      dat_reshape[i,j] = dat[(i-1)*nchannels+j]
    }
  }
  return(dat_reshape)
}

