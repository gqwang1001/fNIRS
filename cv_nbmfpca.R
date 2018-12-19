cv_nbmfpca = function(rejection.index = ind_large, dat = delbox, nsubj = nsubj, nchannels = nchannels, trainSet, tlength = 1){
  
      require(Matrix)
      require(ICSNP)
      require(MASS)  
      
      #### dat -- data input.
      #### subj.selected -- Selected subjects that may be useful in the cross-valadation. 
      #### nsubj -- Number of subjects
      #### nchannels -- Number of channels per subject
      #### rejection.ind -- preselected channels that are not included in the analysis. 
      #### trainSet -- Index set of training data
      
      ind.mat = matrix(1:(nsubj*nchannels), nrow = nsubj, ncol = nchannels, byrow = T)
      train.ind = as.vector(t(ind.mat[trainSet, ])) # indices for training data
      dat.train = dat[train.ind, ]
      dat.test = dat[-train.ind, ]
      
      
      rejection_matrix = matrix(1,nrow = nsubj,ncol = nchannels)
      rejection_matrix[rejection.index] = 0 #0 means reject
      
      rej.mat.train = rejection_matrix[trainSet, ]
      rej.mat.test = rejection_matrix[-trainSet, ]
          
      step1.train = cv.nbmfpca.step1(dat = dat.train, rejection_matrix = rej.mat.train, trainData = T)
      step1.test  = cv.nbmfpca.step1(dat = dat.test, rejection_matrix = rej.mat.test, trainData = F)
      
      score.train = cv.nbmfpca.step2(step1.train.results = step1.train, step1.test.results = step1.test, trainData = T)
      score.test  = cv.nbmfpca.step2(step1.train.results = step1.train, step1.test.results = step1.test, trainData = F)

      dat.results=list(score.train = score.train, score.test = score.test) 
      return(dat.results)
}



cv.nbmfpca.step1 = function(dat = dat, rejection_matrix = rejection_matrix,  trainData = TRUE, tlength = 1){
  #### create functional space using training data, and compute residuals for testing data
  
  
  nsubj = nrow(rejection_matrix)
  nchannels = ncol(rejection_matrix)
  rejection_list = as.vector(t(rejection_matrix))

  dat1 = dat[rejection_list==1, ]
  num_subj = colSums(rejection_matrix)  # number of subjects for each channel
  num_channel = rowSums(rejection_matrix) # number of channels in each subjects
  cum_channel_number = cumsum(num_channel)
  
  ###    Calculate the grand mu
  N = tn = dim(dat)[2]
  mu_subj = matrix(nrow = nsubj, ncol = tn)
  for (i in 1:nsubj) mu_subj[i,] = colMeans(as.matrix(dat1[(cum_channel_number[i]-num_channel[i]+1):cum_channel_number[i],]))
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
  
  resid.temp2 = resid.temp
  resid.temp2[rejection_list==0,] = rep(0,tn) # set the rejected observations to 0s

  if (trainData == TRUE){
      ### compute for training data only
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
        
        min.K1 = 4
        min.K2 = 4
        
        K1 = max( which(cumsum(percent1) < 0.9 | percent1 > 1/tn ) + 1, min.K1 )
        K2 = max( which(cumsum(percent2) < 0.9 | percent2 > 1/tn ) + 1, min.K2 )
        
        fpca1.vectors = e1$vectors[, 1:K1]*sqrt(N/tlength)
        fpca2.vectors = e2$vectors[, 1:K2]*sqrt(N/tlength)
        
        
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
        cross.integral = t(fpca1.vectors)%*%fpca2.vectors*tlength/N

        return(list(fpca1.vectors = fpca1.vectors, fpca2.vectors = fpca2.vectors, cross.integral = cross.integral, 
                    resid = resid.temp2, K1 = K1, K2 = K2, 
                    nsubj = nsubj, nchannels = nchannels, num_channel = num_channel,
                    tlength = tlength, N = N,
                    rejection_matrix = rejection_matrix))
        
  }else{
        return(list(resid = resid.temp2,nsubj = nsubj, nchannels = nchannels,
                    tlength = tlength, N = N,num_channel = num_channel,
                    rejection_matrix = rejection_matrix))
    }
  
}

cv.nbmfpca.step2 = function(step1.train.results, step1.test.results, rejection_matrix, trainData = T){
    ##### project the training and testing data onto the functional space generated by training data

    fpca1.vectors = step1.train.results$fpca1.vectors
    fpca2.vectors = step1.train.results$fpca2.vectors
    cross.integral = step1.train.results$cross.integral
    
    K1 = step1.train.results$K1
    K2 = step1.train.results$K2
    tlength = step1.train.results$tlength
    N = step1.train.results$N
    if (trainData == T) {
      resid.temp2 = step1.train.results$resid
      M = step1.train.results$nsubj
      J = step1.train.results$nchannels 
      rejection.index = which(step1.train.results$rejection_matrix==0, arr.ind = T)
      num_channel = step1.train.results$num_channel
      
    } else {
      resid.temp2 = step1.test.results$resid
      M = step1.test.results$nsubj
      J = step1.test.results$nchannels 
      rejection.index = which(step1.test.results$rejection_matrix==0, arr.ind = T)
      num_channel = step1.train.results$num_channel
      
    }
    
  int1 = matrix(0, M*J, K1)
  int2 = matrix(0, M*J, K2)
  for(i in 1:(M*J))   {
    for(j in 1:K1) int1[ i ,j] = sum( resid.temp2[i,] * fpca1.vectors[,j] ) * tlength /N
    for(j in 1:K2) int2[ i ,j] = sum( resid.temp2[i,] * fpca2.vectors[,j] ) * tlength /N  
  }
  
  s1 = matrix(NA, M*J, K1)
  s2 = matrix(NA, M*J, K2)
  
  design.xi = ginv( diag(rep(1,K1)) - cross.integral %*% t(cross.integral) )
  for(m in 1:M) {
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
  
  s1.reshape = array(NA, c(M,J,K1))
  s2.reshape = array(NA, c(M,J,K2))
  
  s1.reshape1 = matrix(nrow = M,ncol = K1)
  for(i in 1:K1){ 
    s1.reshape[,,i] = fnirs.reshape(s1[,i],nsubj = M,nchannels = J)
    s1.reshape1[,i] = s1.reshape[,1,i]
  }
  
  for(i in 1:K2) {
    s.temp = fnirs.reshape(s2[,i],nsubj = M, nchannels = J)
    s.temp[rejection.index]=NA;
    s2.reshape[,,i]=s.temp
  }

  return(list(scores1 = s1.reshape1, scores2 = s2.reshape))
  
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

