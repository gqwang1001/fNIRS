fnirs.reshape = function(delbox_intg,nsubj,nchannels){
  delbox_intg_reshape = matrix(NA,nrow = nsubj,ncol = nchannels)
  for(i in 1:nsubj){
    for (j in 1:nchannels){
      delbox_intg_reshape[i,j] = delbox_intg[(i-1)*nchannels+j]
    }
  }
  return(delbox_intg_reshape)
}
