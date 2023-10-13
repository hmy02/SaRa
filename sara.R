sara=function(h,y){
  Y=c(rep(0,h),y,rep(0,h-1))
  D_x=c()
  for (x in h+1:(length(Y)-h)) {
    b=(x-h):(x+h-1) #b是6个数
    D=abs(mean(Y[b[1:h]])-mean(Y[b[(h+1):(2*h)]]))
    D_x=c(D_x,D)
  }
  D_x=c(rep(0,h),D_x,rep(0,h))
  index = c()
  for(i in (h+1):(length(D_x)-h)){
    f=max(D_x[(i-h):(i+h)])
    xx=which((D_x==f)&(c(rep(F,i-h-1),rep(T,2*h+1),rep(F,length(D_x)-i-h))))
    index = c(index, xx[1])
  }
  index = unique(index)
  lambda = 2*sqrt(2/h)*sd(y)
  site = index[which(D_x[index]>=lambda)] - h
  return(site)
}


multi_sara<-function(y,h_k){
  # 构建候选者池
  candidate=c()
  for(h in h_k){
    candidate=c(candidate,sara(h,y))
  }
  candidate_pool=sort(unique(candidate))
  sigmaJ2=var(y[candidate_pool])
  last=Inf
  # n/2*log(σ)+Jlog(n)不再减小时停止
  while((length(y)/2*log(sigmaJ2)+length(candidate_pool)*log(length(y)))<last){
    last=length(y)/2*log(sigmaJ2)+length(candidate_pool)*log(length(y))
    # 寻找使sigma变大最小的(least increase)
    min_index=candidate_pool[1]
    for(j in candidate_pool){
      if((sd(y)-sd(y[-j]))<(sd(y)-sd(y[-min_index]))){
        min_index=j
      }
    }
    # 剔除该项并更新数据
    candidate_pool=candidate_pool[candidate_pool!=min_index]
    sigmaJ2=var(y[candidate_pool])
  }
  return(candidate_pool)
}
