
select_depmix =	function(data,epochs = 3, estat = 3, MIC = stats::AIC,SIC =stats::BIC,family = gaussian(),fparam,isobsbinary,prior.inclusion,fobserved = fobserved,ranges = 1,ns =3,initpr = c(0,1,0),a.prior.em = list(min = 100,max = 200), b.prior.em = list(min = 1,max = 10), a.post.em = 10,wcrit = 1,a.Q.prior = list(min = 5, max = 5), b.Q.prior = list(min =15,max = 15), X.min = 3, X.max = 10, a.maxit.prior  = list(min = 15, max = 35) ,b.maxit.prior  = list(min = 2, max = 2) ,s2.post.maxit = 3 ,col.rate.maxit = list(min = 0.97, max = 0.99),tmin = list(min = 0.0000005, max = 0.05),temp = list(min = 20000,max = 2000000000),dt = list(min = 2, max = 6) ,p_id = 1,seeds = runif(4,1,1000),cores = 4) {

  params = NULL
  for(i in 1:cores)
  {
      params[[i]]=list(epochs = epochs, estat = estat,MIC = MIC, SIC = SIC, fparam = fparam,family = family,isobsbinary = isobsbinary,prior.inclusion = prior.inclusion,fobserved = fobserved,ranges = ranges,ns = ns,initpr = initpr,data = data,a.prior.em = runif(n = 1,a.prior.em$min,a.prior.em$max), b.prior.em = runif(n = 1,b.prior.em$min,b.prior.em$max), a.post.em = a.post.em,wcrit = wcrit,a.Q.prior = runif(n = 1, a.Q.prior$min, a.Q.prior$max), b.Q.prior = runif(n = 1, b.Q.prior$min, b.Q.prior$max), X.min = X.min, X.max = X.max,a.maxit.prior  = runif(n = 1,a.maxit.prior$min,a.maxit.prior$max) ,b.maxit.prior  =  runif(n = 1,b.maxit.prior$min,b.maxit.prior$max) ,s2.post.maxit = s2.post.maxit ,col.rate.maxit = runif(n=1,col.rate.maxit$min,col.rate.maxit$max),tmin = runif(n = 1,tmin$min, tmin$max),temp = runif(n = 1,temp$min,temp$max),dt = runif(n = 1,dt$min, dt$max) ,p_id = p_id,seed = seeds[i])
      
      #params[[i]]=list(MIC = stats::AIC,SIC =stats::BIC,fparam = fparam,isobsbinary = c(0,0,rep(1,length(fparam))),prior.inclusion = array(1,c(length(fparam),2)),fobserved = fobserved,ranges = 1,ns = ns,initpr = c(0,1,0),data = X,a.prior.em = runif(n = 1,min = 100,max = 200), b.prior.em = runif(n = 1,min = 1,max = 10), a.post.em = 10,wcrit = 1,a.Q.prior = 5, b.Q.prior = 15, X.min = 3, X.max = 10,a.maxit.prior  = runif(n = 1,min = 15, max = 35) ,b.maxit.prior  = 2 ,s2.post.maxit = 3 ,col.rate.maxit = runif(n=1,min = 0.97, max = 0.99),tmin = runif(n = 1,min = 0.0000005, max = 0.05),temp = runif(n = 1,min = 20000,max = 2000000000),dt = runif(n = 1,min = 2, max = 6) ,p_id = 1,seed = runif(1,1,1000))
  }
  
  results = mclapply(FUN = function(x) do.call(depmixS4pp::simulated_annealing,x),X = params,mc.cores = cores)
  
  
 
  best.mic = -1
  best.sic = -1
  
  res.aic = NULL
  res.bic = NULL
  for(i in 1:M)
  {
    if(length(results[[i]])>1)
    { 
      if(best.mic ==-1)
      {
        best.mic = i
        best.sic = i
      }
      res.aic = c(res.aic,AIC(results[[i]]$model))
      res.bic = c(res.bic,BIC(results[[i]]$model))
      if(MIC(results[[i]]$model)<MIC(results[[best.mic]]$model))
        best.mic = i 
      if(SIC(results[[i]]$model)<SIC(results[[best.sic]]$model))
        best.sic = i 
    }
    else
    {
      res.aic =c(res.aic,Inf)
      res.bic = c(res.bic,Inf)
    }
  }
  
  return(list(best.mic=best.mic,best.sic = best.sic, bics = res.bic, aics = res.aic, results = results))
  
}


