


simulated_annealing =	function(epochs = 3, estat = 3, family = gaussian(), fparam,isobsbinary, MIC = stats::BIC,SIC =stats::AIC, fobserved, ranges, ns,initpr, data,a.prior.em, b.prior.em, a.post.em, wcrit, a.Q.prior, b.Q.prior, X.min, X.max,
                               a.maxit.prior=20,b.maxit.prior=0.1, prior.inclusion, s2.post.maxit,col.rate.maxit,tmin,temp,dt,p_id,seed) {
    
  if(family$family == "binomial")
  {
    osv = paste("cbind(",fobserved[1],",",fobserved[2],")")
  }
  else
    osv  = fobserved
  
  # initialize the set of parameters
  maxit=a.maxit.prior
  stm = proc.time()
  set.seed(seed, kind = NULL, normal.kind = NULL)
  print(paste("begin SA methasearch procedure","with EM local search greedy updates for each neighbourhood!"))
  tcur=temp
  n.obs.em = 0
  sum.obs.em =0
  n.obs.maxit = 0
  sum.obs.maxit =0
  n.obs.k = 0
  sum.obs.k =0
  post.inclusion =  prior.inclusion 
  
  tid=0
  ii=1
  
  globcrit=array(2,data = Inf)
  
  fm = NULL
  varcur.obs = rbinom(n = length(fparam),size = 1,prob = 0.5)
  varcur.trs = rbinom(n = length(fparam),size = 1,prob = 0.5)
  
  N=length(varcur.obs)
  
  covobs = fparam[which(varcur.obs == 1)]
  covtrans = fparam[which(varcur.trs == 1)]
  obsconst=1
  transconst=1
  formula = as.formula(paste(osv, " ~ ",obsconst,ifelse(test = length(covobs)>0, " + ",""), paste(covobs, collapse=" + ")))
  transition =as.formula(paste(" ~ ",transconst,ifelse(test = length(covtrans)>0, " + ",""), paste(covtrans, collapse=" + ")))
  if(length(covtrans)>0 && length(covobs)>0)
  {
    
    formula = as.formula(paste(osv, " ~ ","1 + ", paste(covobs, collapse=" + ")))
    
    transition =as.formula(paste(" ~ ","1 + ", paste(covtrans, collapse=" + ")))
    
    mod = depmix(response = formula,data=data,transition = transition,nstates=ns, instart = initpr,family=family)
    
  }else 
  {
    if(length(covtrans)>0 && length(covobs)==0)
    {
      transition =as.formula(paste(" ~ ","1 + ", paste(covtrans, collapse=" + ")))
      mod = depmix(response = as.formula(paste(osv, " ~ ","1")) ,data=data,transition = transition,nstates=ns, instart = initpr,family=family)
    }
    
    if(length(covtrans)==0 && length(covobs)>0)
    {
      
      formula = as.formula(paste(osv, " ~ ","1 + ", paste(covobs, collapse=" + ")))
      mod = depmix(response = formula,data=data,nstates=ns, instart = initpr,family=family)
    }
    if(length(covtrans)==0 && length(covobs)==0)
    {
      mod = depmix(response = as.formula(paste(osv, " ~ ","1")),data=data,nstates=ns, instart = initpr,family=family)
    }
  }
  curn=rnbinom(n = 1,size = a.prior.em + sum.obs.em,prob = 1-1/(1+b.prior.em + n.obs.em))
  
  #print(curn) 
  #rgamma(n = 1,shape = (a.prior.em + n.obs.em*a.post.em),rate = (b.prior.em + sum.obs.em))
  em.params=em.control(maxit = curn,random.start = TRUE)
  
  #fm=fit(mod,emcontrol = em.params)
  
  capture.output({withRestarts(tryCatch(capture.output({fm=fit(mod,emcontrol = em.params)}),error = function(e){invisible(e)}), 
                               abort = function(){ fm=mod;  onerr=TRUE})}) # fit the modal, get local improvements
  
  
  globcrit[1]=MIC(fm)
  globcrit[2]=SIC(fm)
  varglob.obs=varcur.obs
  varglob.trs=varcur.trs
  glmod=fm
  
  curcrit=globcrit
  
  for(epoch in 1:epochs)
  {
    print(paste0("start epoch ",epoch))
    tcur = temp
    while(tcur>=tmin)
    {
      if(tid %% p_id == 0){
        
        n.obs.maxit = 0
        sum.obs.maxit =0
        maxit = rnbinom(n = 1,size = a.maxit.prior +  sum.obs.maxit,prob = 1-1/(1+b.maxit.prior + n.obs.maxit))
        #rmutil::rbetabinom(n=1,)
        #sigma =1/(1/(s2.prior.maxit) + n.obs.maxit/s2.post.maxit)
        
        print(paste(ii," iterations completed up to now"))
        print(paste("MIC.cur.glob = ", globcrit[1]))
        print(paste("aneal to temp = ",tcur," proceed with this temp. for",
                    maxit, "iterations"))
        
      }
      for(i in 1:maxit)
      {
        withRestarts(tryCatch({
          
          # change the neighbourhood i.e. randomly change two positions in the given vector
          varcurb.obs=varcur.obs
          varcurb.trs=varcur.trs
          
          K=rbinom(n = 1,size = N,prob = rbeta(n = 1,shape1 = (a.Q.prior + sum.obs.k),shape2 = (b.Q.prior + N*n.obs.k - sum.obs.k)))
          
          
          
          
          chids = sample.int(n=length(varcur.obs),size = K,replace = F)
          varcur.obs[chids] = rbinom(n = K,size = 1,prob = post.inclusion[chids,1]/sum(post.inclusion[,1]))
          varcur.trs[chids] = rbinom(n = K,size = 1,prob = post.inclusion[chids,2]/sum(post.inclusion[,2]))
          
          covobs = fparam[which(varcur.obs == 1)]
          covtrans = fparam[which(varcur.trs == 1)]
          obsconst=1
          transconst=1
          formula = as.formula(paste(osv, " ~ ",obsconst,ifelse(test = length(covobs)>0, " + ",""), paste(covobs, collapse=" + ")))
          transition =as.formula(paste(" ~ ",transconst,ifelse(test = length(covtrans)>0, " + ",""), paste(covtrans, collapse=" + ")))
          if(length(covtrans)>0 && length(covobs)>0)
          {  
            
            formula = as.formula(paste(osv, " ~ ","1 + ", paste(covobs, collapse=" + ")))
            transition =as.formula(paste(" ~ ","1 + ", paste(covtrans, collapse=" + ")))
            
            mod = depmix(response = formula,data=data,transition = transition,nstates=ns, instart = initpr,family=family)
            
          }else 
          {
            if(length(covtrans)>0 && length(covobs)==0)
            {
              
              transition =as.formula(paste(" ~ ","1 + ", paste(covtrans, collapse=" + ")))
              mod = depmix(response = as.formula(paste(osv, " ~ ","1")),data=data,transition = transition,nstates=ns, instart = initpr,family=family)
            }
            
            if(length(covtrans)==0 && length(covobs)>0)
            {
              
              formula = as.formula(paste(osv, " ~ ","1 + ", paste(covobs, collapse=" + ")))
              
              mod = depmix(response = formula,data=data,nstates=ns, instart = initpr,family=family)
            }
            if(length(covtrans)==0 && length(covobs)==0)
            {
              mod = depmix(response = as.formula(paste(osv, " ~ ","1")),data=data,nstates=ns, instart = initpr,family=family)
              
            }
          }
          
          # set initial values for the parameters with respect to a given neighbourhood
          
          
          fmb=fm
          onerr=FALSE
          curn= rnbinom(n = 1,size = a.prior.em + sum.obs.em,prob = 1-1/(1+b.prior.em + n.obs.em))
          #print(curn)
          #  rgamma(n = 1,shape = (a.prior.em + n.obs.em*a.post.em),rate = (b.prior.em + sum.obs.em))
          em.params=em.control(maxit = curn,random.start = TRUE)
          
          #finally = print(paste("loop body finished in iteration ",ii))
          # optimize in the new neighbourhood
          capture.output({withRestarts(tryCatch(capture.output({fm=fit(mod,emcontrol = em.params)}),error = function(e){invisible(e)}), 
                                       abort = function(){varcur.obs=varcurb.obs; varcur.trs=varcurb.trs; fm=fmb;  onerr=TRUE})}) # fit the modal, get local improvements
          
          ii=ii+1
          
          #update local values
          
          
          if(MIC(fm)<globcrit[1])
          {
            print(paste("update global optimum with new global MIC = ", MIC(fm)))
            
            globcrit[1]=MIC(fm)
            globcrit[2]=SIC(fm)
            glmod=fm
            varglob.obs=varcur.obs
            varglob.trs=varcur.trs
            if(epoch<estat)
            {
              n.obs.em=(n.obs.em + 1)
              sum.obs.em=(sum.obs.em + curn)
              sum.obs.k=sum.obs.k + K
              n.obs.k=n.obs.k + 1
              n.obs.maxit = (n.obs.maxit +1)
              sum.obs.maxit = sum.obs.maxit + maxit
              post.inclusion[,1] = prior.inclusion[,1] + varglob.obs
              post.inclusion[,2] = prior.inclusion[,2] + varglob.trs
            }
            
          }
          if(((MIC(fm)-curcrit[1]))-(tcur)<=log(runif(n = 1,min = 0,max = 1)) && !onerr)
          {
            
            curcrit[1]=MIC(fm)
            curcrit[2]=SIC(fm)
            
            if(epoch<estat && tcur < 1)
            {
              n.obs.em=(n.obs.em + 1)
              sum.obs.em=(sum.obs.em + curn)
              sum.obs.k=sum.obs.k + K
              n.obs.k=n.obs.k + 1
              n.obs.maxit = (n.obs.maxit +1)
              sum.obs.maxit = sum.obs.maxit + maxit
              post.inclusion[,1] = prior.inclusion[,1] + varcur.obs
              post.inclusion[,2] = prior.inclusion[,2] + varcur.trs
            }
            
            
            
          }else
          {
            varcur.obs=varcurb.obs
            varcur.trs=varcurb.trs
          }
          
        }),abort = function(){varcur.obs=varcurb.obs; varcur.trs=varcurb.trs; fm=fmb;  onerr=TRUE})
      }
      # update temperature
      tcur=tcur*exp(-dt)
      tid=tid+1
      n.obs.k = 0
      sum.obs.k =0
    }
    print(paste(ii," iterations completed in total"))
    print(paste(tcur," is the final temperature"))
    
    print(paste("MIC.cur.glob = ", globcrit[1]))
  }
  
  etm = proc.time()
  tt=(etm[3]-stm[3])*1000
  
  return(list(model=glmod, vars.obs=varglob.obs,vars.trs=varglob.trs , time=(tt)))
	}


