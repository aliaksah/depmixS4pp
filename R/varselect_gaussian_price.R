
library("depmixS4")

library(parallel)

simulated_annealing_gaussian=function(epochs = 3, estat = 3, fparam,isobsbinary, MIC = stats::BIC,SIC =stats::AIC, fobserved, ranges, ns,initpr, data,a.prior.em, b.prior.em, a.post.em, wcrit, a.Q.prior, b.Q.prior, X.min, X.max,
                                      a.maxit.prior=20,b.maxit.prior=0.1, prior.inclusion =array(1,dim = c(length(c(0,0,1,1,1,1,1,1,1,1,1,1)),2)), s2.post.maxit,col.rate.maxit,tmin,temp,dt,p_id,seed){
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
  formula = as.formula(paste(fobserved, " ~ ",obsconst,ifelse(test = length(covobs)>0, " + ",""), paste(covobs, collapse=" + ")))
  transition =as.formula(paste(" ~ ",transconst,ifelse(test = length(covtrans)>0, " + ",""), paste(covtrans, collapse=" + ")))
  if(length(covtrans)>0 && length(covobs)>0)
  {
    
    formula = as.formula(paste(fobserved, " ~ ","1 + ", paste(covobs, collapse=" + ")))
    
    transition =as.formula(paste(" ~ ","1 + ", paste(covtrans, collapse=" + ")))
    
    mod = depmix(response = formula,data=data,transition = transition,nstates=ns, instart = initpr)
                
  }else 
  {
    if(length(covtrans)>0 && length(covobs)==0)
    {
      transition =as.formula(paste(" ~ ","1 + ", paste(covtrans, collapse=" + ")))
      mod = depmix(response = as.formula(paste(fobserved, " ~ ","1")) ,data=data,transition = transition,nstates=ns, instart = initpr)
    }
    
    if(length(covtrans)==0 && length(covobs)>0)
    {
      
      formula = as.formula(paste(fobserved, " ~ ","1 + ", paste(covobs, collapse=" + ")))
      mod = depmix(response = formula,data=data,nstates=ns, instart = initpr)
    }
    if(length(covtrans)==0 && length(covobs)==0)
    {
      mod = depmix(response = as.formula(paste(fobserved, " ~ ","1")),data=data,nstates=ns, instart = initpr)
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
          formula = as.formula(paste(fobserved, " ~ ",obsconst,ifelse(test = length(covobs)>0, " + ",""), paste(covobs, collapse=" + ")))
          transition =as.formula(paste(" ~ ",transconst,ifelse(test = length(covtrans)>0, " + ",""), paste(covtrans, collapse=" + ")))
          if(length(covtrans)>0 && length(covobs)>0)
          {  
            
            formula = as.formula(paste(fobserved, " ~ ","1 + ", paste(covobs, collapse=" + ")))
            transition =as.formula(paste(" ~ ","1 + ", paste(covtrans, collapse=" + ")))
            
            mod = depmix(response = formula,data=data,transition = transition,nstates=ns, instart = initpr)
                         
          }else 
          {
            if(length(covtrans)>0 && length(covobs)==0)
            {
              
              transition =as.formula(paste(" ~ ","1 + ", paste(covtrans, collapse=" + ")))
              mod = depmix(response = as.formula(paste(fobserved, " ~ ","1")),data=data,transition = transition,nstates=ns, instart = initpr)
            }
            
            if(length(covtrans)==0 && length(covobs)>0)
            {
              
              formula = as.formula(paste(fobserved, " ~ ","1 + ", paste(covobs, collapse=" + ")))
              
              mod = depmix(response = formula,data=data,nstates=ns, instart = initpr)
            }
            if(length(covtrans)==0 && length(covobs)==0)
            {
              mod = depmix(response = as.formula(paste(fobserved, " ~ ","1")),data=data,nstates=ns, instart = initpr)
              
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

sp500data = read.table("file:///nr/samba/user/ahu/depmixS4/data/all_stocks_5yr.csv",sep = ",",header = T)

name = "AMZN"
tmp = sp500data[which(sp500data$Name==name),c(1,5)]
tmp[is.na(tmp)]=0
tmp[,2]=tmp[,2]# c(0,diff(log(tmp[,2]), lag=1))
names(tmp)[2] = name
X = tmp
for(name in unique(sp500data$Name))
{
  if(name == "AMZN")
    next
  tmp = sp500data[which(sp500data$Name==name),c(1,5)]
  tmp[is.na(tmp)]=0
  tmp[,2]= tmp[,2]#c(0,diff(log(tmp[,2]), lag=1))
  names(tmp)[2] = name
  X = merge(X,tmp,by = "date",all.x = T) 
}

X = X[-1,]

X[is.na(X)] = 0


cors = cor(X[,-1])
slectid = c(1,which(abs(cors[1,])>0.95)+1)

X = X[,slectid]

X.test = X[983:1258,]
X = X[1:982,]

rm(cors)
rm(tmp)
rm(sp500data)
gc()

fparam = colnames(X)[3:length(colnames(X))]
fobserved = colnames(X)[2]


ns = 3
M=50
params = NULL
for(i in 1:M)
{
  withRestarts(tryCatch({
    params[[i]]=list(MIC = stats::AIC,SIC =stats::BIC,fparam = fparam,isobsbinary = c(0,0,rep(1,length(fparam))),prior.inclusion = array(1,c(length(fparam),2)),fobserved = fobserved,ranges = 1,ns = ns,initpr = c(0,1,0),data = X,a.prior.em = runif(n = 1,min = 100,max = 200), b.prior.em = runif(n = 1,min = 1,max = 10), a.post.em = 10,wcrit = 1,a.Q.prior = 5, b.Q.prior = 15, X.min = 3, X.max = 10,a.maxit.prior  = runif(n = 1,min = 15, max = 35) ,b.maxit.prior  = 2 ,s2.post.maxit = 3 ,col.rate.maxit = runif(n=1,min = 0.97, max = 0.99),tmin = runif(n = 1,min = 0.0000005, max = 0.05),temp = runif(n = 1,min = 20000,max = 2000000000),dt = runif(n = 1,min = 2, max = 6) ,p_id = 1,seed = runif(1,1,1000))
  }),abort = function(){onerr=TRUE})
}

results = mclapply(FUN = function(x) do.call(simulated_annealing_gaussian,x),X = params)


res.aic = NULL
res.bic = NULL
best = 1
for(i in 1:M)
{
  if(length(results[[i]])>1)
  { 
    res.aic = c(res.aic,AIC(results[[best]]$model))
    res.bic = c(res.bic,BIC(results[[best]]$model))
    if(AIC(results[[i]]$model)<AIC(results[[best]]$model))
      best = i 
  }
}

ids = sample.int(50,50,F)
png(file=paste("modsFIN.png",sep = ""),width     = 10,
    height    = 5,
    units     = "in",
    res       = 500)
plot(res.bic[ids],col = "red",ylim = c(min(c(res.aic,res.bic)),max(c(res.aic,res.bic))),type = "p",ylab = "Criterion",xlab = "Thread")
lines(rep(min(res.bic),M),col = "red")
points(res.aic[ids],col = "blue",type = "p")
lines(rep(min(res.aic),M),col = "blue")
dev.off()


results[[best]]$model
write.table(depmixS4::getpars(results[[best]]$model),"bestmodel1")

depmixS4::getmodel(results[[best]]$model)

posts = depmixS4::posterior(results[[best]]$model)

plot(X$AMZN,type = "p", col = "red")
#plot(exp(cumsum(X$AMZN)))
lines(posts[,1],col = "red")

summary(results[[best]]$model)


depmixS4::getpars(results[[best]]$model)


save.image("d:/filenameal.RData")