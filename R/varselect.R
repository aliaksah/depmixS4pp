library("depmixS4")
library(RCurl)
library(parallel)

rm(list = ls(all = TRUE))


M=5
size=1

data.example2 = read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/Epigenetic%20Data/epigen.txt"),sep = ",",header = T)[,2:30]
data.example2$strfact=0
data.example2$strfact[which(data.example2$strand=="+")]=1
names(data.example2)[c(21,23,24)] = c("Ma","Mg","Md")
data.example2$maexpr=data.example2$Ma*data.example2$express
data.example2$mgexpr=data.example2$Mg*data.example2$express
data.example2$mdexpr=data.example2$Md*data.example2$express



fparam = colnames(data.example2)[c(8:10,12:17,21,23,26,30,31,32,33)]
fobserved = colnames(data.example2)[5:6]


ns = 2
M=30
params = NULL
for(i in 1:M)
{
  withRestarts(tryCatch({
    params[[i]]=list(fparam = fparam,isobsbinary = c(0,0,rep(1,length(fparam))),prior.inclusion = array(1,c(length(fparam),2)),fobserved = fobserved,ranges = 1,ns = ns,initpr = c(0,1),data = data.example2,a.prior.em = runif(n = 1,min = 100,max = 200), b.prior.em = runif(n = 1,min = 1,max = 10), a.post.em = 10,wcrit = 1,a.Q.prior = 5, b.Q.prior = 15, X.min = 3, X.max = 10,a.maxit.prior  = runif(n = 1,min = 15, max = 35) ,b.maxit.prior  = 2 ,s2.post.maxit = 3 ,col.rate.maxit = runif(n=1,min = 0.97, max = 0.99),tmin = runif(n = 1,min = 0.0000005, max = 0.05),temp = runif(n = 1,min = 20000,max = 2000000000),dt = runif(n = 1,min = 2, max = 6) ,p_id = 1,seed = runif(1,1,1000))
  }),abort = function(){onerr=TRUE})
}
results = mclapply(FUN = function(x) do.call(simulated_annealing_binomial,x),X = params)

best = 2
for(i in 1:M)
{
  if(length(results[[i]])>1)
    if(BIC(results[[i]]$model)<BIC(results[[best]]$model))
      best = i 
}

results[[best]]$model

posts = depmixS4::posterior(results[[best]]$model)

png(file=paste("classifyhmm.png",sep = ""),width     = 10,
    height    = 5,
    units     = "in",
    res       = 500)
plot(data.example2$total_bases, col = 4, pch = 16,ylim = c(-13,30),xlab = "",ylab = "Data and methylation probabilities",xaxt = "n",yaxt = "n")
axis(1, at=seq(1,length(data.example2$pos), by = 50), labels=data.example2$pos[seq(1,length(data.example2$pos), by = 50)],las=2)
axis(2, at=c(seq(0,max(data.example2$total_bases)+5,5)), labels = c(seq(0,max(data.example2$total_bases)+5,5)))
axis(2, at=c(-13,-3), labels = c(0,1), las = 2)
points(data.example2$methylated_bases,col=2, pch = 15)
#lines(rep(2.6,length(data.example2$total_bases)),col=3,lwd=4,lty = 1)
lines(rep(- 2.1,length(data.example2$total_bases)))
#points(testtot,col=3, pch = 17)
#points(testmeth,col=7, pch = 17)
lines(10*data.example2$methylated_bases/data.example2$total_bases-13,col = 5,lwd=2)
#lines(10*res1[1,] - 11,col=8)
lines(10*(posts$S2)-13,lwd=2,col = "green")
lines(10*(posts$state-1)-13,lwd=2,col = "sienna4")
#lines(10*g(fm5$summary.random[[2]]$mode)-13,lwd=4,col = 1)
#lines(10*g(fm4$summary.random[[1]]$mean)-13,col =6,lwd=2)
dev.off()


#system("shutdown -h now")
simulated_annealing_binomial=function(epochs = 3, estat = 3, fparam,isobsbinary, MIC = stats::BIC,SIC =stats::AIC, fobserved, ranges, ns,initpr, data,a.prior.em, b.prior.em, a.post.em, wcrit, a.Q.prior, b.Q.prior, X.min, X.max,
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
  formula = as.formula(paste(paste("cbind(",fobserved[1],",",fobserved[2],")"), " ~ ",obsconst,ifelse(test = length(covobs)>0, " + ",""), paste(covobs, collapse=" + ")))
  transition =as.formula(paste(" ~ ",transconst,ifelse(test = length(covtrans)>0, " + ",""), paste(covtrans, collapse=" + ")))
  if(length(covtrans)>0 && length(covobs)>0)
  {
    
    formula = as.formula(paste(paste("cbind(",fobserved[1],",",fobserved[2],")"), " ~ ","1 + ", paste(covobs, collapse=" + ")))

    transition =as.formula(paste(" ~ ","1 + ", paste(covtrans, collapse=" + ")))
    
    mod = depmix(response = formula,data=data,transition = transition,nstates=ns, instart = initpr,
                 family=binomial())
  }else 
  {
    if(length(covtrans)>0 && length(covobs)==0)
    {
      transition =as.formula(paste(" ~ ","1 + ", paste(covtrans, collapse=" + ")))
      mod = depmix(response = as.formula(paste("cbind(",fobserved[1],",",fobserved[2],")", " ~ ","1")) ,data=data,transition = transition,nstates=ns, instart = initpr,family=binomial())
    }
    
    if(length(covtrans)==0 && length(covobs)>0)
    {
      
      formula = as.formula(paste(paste("cbind(",fobserved[1],",",fobserved[2],")"), " ~ ","1 + ", paste(covobs, collapse=" + ")))
      mod = depmix(response = formula,data=data,nstates=ns, instart = initpr,
                   family=binomial())
    }
    if(length(covtrans)==0 && length(covobs)==0)
    {
      mod = depmix(response = as.formula(paste("cbind(",fobserved[1],",",fobserved[2],")", " ~ ","1")),data=data,nstates=ns, instart = initpr, family=binomial())
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
          formula = as.formula(paste(paste("cbind(",fobserved[1],",",fobserved[2],")"), " ~ ",obsconst,ifelse(test = length(covobs)>0, " + ",""), paste(covobs, collapse=" + ")))
          transition =as.formula(paste(" ~ ",transconst,ifelse(test = length(covtrans)>0, " + ",""), paste(covtrans, collapse=" + ")))
          if(length(covtrans)>0 && length(covobs)>0)
          {  
          
            formula = as.formula(paste(paste("cbind(",fobserved[1],",",fobserved[2],")"), " ~ ","1 + ", paste(covobs, collapse=" + ")))
            transition =as.formula(paste(" ~ ","1 + ", paste(covtrans, collapse=" + ")))
            
            mod = depmix(response = formula,data=data,transition = transition,nstates=ns, instart = initpr,
                         family=binomial())
          }else 
          {
            if(length(covtrans)>0 && length(covobs)==0)
            {
             
              transition =as.formula(paste(" ~ ","1 + ", paste(covtrans, collapse=" + ")))
              mod = depmix(response = as.formula(paste("cbind(",fobserved[1],",",fobserved[2],")", " ~ ","1")),data=data,transition = transition,nstates=ns, instart = initpr,family=binomial())
            }
            
            if(length(covtrans)==0 && length(covobs)>0)
            {
              
              formula = as.formula(paste(paste("cbind(",fobserved[1],",",fobserved[2],")"), " ~ ","1 + ", paste(covobs, collapse=" + ")))
              
              mod = depmix(response = formula,data=data,nstates=ns, instart = initpr,
                           family=binomial())
            }
            if(length(covtrans)==0 && length(covobs)==0)
            {
              mod = depmix(response = as.formula(paste("cbind(",fobserved[1],",",fobserved[2],")", " ~ ","1")),data=data,nstates=ns, instart = initpr,family=binomial())
              
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





print(formula)

View(vars)

dim(X)
View(X[X$chrom == 3,][2:100000,1:9])
matplot(X[1000:7000,5:7],type="l")





ns=2

fparam = colnames(X)[9:18]
fobserved = colnames(X)[6:7]

simulated_annealing(fparam = fparam,fobserved = fobserved,ranges = 1,ns = ns,initpr = c(0,1),data = X[X$chrom == 1,][2:1000,1:18],em.params = em.control(maxit = 50,random.start = TRUE),wcrit = 1,maxit = 100,tmin = 0.5,temp = 1000,dt = 1.5,p_id = 100,seed = 5)

varcur=rbinom(n = length(fparam)+2,size = 2,prob = 0.5)
varcur[2]=as.integer(varcur[2]!=varcur[1])*varcur[2]
covobs = fparam[which(varcur[-c(1,2)] == 1)]
covtrans = fparam[which(varcur[-c(1,2)] == 2)]
obsconst=as.integer((varcur[1]==1||varcur[2]==1))
transconst=as.integer((varcur[1]==2||varcur[2]==2))
formula = as.formula(paste(paste("cbind(",fobserved[1],",",fobserved[2],")"), " ~ ",obsconst," + ", paste(covobs, collapse=" + ")))
transition =as.formula(paste(" ~ ",transconst,"+ ", paste(covtrans, collapse=" + ")))
mod = depmix(response = formula,data=X[X$chrom == 1,][2:1000,1:18],transition = transition,nstates=ns, instart = c(0,1),
             family=binomial())

set.seed(1)


fm= fit(mod)
sum=as.data.frame(summary(fm))
sum=summary(fm)[]
sum


AIC(fm)
BIC(fm)


globobs=array(data = 0,dim = c(ns,length(fparam)+2))
globtrans=array(data = 0,dim = c(ns,length(fparam)+2))

params=getpars(fm)
params
params[1]=0.11

tmp1=params[(ns + ns*(length(covtrans)+transconst)*2 + 1):(ns + ns*(length(covtrans)+transconst)*2 + ns*(length(covobs)+obsconst))]
globobs[1,which(varcur == 1)] =tmp1[1:(length(tmp1)/2)]
globobs[2,which(varcur == 1)] =tmp1[(length(tmp1)/2+1):length(tmp1)]
tmp1=params[ns+1:(ns*(length(covtrans)+transconst)*2)]
tmp1
globtrans[1,which(varcur == 2)] =tmp1[seq(2, ((length(tmp1)/2)), by=2)]
globtrans[2,which(varcur == 2)] =tmp1[seq((length(tmp1)/2+2),((length(tmp1))), by=2)]

globobs[1,which(varcur == 1)] 
# summary of the model
summary(fm)

globtrans[1,5]

print(fm)



getpars(fm)

mod=setpars(object = mod,values = getpars(fm))

params=getpars(fm)

# for observed process
params=getpars(mod)
tmp1=params[(ns + ns*(length(covtrans)+transconst)*2 + 1):(ns + ns*(length(covtrans)+transconst)*2 + ns*(length(covobs)+obsconst))]
tmp1[1:(length(tmp1)/2)]=curobs[1,which(varcur == 1)]
tmp1[(length(tmp1)/2+1):length(tmp1)]=curobs[2,which(varcur == 1)]
params[(ns + ns*(length(covtrans)+transconst)*2 + 1):(ns + ns*(length(covtrans)+transconst)*2 + ns*(length(covobs)+obsconst))]=tmp1
# for transitions of the hidden process
tmp1=params[ns+1:(ns*(length(covtrans)+transconst)*2)]
tmp1[seq(2, ((length(tmp1)/2)), by=2)]=curtrans[1,which(varcur == 2)] 
tmp1[seq((length(tmp1)/2+2),((length(tmp1))), by=2)]=curtrans[2,which(varcur == 2)]
params[ns+1:(ns*(length(covtrans)+transconst)*2)]=tmp1
# for initial probabilities of the hidden process
params[1:ns]=curinit
# write them down      
mod=setpars(object = mod,values = params)

print(params)

fit(mod)

AIC(mod)

curpar

logLik(object = mod,method = "classification") # FB algorithm or LH algorithm
logLik(object = mod,method = "fb")
logLik(object = mod,method = "lystig")

getpars(mod)


# update the parameters and likelihood:
updLogLik = function(p, model)
{
  lf = setpars(object = model,values = p)
}

res=nlm()

mod = setpars(object = mod,values = c(0.5,0.5,0.3,0.7,1,0,0.6,2.1,0.5,0.99))

getpars(mod)

logLik(object = mod,method = "classification") # FB algorithm or LH algorithm or given the Viterbi sequence
logLik(object = mod,method = "fb")
logLik(object = mod,method = "lystig")

post = posterior(fm)
dim(post)
View(post)

print("took 20 minutes")


FREQ = array(0,dim = c(2,2))
dim(X)[1]

for(i in 1:100000)
{
  
  if(X$methylation_call[i]==(post$state[i]-1) && X$methylation_call[i]==0)
  {
    FREQ[1,1]=FREQ[1,1]+1
  }
  if(X$methylation_call[i]==(post$state[i]-1) && X$methylation_call[i]==1)
  {
    FREQ[2,2]=FREQ[2,2]+1
  }
  if(X$methylation_call[i]!=(post$state[i]-1) && X$methylation_call[i]==0)
  {
    FREQ[2,1]=FREQ[2,1]+1
  }
  if(X$methylation_call[i]!=(post$state[i]-1) && X$methylation_call[i]==1)
  {
    FREQ[1,2]=FREQ[1,2]+1
  }
}

FREQ

