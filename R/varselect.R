library("depmixS4")
library(RCurl)

rm(list = ls(all = TRUE))


M<-5
size<-1

data.example2 <- read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/Epigenetic%20Data/epigen.txt"),sep = ",",header = T)[,2:30]
data.example2$strfact<-0
data.example2$strfact[which(data.example2$strand=="+")]<-1
names(data.example2)[c(21,23,24)] = c("Ma","Mg","Md")
data.example2$maexpr<-data.example2$Ma*data.example2$express
data.example2$mgexpr<-data.example2$Mg*data.example2$express
data.example2$mdexpr<-data.example2$Md*data.example2$express



fparam <- colnames(data.example2)[c(8:10,12:17,21,23,26,30,31,32,33)]
fobserved <- colnames(data.example2)[5:6]


ns = 2
M<-1
globobs<-array(data = 0,dim = c(ns,length(fparam)+2))
globtrans<-array(data = 0,dim = c(ns,length(fparam)+2))
opt<-array(data = "",dim = 4*M+1)
opt[1]<-as.character(paste(c("MOD_ID","C1","C2",fparam,"P0","P1","AIC","BIC","TIME"), collapse = ";"))
for(i in 1:M)
{
  withRestarts(tryCatch({
    result<-simulated_annealing(fparam = fparam,isobsbinary = c(0,0,rep(1,length(fparam))),fobserved = fobserved,ranges = 1,ns = ns,initpr = c(0,1),data = data.example2,a.prior.em = runif(n = 1,min = 40,max = 70), b.prior.em = runif(n = 1,min = 9,max = 25), a.post.em = 10,wcrit = 1,a.Q.prior = 5, b.Q.prior = 15, X.min = 3, X.max = 10,mu.prior.maxit = runif(n = 1,min = 15, max = 35) ,s2.prior.maxit = 2 ,s2.post.maxit = 3 ,col.rate.maxit = runif(n=1,min = 0.97, max = 0.99),tmin = runif(n = 1,min = 0.0000005, max = 0.05),temp = runif(n = 1,min = 20000,max = 2000000000),dt = runif(n = 1,min = 2, max = 6) ,p_id = 1,seed = runif(1,1,1000))
    params<-getpars(result$model)  # get the current values
    #update global values
    covobs <- fparam[which(result$vars[-c(1,2)] == 1)]
    covtrans <- fparam[which(result$vars[-c(1,2)] == 2)]
    obsconst<-as.integer((result$vars[1]==1||result$vars[2]==1))
    transconst<-as.integer((result$vars[1]==2||result$vars[2]==2))
    globobs<-array(data = 0,dim = c(ns,length(fparam)+2))
    globtrans<-array(data = 0,dim = c(ns,length(fparam)+2))
    tmp1<-params[(ns + ns*(length(covtrans)+transconst)*2 + 1):(ns + ns*(length(covtrans)+transconst)*2 + ns*(length(covobs)+obsconst))]
    globobs[1,which(result$vars == 1)] <-tmp1[1:(length(tmp1)/2)]
    globobs[2,which(result$vars == 1)] <-tmp1[(length(tmp1)/2+1):length(tmp1)]
    tmp2<-params[ns+1:(ns*(length(covtrans)+transconst)*2)]
    if(length(tmp2)>2)
    {
      globtrans[1,which(result$vars== 2)] <-tmp2[seq(2, ((length(tmp2)/2)), by=2)]
      globtrans[2,which(result$vars == 2)] <-tmp2[seq((length(tmp2)/2+2),((length(tmp2))), by=2)]
    }
    
    opt[4*(i-1)+2]<-as.character(paste(c(paste(i,"obs"),globobs[1,],params[1:2],AIC(result$model),BIC(result$model),result$time), collapse = ";"))
    opt[4*(i-1)+3]<-as.character(paste(c(paste(i,"obs"),globobs[2,],params[1:2],AIC(result$model),BIC(result$model),result$time), collapse = ";"))
    opt[4*(i-1)+4]<-as.character(paste(c(paste(i,"trs"),globtrans[1,],params[1:2],AIC(result$model),BIC(result$model),result$time), collapse = ";"))
    opt[4*(i-1)+5]<-as.character(paste(c(paste(i,"trs"),globtrans[2,],params[1:2],AIC(result$model),BIC(result$model),result$time), collapse = ";"))
    print(i)                 
    
    
  }),abort = function(){onerr<-TRUE})
}

fileConn<-file("/mn/anatu/ansatte-u3/aliaksah/Desktop/Untitled Folder/output1.csv")
for(i in 2:M)
{
  write(opt[i],file=fileConn,append=TRUE)
}
#writeLines(opt, fileConn)
close(fileConn)

#system("shutdown -h now")


simulated_annealing(fparam = fparam,isobsbinary = c(0,0,1,1,1,1,1,1,1,1,1,1),fobserved = fobserved,ranges = 1,ns = ns,initpr = c(0,1),data = X[X$chrom == 1,][2:1000,1:18],a.prior.em = 50, b.prior.em = 15, a.post.em = 10,wcrit = 1,a.Q.prior = 5, b.Q.prior = 15, X.min = 3, X.max = 10,mu.prior.maxit = 25 ,s2.prior.maxit = 2 ,s2.post.maxit = 3 ,col.rate.maxit = 0.99,tmin = 0.0005,temp = 200000,dt = 3 ,p_id = 1,seed = runif(1,1,1000))
## simulated annealing optimization with limited EM greedy improvements
simulated_annealing<-function(fparam,isobsbinary, fobserved, ranges, ns,initpr, data,a.prior.em, b.prior.em, a.post.em, wcrit, a.Q.prior, b.Q.prior, X.min, X.max,
                              mu.prior.maxit, s2.prior.maxit, s2.post.maxit,col.rate.maxit,tmin,temp,dt,p_id,seed){
  # initialize the set of parameters
  maxit<-mu.prior.maxit
  stm <- proc.time()
  set.seed(seed, kind = NULL, normal.kind = NULL)
  print(paste("begin SA methasearch procedure","with EM local search greedy updates for each neighbourhood!"))
  tcur<-temp
  n.obs.em <- 0
  sum.obs.em <-0
  n.obs.maxit <- 0
  sum.obs.maxit <-0
  n.obs.k <- 0
  sum.obs.k <-0
  
  tid<-0
  ii<-1
  
  globcrit<-array(2,data = Inf)
  
  fm = NULL
  varcur<-rbinom(n = length(fparam)+2,size = 2,prob = 0.5)
  N<-length(varcur)
  varcur[2]=as.integer(varcur[2]!=varcur[1])*varcur[2]
  covobs <- fparam[which(varcur[-c(1,2)] == 1)]
  covtrans <- fparam[which(varcur[-c(1,2)] == 2)]
  obsconst<-as.integer((varcur[1]==1||varcur[2]==1))
  transconst<-as.integer((varcur[1]==2||varcur[2]==2))
  formula <- as.formula(paste(paste("cbind(",fobserved[1],",",fobserved[2],")"), " ~ ",obsconst,ifelse(test = length(covobs)>0, " + ",""), paste(covobs, collapse=" + ")))
  transition <-as.formula(paste(" ~ ",transconst,ifelse(test = length(covtrans)>0, " + ",""), paste(covtrans, collapse=" + ")))
  if(length(covtrans)>0 && length(covobs)>0)
  {
    if(sum(varcur[which(varcur == 1)])==sum(isobsbinary[which(varcur == 1)]))
      formula <- as.formula(paste(paste("cbind(",fobserved[1],",",fobserved[2],")"), " ~ ","1 + ", paste(covobs, collapse=" + ")))
    if(sum(varcur[which(varcur == 2)])==2*sum(isobsbinary[which(varcur == 1)]))
      transition <-as.formula(paste(" ~ ","1 + ", paste(covtrans, collapse=" + ")))
    
    mod <- depmix(response = formula,data=data,transition = transition,nstates=ns, instart = initpr,
                  family=binomial())
  }else 
  {
    if(length(covtrans)>0 && length(covobs)==0)
    {
      if(sum(varcur[which(varcur == 2)])==2*sum(isobsbinary[which(varcur == 1)]))
        transition <-as.formula(paste(" ~ ","1 + ", paste(covtrans, collapse=" + ")))
      mod <- depmix(response = as.formula(paste("cbind(",fobserved[1],",",fobserved[2],")", " ~ ","1")) ,data=data,transition = transition,nstates=ns, instart = initpr,family=binomial())
    }
    
    if(length(covtrans)==0 && length(covobs)>0)
    {
      if(sum(varcur[which(varcur == 1)])==sum(isobsbinary[which(varcur == 1)]))
        formula <- as.formula(paste(paste("cbind(",fobserved[1],",",fobserved[2],")"), " ~ ","1 + ", paste(covobs, collapse=" + ")))
      mod <- depmix(response = formula,data=data,nstates=ns, instart = initpr,
                    family=binomial())
    }
    if(length(covtrans)==0 && length(covobs)==0)
    {
      mod <- depmix(response = as.formula(paste("cbind(",fobserved[1],",",fobserved[2],")", " ~ ","1")),data=data,nstates=ns, instart = initpr, family=binomial())
    }
  }
  curn<-rgamma(n = 1,shape = (a.prior.em + n.obs.em*a.post.em),rate = (b.prior.em + sum.obs.em))
  em.params<-em.control(maxit = curn,random.start = TRUE)
  
  #fm<-fit(mod,emcontrol = em.params)
  
  capture.output({withRestarts(tryCatch(capture.output({fm<-fit(mod,emcontrol = em.params)}),error = function(e){invisible(e)}), 
                               abort = function(){ fm<-mod;  onerr<-TRUE})}) # fit the modal, get local improvements
  
  
  globcrit[1]<-AIC(fm)
  globcrit[2]<-BIC(fm)
  varglob<-varcur
  glmod<-fm
  
  curcrit<-globcrit
  
  
  
  while(tcur>=tmin)
  {
    if(tid %% p_id == 0){
      
      n.obs.maxit <- 0
      sum.obs.maxit <-0
      maxit <- abs(as.integer((mu.prior.maxit/(s2.prior.maxit) +sum.obs.maxit/s2.post.maxit)/(1/(s2.prior.maxit) + n.obs.maxit/s2.post.maxit)))
      sigma <-1/(1/(s2.prior.maxit) + n.obs.maxit/s2.post.maxit)
      print(paste(ii," iterations completed up to now"))
      print(paste("AIC.cur.glob = ", globcrit[1]))
      print(paste("aneal to temp = ",tcur," proceed with this temp. for",
                  maxit, "iterations"))
      
    }
    for(i in 1:maxit)
    {
      withRestarts(tryCatch({
        
        # change the neighbourhood i.e. randomly change two positions in the given vector
        varcurb<-varcur
        
        K<-rbinom(n = 1,size = N,prob = rbeta(n = 1,shape1 = (a.Q.prior + sum.obs.k),shape2 = (b.Q.prior + N*n.obs.k - sum.obs.k)))
        for(t in 1:K)
        {
          varcur[as.integer(runif(n = 1,min = 1,max = length(varcur)))] = rbinom(n = 1,size = 2,prob = 0.5)
        }     
        varcur[2]=as.integer(varcur[2]!=varcur[1])*varcur[2]
        covobs <- fparam[which(varcur[-c(1,2)] == 1)]
        covtrans <- fparam[which(varcur[-c(1,2)] == 2)]
        obsconst<-as.integer((varcur[1]==1||varcur[2]==1))
        transconst<-as.integer((varcur[1]==2||varcur[2]==2))
        formula <- as.formula(paste(paste("cbind(",fobserved[1],",",fobserved[2],")"), " ~ ",obsconst,ifelse(test = length(covobs)>0, " + ",""), paste(covobs, collapse=" + ")))
        transition <-as.formula(paste(" ~ ",transconst,ifelse(test = length(covtrans)>0, " + ",""), paste(covtrans, collapse=" + ")))
        if(length(covtrans)>0 && length(covobs)>0)
        {  
          if(sum(varcur[which(varcur == 1)])==sum(isobsbinary[which(varcur == 1)]))
            formula <- as.formula(paste(paste("cbind(",fobserved[1],",",fobserved[2],")"), " ~ ","1 + ", paste(covobs, collapse=" + ")))
          
          if(sum(varcur[which(varcur == 2)])==2*sum(isobsbinary[which(varcur == 1)]))
            transition <-as.formula(paste(" ~ ","1 + ", paste(covtrans, collapse=" + ")))
          
          mod <- depmix(response = formula,data=data,transition = transition,nstates=ns, instart = initpr,
                        family=binomial())
        }else 
        {
          if(length(covtrans)>0 && length(covobs)==0)
          {
            if(sum(varcur[which(varcur == 2)])==2*sum(isobsbinary[which(varcur == 1)]))
              transition <-as.formula(paste(" ~ ","1 + ", paste(covtrans, collapse=" + ")))
            mod <- depmix(response = as.formula(paste("cbind(",fobserved[1],",",fobserved[2],")", " ~ ","1")),data=data,transition = transition,nstates=ns, instart = initpr,family=binomial())
          }
          
          if(length(covtrans)==0 && length(covobs)>0)
          {
            if(sum(varcur[which(varcur == 1)])==sum(isobsbinary[which(varcur == 1)]))
              formula <- as.formula(paste(paste("cbind(",fobserved[1],",",fobserved[2],")"), " ~ ","1 + ", paste(covobs, collapse=" + ")))
            
            mod <- depmix(response = formula,data=data,nstates=ns, instart = initpr,
                          family=binomial())
          }
          if(length(covtrans)==0 && length(covobs)==0)
          {
            mod <- depmix(response = as.formula(paste("cbind(",fobserved[1],",",fobserved[2],")", " ~ ","1")),data=data,nstates=ns, instart = initpr,family=binomial())
            
          }
        }
        
        # set initial values for the parameters with respect to a given neighbourhood
        
        
        fmb<-fm
        onerr<-FALSE
        curn<-rgamma(n = 1,shape = (a.prior.em + n.obs.em*a.post.em),rate = (b.prior.em + sum.obs.em))
        em.params<-em.control(maxit = curn,random.start = TRUE)
        
        #finally = print(paste("loop body finished in iteration ",ii))
        # optimize in the new neighbourhood
        capture.output({withRestarts(tryCatch(capture.output({fm<-fit(mod,emcontrol = em.params)}),error = function(e){invisible(e)}), 
                                     abort = function(){varcur<-varcurb; fm<-fmb;  onerr<-TRUE})}) # fit the modal, get local improvements
        
        ii<-ii+1
        
        #update local values
        
        
        if(AIC(fm)<globcrit[1])
        {
          print(paste("update global optimum with new global AIC = ", AIC(fm)))
          
          globcrit[1]<-AIC(fm)
          globcrit[2]<-BIC(fm)
          glmod<-fm
          varglob<-varcur
          n.obs.em<-(n.obs.em + 1)
          sum.obs.em<-(sum.obs.em + curn)
          sum.obs.k<-n.obs.k*X.min
          
          
        }
        if(((AIC(fm)-curcrit[1]))/exp(tcur)<=runif(n = 1,min = 0,max = 1) && !onerr)
        {
          
          curcrit[1]<-AIC(fm)
          curcrit[2]<-BIC(fm)
          
          sum.obs.k<-sum.obs.k + N*a.Q.prior/(a.Q.prior + b.Q.prior)
          n.obs.k<-n.obs.k + 1
          
        }else
        {
          varcur<-varcurb
          n.obs.maxit <- (n.obs.maxit +1)
          sum.obs.maxit <-(sum.obs.maxit + rnorm(n = 1,mean = maxit,sd = sqrt(sigma)))*col.rate.maxit
          maxit <- abs(as.integer((mu.prior.maxit/(s2.prior.maxit) +sum.obs.maxit/s2.post.maxit)/(1/(s2.prior.maxit) + n.obs.maxit/s2.post.maxit)))
          sigma <-1/(1/(s2.prior.maxit) + n.obs.maxit/s2.post.maxit)
          
          K<-K+1
          if(K>X.max)
            K=X.max
          sum.obs.k<-sum.obs.k + K
          n.obs.k<-n.obs.k + 1
        }
        
      }),abort = function(){varcur<-varcurb; fm<-fmb;  onerr<-TRUE})
    }
    # update temperature
    tcur<-tcur*exp(-dt)
    tid<-tid+1
    n.obs.k <- 0
    sum.obs.k <-0
  }
  print(paste(ii," iterations completed in total"))
  print(paste(tcur," is the final temperature"))
  
  print(paste("AIC.cur.glob = ", globcrit[1]))
  # update the selected model to perfection
  withRestarts(tryCatch({
    capture.output({fm<-fit(glmod)}) # fit the modal, get local improvements
  }),abort = function(){capture.output({fm<-glmod});  onerr<-TRUE})
  etm <- proc.time()
  tt<-(etm[3]-stm[3])*1000
  
  return(list(model=fm, vars=varglob, time=(tt)))
}

#d <- read.table("../data/160215_GSM1085222_mC_calls_Col_0.tsv",nrows=100000,header=T)
#fit.bin <- glm(cbind(d$methylated_bases,d$total_bases)~mc_class,data=d)
#d$beta <- d$methylated_bases/d$total_bases
#boxplot(d$beta~d$mc_class)
#fit <- lm(beta~mc_class,data=d)
#library(INLA)
#fit.inla <- inla(methylation_call~mc_class+f(pos,model="ar1"),family="binomial",data=d[1:1000,])
#fit.inla <- inla(methylation_call~mc_class+f(pos,model="ar1",hyper=list(theta1=list(initial=0.01))),family="binomial",data=d[1:50000,])

#data <- read.table("http://folk.uio.no/egrini/CELS/GSM1085222_mC_calls_Col_0%202.txt")




print(formula)

View(vars)

dim(X)
View(X[X$chrom == 3,][2:100000,1:9])
matplot(X[1000:7000,5:7],type="l")





ns<-2

fparam <- colnames(X)[9:18]
fobserved <- colnames(X)[6:7]

simulated_annealing(fparam = fparam,fobserved = fobserved,ranges = 1,ns = ns,initpr = c(0,1),data = X[X$chrom == 1,][2:1000,1:18],em.params = em.control(maxit = 50,random.start = TRUE),wcrit = 1,maxit = 100,tmin = 0.5,temp = 1000,dt = 1.5,p_id = 100,seed = 5)

varcur<-rbinom(n = length(fparam)+2,size = 2,prob = 0.5)
varcur[2]=as.integer(varcur[2]!=varcur[1])*varcur[2]
covobs <- fparam[which(varcur[-c(1,2)] == 1)]
covtrans <- fparam[which(varcur[-c(1,2)] == 2)]
obsconst<-as.integer((varcur[1]==1||varcur[2]==1))
transconst<-as.integer((varcur[1]==2||varcur[2]==2))
formula <- as.formula(paste(paste("cbind(",fobserved[1],",",fobserved[2],")"), " ~ ",obsconst," + ", paste(covobs, collapse=" + ")))
transition <-as.formula(paste(" ~ ",transconst,"+ ", paste(covtrans, collapse=" + ")))
mod <- depmix(response = formula,data=X[X$chrom == 1,][2:1000,1:18],transition = transition,nstates=ns, instart = c(0,1),
              family=binomial())

set.seed(1)


fm<- fit(mod)
sum<-as.data.frame(summary(fm))
sum<-summary(fm)[]
sum


AIC(fm)
BIC(fm)


globobs<-array(data = 0,dim = c(ns,length(fparam)+2))
globtrans<-array(data = 0,dim = c(ns,length(fparam)+2))

params<-getpars(fm)
params
params[1]=0.11

tmp1<-params[(ns + ns*(length(covtrans)+transconst)*2 + 1):(ns + ns*(length(covtrans)+transconst)*2 + ns*(length(covobs)+obsconst))]
globobs[1,which(varcur == 1)] <-tmp1[1:(length(tmp1)/2)]
globobs[2,which(varcur == 1)] <-tmp1[(length(tmp1)/2+1):length(tmp1)]
tmp1<-params[ns+1:(ns*(length(covtrans)+transconst)*2)]
tmp1
globtrans[1,which(varcur == 2)] <-tmp1[seq(2, ((length(tmp1)/2)), by=2)]
globtrans[2,which(varcur == 2)] <-tmp1[seq((length(tmp1)/2+2),((length(tmp1))), by=2)]

globobs[1,which(varcur == 1)] 
# summary of the model
summary(fm)

globtrans[1,5]

print(fm)



getpars(fm)

mod<-setpars(object = mod,values = getpars(fm))

params<-getpars(fm)

# for observed process
params<-getpars(mod)
tmp1<-params[(ns + ns*(length(covtrans)+transconst)*2 + 1):(ns + ns*(length(covtrans)+transconst)*2 + ns*(length(covobs)+obsconst))]
tmp1[1:(length(tmp1)/2)]<-curobs[1,which(varcur == 1)]
tmp1[(length(tmp1)/2+1):length(tmp1)]<-curobs[2,which(varcur == 1)]
params[(ns + ns*(length(covtrans)+transconst)*2 + 1):(ns + ns*(length(covtrans)+transconst)*2 + ns*(length(covobs)+obsconst))]<-tmp1
# for transitions of the hidden process
tmp1<-params[ns+1:(ns*(length(covtrans)+transconst)*2)]
tmp1[seq(2, ((length(tmp1)/2)), by=2)]<-curtrans[1,which(varcur == 2)] 
tmp1[seq((length(tmp1)/2+2),((length(tmp1))), by=2)]<-curtrans[2,which(varcur == 2)]
params[ns+1:(ns*(length(covtrans)+transconst)*2)]<-tmp1
# for initial probabilities of the hidden process
params[1:ns]<-curinit
# write them down      
mod<-setpars(object = mod,values = params)

print(params)

fit(mod)

AIC(mod)

curpar

logLik(object = mod,method = "classification") # FB algorithm or LH algorithm
logLik(object = mod,method = "fb")
logLik(object = mod,method = "lystig")

getpars(mod)


# update the parameters and likelihood:
updLogLik <- function(p, model)
{
  lf <- setpars(object = model,values = p)
}

res<-nlm()

mod <- setpars(object = mod,values = c(0.5,0.5,0.3,0.7,1,0,0.6,2.1,0.5,0.99))

getpars(mod)

logLik(object = mod,method = "classification") # FB algorithm or LH algorithm or given the Viterbi sequence
logLik(object = mod,method = "fb")
logLik(object = mod,method = "lystig")

post <- posterior(fm)
dim(post)
View(post)

print("took 20 minutes")


FREQ <- array(0,dim = c(2,2))
dim(X)[1]

for(i in 1:100000)
{
  
  if(X$methylation_call[i]==(post$state[i]-1) && X$methylation_call[i]==0)
  {
    FREQ[1,1]<-FREQ[1,1]+1
  }
  if(X$methylation_call[i]==(post$state[i]-1) && X$methylation_call[i]==1)
  {
    FREQ[2,2]<-FREQ[2,2]+1
  }
  if(X$methylation_call[i]!=(post$state[i]-1) && X$methylation_call[i]==0)
  {
    FREQ[2,1]<-FREQ[2,1]+1
  }
  if(X$methylation_call[i]!=(post$state[i]-1) && X$methylation_call[i]==1)
  {
    FREQ[1,2]<-FREQ[1,2]+1
  }
}

FREQ

