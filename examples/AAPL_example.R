#load the libraries
library(depmixS4pp)
library(parallel)
library(forecast)
library(glmnet)
library(HDeconometrics)
library(RCurl)
library(depmixS4)
#prepare the data
sp500data = read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/depmixS4pp/master/data/all_stocks_5yr.csv"),sep = ",",header = T)
name1 = "AAPL"
tmp = sp500data[which(sp500data$Name==name1),c(1,5)]
spp = tmp$close[1]
sp = tmp$close[983]
tmp[is.na(tmp)] = 0
tmp[,2] = c(0,diff(log(tmp[,2]), lag=1))
names(tmp)[2] = name1
X = tmp
for(name in unique(sp500data$Name))
{
  if(name == name1)
    next
  tmp = sp500data[which(sp500data$Name==name),c(1,5)]
  tmp[is.na(tmp)] = 0
  tmp[,2]= c(0,diff(log(tmp[,2]), lag=1))
  names(tmp)[2] = name
  X = merge(X,tmp,by = "date",all.x = T) 
}
X = X[-1,]
X[is.na(X)] = 0
cors = cor(X[,-1])
slectid = c(1,which(abs(cors[1,])>0.333)+1)
X = X[,slectid]
X.test = X[983:1258,]
X = X[1:982,]

rm(cors)
rm(tmp)
rm(sp500data)
gc()


#run the parallel ASA-EM
fparam = colnames(X)[3:length(colnames(X))]
fobserved = colnames(X)[2]
ns = 3
M=30
results = select_depmix(epochs =3,estat = 3,data = X,MIC = stats::AIC,SIC =stats::BIC,family = gaussian(),fparam = fparam,fobserved = fobserved,isobsbinary = c(0,0,rep(1,length(fparam))),prior.inclusion = array(1,c(length(fparam),2)),ranges = 1,ns = ns,initpr =  c(0,1,0),seeds = runif(M,1,1000),cores = M)


#look at the best found model
results$results[[results$best.mic]]$model
write.table(depmixS4pp::getpars(results$results[[results$best.mic]]$model),"bestmodel.AAPL")
summary(depmixS4pp::getmodel(results$results[[results$best.mic]]$model))
posts = depmixS4::posterior(results$results[[results$best.mic]]$model)


#make inference
prc = spp*exp(cumsum(X$AAPL))


y.pred = depmixS4pp::predict_depmix(X.test,object = results$results[[results$best.mic]]$model,mode = T)


plot(y.pred$y, col = 2)
points(X.test$AAPL, col =3)

#fit and predict for the competing models
fit_autoarima = function(vect)
{
  tmp = cor(vect$x)
  tmp[upper.tri(tmp)] = 0
  diag(tmp) = 0
  xreg = vect$x[,!apply(tmp,2,function(x) any(abs(x) > 0.999))]
  fit = auto.arima(y = vect$y,xreg = xreg,max.p = 10, max.q = 10, max.d = 10,stepwise = T, seasonal = F, ic = "aicc")
  return(list(fit = fit, nxreg = colnames(xreg),vars = array(1,dim(vect$x)[2]))) 
}

predict_autoarima =function(infer,x.new)
{
  preds = forecast(infer$fit, y = x.new$y, h = 1, xreg = x.new$x[,which(colnames(x.new$x)%in%infer$nxreg)])$mean
  return(preds)
}
custfit.lasso = function(vect)
{
  lasso=ic.glmnet(x = vect$x,y=vect$y,family = "gaussian",alpha = 1)
  
  return(list(lasso = lasso,vars =  as.integer(lasso$glmnet$beta[,lasso$glmnet$dim[2]]!=0)))
}
custpredict.lasso = function(infer, x.new)
{
  return(predict(infer$lasso$glmnet,newx = x.new,type = "response")[,which(infer$lasso$glmnet$lambda == infer$lasso$lambda)])#
}
custfit.ridge = function(vect)
{
  ridge = ic.glmnet(x = vect$x,y=vect$y,family = "gaussian",alpha = 0)
  return(list(ridge = ridge, vars = as.integer(ridge$glmnet$beta[,ridge$glmnet$dim[2]]!=0)))
}
custpredict.ridge = function(infer, x.new)
{
  return(predict(infer$ridge$glmnet,newx = x.new,type = "response")[,which(infer$ridge$glmnet$lambda == infer$ridge$lambda)])
}

vect = NULL
vect$x = as.matrix(X[,-c(1,2)])
vect$y = X$AAPL
infer.ar = fit_autoarima(vect)
infer.lasso = custfit.lasso(vect)
infer.ridge = custfit.ridge(vect)


library(CausalImpact)
fit_causalimp = function(vect)
{
  train = as.matrix(vect$data.example)
  fm = CausalImpact(train , vect$trind, vect$teind, model.args = list(niter = 1000, nseasons = 1,standardize.data = F))
  coef = plot(fm$model$bsts.model, "coefficients")
  vars = array(0,length(coef$permutation))
  vars[coef$permutation]=(coef$inclusion.prob>0.5)
  return(list(vars = vars[-1],preds = as.vector(fm$series$point.pred[vect$teind[1]:vect$teind[2]])))
}
predict_causalimp = function(infer,x.new)
{
  return(infer$preds)
}

vect = NULL
vect$x = as.matrix(X.test[,-c(1,2)])
vect$y = X.test$AAPL
pred.ar = predict_autoarima(infer.ar,vect)
pred.lasso = custpredict.lasso(infer.lasso,as.matrix(X.test[,-c(1,2)]))
pred.ridge = custpredict.ridge(infer.ridge,as.matrix(X.test[,-c(1,2)]))
vect = list(data.example = rbind(X,X.test)[,-1], trind =c(1,982),teind = c(983,(983+dim(X.test)[1]-1)))
infer.cimpa = fit_causalimp(vect)
infer.cimpa$preds



prc = sp*exp(cumsum(X.test$AAPL))
prc.pred = sp*exp(cumsum(y.pred$y))
prc.cimpa = sp*exp(cumsum(infer.cimpa$preds))
prc.ar = sp*exp(cumsum(pred.ar))
prc.lasso = sp*exp(cumsum(pred.lasso))
prc.ridge = sp*exp(cumsum(pred.ridge))


library(prophet)
df = X[,c(1,2)]
names(df) = c('ds','y' )
m = prophet(df)
future = make_future_dataframe(m, periods = 276)
future$y = c(X$AAPL,X.test$AAPL)
tail(future)
forecast.prophet = predict(m, future)$yhat[1:276]
prc.proph =sp*exp(cumsum(forecast.prophet))


fit_depmix=function(vect)
{
  train = vect
  fla = as.formula(paste(colnames(train)[1],"~ 1 +",paste0(colnames(train)[-1],collapse = "+")))
  trs= as.formula(paste("~ 1 +",paste0(colnames(train)[-1],collapse = "+")))
  ns = 3
  mod = depmix(response = fla,data=train,transition =trs, nstates=ns,
               family=gaussian())
  fm=fit(mod,emcontrol = em.control(maxit = 200,random.start = TRUE))
  return(list(fm = fm, fla=fla, trs = trs,ns=ns))
}
hmm.vanilla = fit_depmix(X[,-1])
y.pred.vanilla = depmixS4pp::predict_depmix(test = X.test,object =hmm.vanilla$fm, mode = T)


prc.vanil = sp*exp(cumsum(y.pred.vanilla$y))


#print the results
mse.ar.pr = sqrt(mean((prc.ar-prc)^2))
mse.lasso.pr = sqrt(mean((prc.lasso-prc)^2))
mse.ridge.pr = sqrt(mean((prc.ridge-prc)^2))
mse.hmm.pr = sqrt(mean((prc.pred-prc)^2))
mse.cimpa.pr = sqrt(mean((prc.cimpa-prc)^2))
mse.prop.pr = sqrt(mean((prc.proph-prc)^2))
mse.vanilla.pr = sqrt(mean((prc.vanil-prc)^2))

print(paste(mse.hmm.pr,mse.lasso.pr,mse.ridge.pr,mse.ar.pr,mse.cimpa.pr,mse.prop.pr,mse.vanilla.pr))

forecast::dm.test((prc.ar-prc),(prc.pred-prc),h=1,alternative = "greater")
forecast::dm.test((prc.lasso-prc),(prc.pred-prc),h=1,alternative = "greater")
forecast::dm.test((prc.ridge-prc),(prc.pred-prc),h=1,alternative = "greater")
forecast::dm.test((prc.cimpa-prc),(prc.pred-prc),h=1,alternative = "greater")
forecast::dm.test((prc.proph-prc),(prc.pred-prc),h=1,alternative = "greater")
forecast::dm.test((prc.vanil-prc),(prc.pred-prc),h=1,alternative = "greater")


coint = function(phat,p)
{
  e = (p-phat)
  ars.p = arima(x = e,order = c(1,0,0),method = "ML")
  r.p = ars.p$coef[1]
  print(r.p)
  sdr.p = sqrt(ars.p$var.coef[1,1])
  max(abs(r.p - 1.96*sdr.p),abs(r.p + 1.96*sdr.p))
}

coint(prc,prc.pred)
coint(prc,prc.lasso)
coint(prc,prc.ridge)
coint(prc,prc.ar)
coint(prc,prc.cimpa)
coint(prc,prc.proph)
coint(prc,prc.vanil)


save.image("temp.rdata")
