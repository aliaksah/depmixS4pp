#load the libraries
library(depmixS4pp)
library(RCurl)
library(parallel)


#prepare the data
X = read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/Epigenetic%20Data/epigen.txt"),sep = ",",header = T)[,2:30]
X$strfact = 0
X$strfact[which(X$strand=="+")] = 1
names(X)[c(21,23,24)] = c("Ma","Mg","Md")
X$maexpr = X$Ma*X$express_noisy
X$mgexpr = X$Mg*X$express_noisy
X$mdexpr = X$Md*X$express_noisy
X$ppos = X$pos 
fparam = colnames(X)[c(8:10,12:17,21,23,26,29,30,31,32,33)]
fobserved = colnames(X)[5:6]

#run parallel ASA-EM
ns = 2
M=30
results = select_depmix(epochs = 4,estat = 3,data = X,MIC = stats::BIC,SIC =stats::AIC,family = binomial(),fparam = fparam,fobserved = fobserved,isobsbinary = c(0,0,rep(1,length(fparam))),prior.inclusion = array(1,c(length(fparam),2)),ranges = 1,ns = ns,initpr =  c(0,1),seed = runif(M,1:1000),cores = M)

#analyse the results
write.table(depmixS4pp::getpars(results$results[[results$best.mic]]$model),"bestmodel.epigenetic")
summary(results$results[[results$best.mic]]$model)
posts = depmixS4pp::posterior(results$results[[results$best.mic]]$model)

#plot methylations along the genome
plot(X$total_bases, col = 4, pch = 16,ylim = c(-13,30),xlab = "",ylab = "Data and methylation probabilities",xaxt = "n",yaxt = "n")
axis(1, at=seq(1,length(X$pos), by = 50), labels=X$pos[seq(1,length(X$pos), by = 50)],las=2)
axis(2, at=c(seq(0,max(X$total_bases)+5,5)), labels = c(seq(0,max(X$total_bases)+5,5)))
axis(2, at=c(-13,-3), labels = c(0,1), las = 2)
points(X$methylated_bases,col=2, pch = 15)
lines(rep(- 2.1,length(X$total_bases)))
lines(10*X$methylated_bases/X$total_bases-13,col = 5,lwd=2)
lines(10*(posts$S1)-13,lwd=3,col = "sienna4")
lines(10*(1-(posts$state-1))-13,lwd=1,col = "green")



