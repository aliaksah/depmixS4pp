

setMethod("predict_depmix","depmix",
	function(object,test, which="pars",...) {
	  fb = forwardbackward(object)[["alpha"]]
	  predStates  = unlist(lapply(X=array(1:(dim(test)[1])),FUN = function(x){which.max(fb[x,])}))
	  if(length(predStates)==0)
	    predStates = rep(1,(dim(test)[1]))
	  dres = summary(object)
	  dres[is.na(dres)] = 0
	  preds = unlist(lapply(1:length(predStates), function(i){sum(unlist(c(1,test[i,-1]))*dres[predStates[i],1:dim(dres)[2]-1])}))
	  y.hat = preds
	  return(list(y= y.hat,s = predStates))
	}
)


setMethod("predict_depmix_mode","depmix",
          function(object,test, which="pars",...) {
            fb = forwardbackward(object)[["alpha"]]
            predStates  = unlist(lapply(X=array(1:(dim(test)[1])),FUN = function(x){which.max(fb[x,])}))
            if(length(predStates)==0)
              predStates = rep(1,(dim(test)[1]))
            predStates = rep(as.numeric(names(which.max(table(predStates)))),(dim(test)[1]))
            dres = summary(object)
            dres[is.na(dres)] = 0
            preds = unlist(lapply(1:length(predStates), function(i){sum(unlist(c(1,test[i,-1]))*dres[predStates[i],1:dim(dres)[2]-1])}))
            y.hat = preds
            return(list(y= y.hat,s = predStates))
          }
)

