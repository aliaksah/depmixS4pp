


predict_depmix =	function(object,test, mode = FALSE,infer = NULL, which="pars",...) {
    
  if(length(object)==0&length(infer)==0)
  {
    print("Error: no object or parameters specified!")
    return(NULL)
  }
  if(length(object)>0&length(infer)>0)
  {
    print("Warning: both object and parameters are specified! Using object!")
    infer = NULL
  }
  
    if(length(object)==0&length(infer)>0)
    {
      object = depmix(response = infer$fla,data=test,transition =infer$trs, nstates=infer$ns,
                    family=infer$family)
    
      object = setpars(object,infer$params)
    } 
  
	  fb = forwardbackward(object)[["alpha"]]
	  predStates  = unlist(lapply(X=array(1:(dim(test)[1])),FUN = function(x){which.max(fb[x,])}))
	  if(length(predStates)==0)
	    predStates = rep(1,(dim(test)[1]))
	  if(mode)
	    predStates = rep(as.numeric(names(which.max(table(predStates)))),(dim(test)[1]))
	  dres = summary(object)
	  dres[is.na(dres)] = 0
	  preds = unlist(lapply(1:length(predStates), function(i){sum(unlist(c(1,test[i,-1]))*dres[predStates[i],1:dim(dres)[2]-1])}))
	  y.hat = preds
	  return(list(y= y.hat,s = predStates))
	}


