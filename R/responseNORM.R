
setClass("NORMresponse",contains="GLMresponse")

# method 'fit'
# use: in EM (M step)
# returns: (fitted) response with (new) estimates of parameters

setMethod("fit","NORMresponse",
	function(object,w) {
	  x=as.matrix(object@x)
	  y = as.matrix(object@y)
	  x[x == Inf] <- .Machine$double.max.exp
	  x[x == -Inf] <- .Machine$double.min.exp
	  y[y == Inf] <- .Machine$double.max.exp
	  y[y == -Inf] <- .Machine$double.min.exp
	  xna = mean(x,na.rm = T)
	  if(is.nan(xna))
	    x[is.na(x)] <- 0
	  else
	    x[is.na(x)] <- xna
	  xna = mean(y,na.rm = T)
	  if(is.nan(xna))
	    y[is.na(y)] <- 0
	  else
	    y[is.na(y)] <- xna
	  #y <- round(y) # delete me
	  if(!is.null(w))
	  {
	    w[w == Inf] <- .Machine$double.max.exp
	    w[w == -Inf] <- .Machine$double.min.exp
	    xna = mean(w,na.rm = T)
	    if(is.nan(xna))
	      w[is.na(w)] <- 1
	    else
	      w[is.na(w)] <- xna
	  }
		pars <- object@parameters
		
		if((sum(is.na(x))+sum(is.nan(x))+sum(is.infinite(x))+sum(is.na(y))+sum(is.nan(y))+sum(is.infinite(y)))==0&dim(x)[1]>0)
		{
  		if(!is.null(w)) {
  			fit <- lm.wfit(x=x,y=y,w=w)
  		} else {
  			fit <- lm.fit(x=x,y=y)
  		}
  		pars$coefficients <- fit$coefficients
  		if(!is.null(w)) {
  			pars$sd <- sqrt(sum(w*fit$residuals^2/sum(w)))
  		} else {
  			pars$sd <- sd(fit$residuals)
  		}
		}
		object <- setpars(object,unlist(pars))
		object
	}
)

setMethod("logDens","NORMresponse",
	function(object) {
		dnorm(x=object@y,mean=predict(object),sd=object@parameters$sd,log=TRUE)
	}
)

setMethod("dens","NORMresponse",
	function(object,log=FALSE) {
		dnorm(x=object@y,mean=predict(object),sd=object@parameters$sd,log=log)
	}
)

setMethod("predict","NORMresponse",
	function(object) {
		object@x%*%object@parameters$coefficients
	}
)

setMethod("simulate",signature(object="NORMresponse"),
  function(object,nsim=1,seed=NULL,times) {
    if(!is.null(seed)) set.seed(seed)
    if(missing(times)) {
      # draw in one go
      mu <- predict(object)
    } else {
      mu <- predict(object)[times]
    }  
    nt <- length(mu)
    sd <- object@parameters$sd
    response <- rnorm(nt*nsim,mean=mu,sd=sd)
    #if(nsim > 1) response <- matrix(response,ncol=nsim)
    response <- as.matrix(response)
    return(response)
  }
)
