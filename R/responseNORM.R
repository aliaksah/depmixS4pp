
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
	  if(is.matrix(y)) na <- unlist(apply(y,2,function(x) which(is.na(x)))) else na <- which(is.na(y))
	  if(is.matrix(x)) na <- c(na,unlist(apply(x,2,function(x) which(is.na(x))))) else na <- c(na,which(is.na(x)))
	  if(!is.null(w)) na <- c(na,which(is.na(w)))
	  #y <- as.matrix(y)
	  #x <- as.matrix(x)
	  na <- unique(na)
	  if(length(na)>0&length(na)<nrow(y)/1.2) {
	    x <- x[-na,]
	    y <- y[-na,]
	    #y <- round(y) # delete me
	    if(!is.null(w)) w <- w[-na]
	  }
	  else if(length(na)==nrow(y)&length(na)>0)
	  {
	    
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
	      xna = mean(w,na.rm = T)
	      if(is.nan(xna))
	        w[is.na(w)] <- 1
	      else
	        w[is.na(w)] <- xna
	    }
	  }
		pars <- object@parameters
		
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
