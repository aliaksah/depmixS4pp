\name{select_depmix}

\alias{select}


\title{select 'depmix' or 'mix' models by running simulated_annealing on multiple cores}

\description{

	\code{select_depmix} select based on simulated annealing \code{\link{depmix}} or
	\code{\link{mix}} models. 

}

\usage{

	select_depmix(object, test, mode, infer)

}

\arguments{

	\item{object}{An object of class \code{(dep-)mix}.}
	
	\item{object}{Parameters of NHHMM (if object is not present) as a list() with components: ns (number of states), fla (formula for observed process), trs (formula for process of transitions), family (family of the distributions of the responses), params (parameters of the fitted NHHMM).}
	
	\item{test}{Data set for predictions.}

	\item{mode}{Indicates whether the mode state will be used for predictions (most probable states' sequence will be used otherwise)}

}

\details{ 

	Prediction function.
	
}

\value{
	
	\code{predict_depmix} returns a vector of values of linear predictor (for glm you should further transform them through the appropriate link function)
	and the corresponding states of the latent process:
	
	
	\describe{
		\item{\code{y}:}{Values of the linear predictor.}
	
		\item{\code{s }:}{States corresponding to the values of the linear predictor.}
	
		
	}
	
	

}


\examples{

}

\author{Aliaksandr Hubin}

\keyword{methods}
