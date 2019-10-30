\name{simulated_annealing}

\alias{select_depmix}


\title{Predict 'depmix' or 'mix' models}

\description{

	\code{predict_depmix} predict based on \code{\link{depmix}} or
	\code{\link{mix}} models, optionally based on a mode state. 

}

\usage{

	\S4method{fit}{mix}(object, test, mode)

}

\arguments{

	\item{object}{An object of class \code{(dep-)mix}.}
	
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
