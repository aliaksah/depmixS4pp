**An Adaptive Simulated Annealing EM Algorithm for Inference on Non-homogeneous Hidden Markov Models**

<p align="justify">
In this package we are adding adaptive simulated variable selection, and various prediction options for the standard depmixS4 (https://cran.r-project.org/web/packages/depmixS4/) for inference on non-homogeneous hidden Markov models (NHHMM).
NHHMM are a subclass of dependent mixture models used for semi-supervised learning, where both transition probabilities between the latent states and the mean parameter of the probability distribution of the responses (for a given state) depend on the set of up to p covariates. A priori we do not know which (and how) covariates influence the transition probabilities and the mean parameters. This induces a complex combinatorial optimization problem for model selection with 4 to the power of p potential configurations. To address the problem, in this article we propose an adaptive (A) simulated annealing (SA) expectation maximization (EM) algorithm (ASA-EM) for joint optimization of models and their parameters with respect to a criterion of interest.
</p>

***

***

* Full text of the paper introducing An Adaptive Simulated Annealing EM Algorithm for Inference on Non-homogeneous Hidden Markov Models: [ACM](https://dl.acm.org/icps.cfm)
* Some applications of the package are available on: [GitHub](https://github.com/aliaksah/depmixS4pp/tree/master/examples)  


***

* Install binary on Linux or Mac Os:
```R 
install.packages("https://github.com/aliaksah/depmixS4pp/blob/master/depmixS4pp_1.0.tar.gz?raw=true", repos = NULL, type="source")
```
* An expert parallel call of parallel ASA-EM (see [select_depmix](https://rdrr.io/github/aliaksah/depmixS4pp/man/select_depmix.html)): 
```R 
results = depmixS4pp::select_depmix(epochs =3,estat = 3,data = X,MIC = stats::AIC,SIC =stats::BIC,family = gaussian(),fparam = fparam,fobserved = fobserved,isobsbinary = c(0,0,rep(1,length(fparam))),prior.inclusion = array(1,c(length(fparam),2)),ranges = 1,ns = ns,initpr =  c(0,1,0),seeds = runif(M,1,1000),cores = M)
```
* A simple call of prediction function (see [predict_depmix](https://rdrr.io/github/aliaksah/depmixS4pp/man/predict_depmix.html)): 
```R 
y.pred = depmixS4pp::predict_depmix(X.test,object = results$results[[results$best.mic]]$model,mode = F)
```
***

**Consider an example of Apple (AAPL) stock price predictions**

<p align="justify">
The example is focused upon Aple stock (AAPL) prediction with respect to log-returns of p=29 stocks from the S&P500 listing with the highest correlations to AAPL. The addressed data  consists of 1258 observations of log returns for 30 stocks based on the daily close price. Th predictors are: "ADI", "AMAT", "AMP", "AOS", "APH", "AVGO", "BLK", "BRK.B", "CSCO", "FISV", "GLW", "GOOGL", "GOOG", "HON"   "INTC", "ITW", "IVZ", "LRCX", "MA", "MCHP", "MSFT", "QRVO", "ROK", "SPGI", "SWKS", "TEL", "TSS", "TXN", "TXT", where the underscripts are representing the official S&P500 acronyms of the stocks. In terms of physical time, the observations are ranged from 11.02.2013 to 07.02.2018. The focus of this example is in both inference and predictions, hence we divided the data into a training data set (before 01.01.2017) and a testing data set (after 01.01.2017). The full data processing script is available at https://github.com/aliaksah/depmixS4pp/blob/master/examples/AAPL_example.R, but here let:
</p>

```R 
X #be the training data on 30 stocks
X.test #be the test data on 30 stocks
```

The primary model selection criterion addressed in this example is AIC due to that we are mainly interested in predictions, whilst BIC is the secondary reported criterion.

Then we specify initial formulas with all possible covariates for but the observed and latent processes as:

```R 
fparam = c("ADI", "AMAT", "AMP", "AOS", "APH", "AVGO", "BLK", "BRK.B", "CSCO", "FISV", "GLW", "GOOGL", "GOOG", "HON"   "INTC", "ITW", "IVZ", "LRCX", "MA", "MCHP", "MSFT", "QRVO", "ROK", "SPGI", "SWKS", "TEL", "TSS", "TXN", "TXT")
```
and the observations as


```R 
fobserved = "AAPL"
```
We will run the ASA-EM with gaussian observations with 3 latent states on 30 cores. Specify:

```R 
ns = 3
M=30
```
And run the inference with the stated in the call choice of tuning parameters of the algorithm:

```R
results = depmixS4pp::select_depmix(epochs =3,estat = 3,data = X,MIC = stats::AIC,SIC =stats::BIC,family = gaussian(),fparam = fparam,fobserved = fobserved,isobsbinary = c(0,0,rep(1,length(fparam))),prior.inclusion = array(1,c(length(fparam),2)),ranges = 1,ns = ns,initpr =  c(0,1,0),seeds = runif(M,1,1000),cores = M)
```

Now we can make the predictions (conditional on the mode state of NHHMM) as:

```R 
y.pred = depmixS4pp::predict_depmix(X.test,object = results$results[[results$best.mic]]$model,mode = T)
```

And plot the predictions for the log-returns (red) and actual test log-returns (green):


```R 
plot(X.test$AAPL, col =3,xaxt = "n")
points(y.pred$y, col = 2)
axis(1, at=seq(1,length(X.test$AAPL), by = 50), labels=X.test$date[seq(1,length(X.test$date), by = 50)],las=2)
```
<p align="center">
<img src="https://raw.githubusercontent.com/aliaksah/depmixS4pp/master/examples/results/logreturns.png" width="600" height="350">
</p>
Finally we can plot the price predictions. To do that we first trasform the actual and predicted log-returns via:

```R 
prc = sp*exp(cumsum(X.test$AAPL))
prc.pred = sp*exp(cumsum(y.pred$y))
``` 

And plot the predictions for the price (blue) and actual test price (brown) via:

```R 
plot(prc,col = "brown",type = "l",ylim = c(100,200),xaxt = "n")
lines(prc.pred, col = "blue")
axis(1, at=seq(1,length(X.test$AAPL), by = 50), labels=X.test$date[seq(1,length(X.test$date), by = 50)],las=2)
```
<p align="center">
<img src="https://raw.githubusercontent.com/aliaksah/depmixS4pp/master/examples/results/prices.png" width="600" height="350">
</p>
Here, we have the RMSE of **10.816** obtained via:

```R 
mse.hmm.pr = sqrt(mean((prc.pred-prc)^2))
```

Benchmarks against other methods can be found [here](https://github.com/aliaksah/depmixS4pp/blob/master/examples/AAPL_example.R).  

**Additionally the research was presented via the following selected contributions:**


*Scientific lectures*

1. Hubin, Aliaksandr; Storvik, Geir Olve. On model selection in Hidden Markov Models with covariates. Workshop, Klækken Workshop 2015; Klækken, 29.05.2015.

***

## Developed by [Aliaksandr Hubin](https://scholar.google.com/citations?user=Lx-G8ckAAAAJ&hl=en/)
 
 ***
