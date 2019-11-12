*An Adaptive Simulated Annealing EM Algorithm for Inference on Non-homogeneous Hidden Markov Models*

<p align="justify">
In this package we are adding adaptive simulated variable selection, and various prediction options for the standard [depmixS4](https://dl.acm.org/icps.cfm) for inference on-homogeneous hidden Markov models (NHHMM)
NHHMM are a subclass of dependent mixture models used for semi-supervised learning, where both transition probabilities between the latent states and the mean parameter of the probability distribution of the responses (for a given state) depend on the set of up to p covariates. A priori we do not know which (and how) covariates influence the transition probabilities and the mean parameters. This induces a complex combinatorial optimization problem for model selection with 4^p potential configurations. To address the problem, in this article we propose an adaptive (A) simulated annealing (SA) expectation maximization (EM) algorithm (ASA-EM) for joint optimization of models and their parameters with respect to a criterion of interest.
</p>

***

***

* Full text of the paper introducing An Adaptive Simulated Annealing EM Algorithm for Inference on Non-homogeneous Hidden Markov Models: [ACM](https://dl.acm.org/icps.cfm)
* Some applications of the package are available on [GitHub](https://github.com/aliaksah/depmixS4pp/tree/master/examples)  


***

* Install binary on Linux or Mac Os:
```R 
install.packages("https://github.com/aliaksah/depmixS4pp/blob/master/depmixS4_1.0.tar.gz?raw=true", repos = NULL, type="source")
```
* An expert parallel call of ASA-EM (see [pinferunemjmcmc](https://rdrr.io/github/aliaksah/EMJMCMC2016/man/pinferunemjmcmc.html) for details): 
```R 
pinferunemjmcmc(n.cores =30, report.level =  0.8 , num.mod.best = NM,simplify = T, predict = T,test.data = as.data.frame(test),link.function = g, runemjmcmc.params =list(formula = formula1,data = data.example,gen.prob = c(1,1,1,1,0),estimator =estimate.bas.glm.cpen,estimator.args =  list(data = data.example,prior = aic.prior(),family = binomial(),yid=31, logn = log(143),r=exp(-0.5)),recalc_margin = 95, save.beta = T,interact = T,relations = c("gauss","tanh","atan","sin"),relations.prob =c(0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=4,mutation_rate = 100,last.mutation=1000, max.tree.size = 6, Nvars.max = 20,p.allow.replace=0.5,p.allow.tree=0.4,p.nor=0.3,p.and = 0.9),n.models = 7000,unique =T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(max.N.glob=as.integer(10), min.N.glob=as.integer(5), max.N=as.integer(3), min.N=as.integer(1), printable = F)))
```
* A simple call of parallel inference on Bayesian logic regression is (see [LogicRegr](https://rdrr.io/github/aliaksah/EMJMCMC2016/man/LogicRegr.html) for details): 
```R 
LogicRegr(formula = formula1,data = data.example,family = "Gaussian",prior = "G",report.level = 0.5,d = 15,cmax = 2,kmax = 15,p.and = 0.9,p.not = 0.01,p.surv = 0.2,ncores = 32)
```
***


**Additionally the research was presented via the following selected contributions:**


*Scientific lectures*

1. Hubin, Aliaksandr; Storvik, Geir Olve. On model selection in Hidden Markov Models with covariates. Workshop, Klækken Workshop 2015; Klækken, 29.05.2015.

***

## Developed by [Aliaksandr Hubin](https://scholar.google.com/citations?user=Lx-G8ckAAAAJ&hl=en/)
 
 ***
