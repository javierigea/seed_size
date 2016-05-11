library(BAMMtools)
library(coda)
library(ape)
library(mvtnorm)
library(lattice)
#mcmcmout:mcmcoutfile
analyse_BAMM_convergence<-function(mcmcout,burnin){
  mcmcout.div <- read.csv(mcmcout, header=T)
  #see if convergence is ok
  plot(mcmcout.div$logLik ~ mcmcout.div$generation)
  #remove initial 10% (burnin)
  burnstart.div <- floor(burnin * nrow(mcmcout.div))
  postburn.div <- mcmcout.div[burnstart.div:nrow(mcmcout.div), ]
  plot(postburn.div$logLik ~ postburn.div$generation)
  plot(postburn.div$logLik ~ postburn.div$N_shifts)
  #check ESS (should be >200)
  print(paste('ESS Lik',effectiveSize(postburn.div$logLik),sep=': '))
  print(paste('ESS Nshifts',effectiveSize(postburn.div$N_shifts),sep=': '))
  #compute posterior probabilities of models
  post_probs.div <- table(postburn.div$N_shifts) / nrow(postburn.div)
  names(post_probs.div)
  #plot the posterior distribution of Nrateshifts
  plot(names(post_probs.div),post_probs.div)
}