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

#this plots prior vs posterior with lines instead of barplots
#warning: xlim for the graph is now set as (0,300): change inside the function if different

plotPriormod<-function (mcmc, expectedNumberOfShifts = 1, burnin = 0.15, priorCol = "light blue", postCol = "red", legendPos = "topright", ...) 
{
  if (!class(mcmc) %in% c("character", "data.frame", "matrix")) {
    stop("mcmc must be either a dataframe or the path to the mcmc_out file.")
  }
  if (is.character(mcmc)) {
    mcmc <- read.csv(mcmc, stringsAsFactors = FALSE)
  }
  mcmc2 <- mcmc[floor(burnin * nrow(mcmc)):nrow(mcmc), ]
  obsK <- seq(from = 0, to = max(mcmc2[, "N_shifts"]), by = 1)
  prior <- sapply(obsK, prob.k, poissonRatePrior = 1/expectedNumberOfShifts)
  prior <- data.frame(N_shifts = obsK, prob = prior)
  posterior <- sapply(obsK, function(x) length(which(mcmc2[, 
                                                           "N_shifts"] == x)))/nrow(mcmc2)
  names(posterior) <- obsK
  posterior <- data.frame(N_shifts = names(posterior), prob = posterior)
  plot(c(0,0),xlim=c(0,300),ylim=c(0,0.2),xlab='n shifts',type='n', xaxs='i',yaxs='i',xaxt='n',ylab='',yaxt='n',border=NA)
  #barplot(prior[, 2], xlim=c(0,300),names.arg = prior[, 1], ylim = c(0,0.2),type='n', border = NA, col = priorCol, 
  #        xlab = "n shifts", xaxs='i',yaxs='i',xaxt='n',yaxt='n', ...)
  lines(prior[, 2], xlim=c(0,300),names.arg = prior[, 1], ylim = c(0,0.2), border = NA, col = priorCol, 
        xlab = "n shifts", xaxs='i',yaxs='i',xaxt='n',yaxt='n', lwd=3, ...)
  axis(1, at=seq(0,300,by=25),cex.axis=.7)
  axis(2, at=seq(0,0.2,by=0.05),las=2,cex.axis=.7)
  lines(posterior[, 2], add = TRUE, border = NA, col = BAMMtools::transparentColor(postCol, 
                                                                                   0.4), lwd=3,axes = FALSE)
  
  
  legend(x = legendPos, y = NULL, legend = c("prior", "posterior"), 
         fill = c(priorCol, BAMMtools::transparentColor(postCol, 
                                                        0.4)), bty = "n", cex = .7)
  invisible(cbind(N_shifts = prior$N_shifts, priorProbs = prior$prob, 
                  postProbs = posterior$prob))
}
prob.k <- function(k, poissonRatePrior=1)	{
  Denom <- (poissonRatePrior + 1)^(k+1)
  Prob <- poissonRatePrior / Denom
  return(Prob)
}
