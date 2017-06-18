#################################################################################
# (c) Javier Igea
# igea.javier@gmail.com
# Last modified: 18/06/2017

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#################################################################################
library(BAMMtools)
library(coda)
library(ape)
library(mvtnorm)
library(lattice)
library(plyr)
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

########this function plots multiple BAMM runs (paths in BAMMrunfiles)
plot_multiple_BAMM_runs<-function(BAMMrunfiles){
  palette<-rainbow(n=length(BAMMrunfiles))
  for (i in 1:length(BAMMrunfiles)){
    mcmcout<- read.csv(BAMMrunfiles[i], header=T)
    if(i==1){
      plot(mcmcout$logLik ~ mcmcout$generation,col=palette[i],pch=16)
    }
    else{
      points(mcmcout$logLik ~ mcmcout$generation,col=palette[i],pch=16)  
    }
    
  }
  legend('bottomright',legend=as.character(c(seq(1:length(BAMMrunfiles)))),col=palette,lty=1,bty='n',lwd=2)
}
########this function plots multiple BAMM runs (paths in BAMMrunfiles)
combine_multiple_BAMM_runs<-function(BAMMrunfiles,burnin.list,name){
  palette<-rainbow(n=length(BAMMrunfiles))
  list.mcmc<-list()
  for (i in 1:length(BAMMrunfiles)){
    mcmcout<- read.csv(BAMMrunfiles[i], header=T)
    burnin<-burnin.list[i]
    burnstart <- floor(burnin * nrow(mcmcout))
    postburn <- mcmcout[burnstart:nrow(mcmcout), ]
    list.mcmc[[i]]<-postburn
  }
  
  for (i in 2:length(list.mcmc)){
    end.previous<-list.mcmc[[i-1]][nrow(list.mcmc[[i-1]]),'generation']
    interval<-list.mcmc[[i]][2,'generation']-list.mcmc[[i]][1,'generation']
    start<-end.previous+interval
    list.mcmc[[i]]$generation<-seq(from=start,to=(nrow(list.mcmc[[i]])*interval+start)-interval,by=interval)
  }
  for (i in 1:length(list.mcmc)){
    if(i == 1){
      plot(list.mcmc[[i]]$logLik~list.mcmc[[i]]$generation,pch=16,xlim=c(0,list.mcmc[[length(list.mcmc)]][nrow(list.mcmc[[length(list.mcmc)]]),'generation']),col=palette[i])
    }
    else{
      points(list.mcmc[[i]]$logLik~list.mcmc[[i]]$generation,pch=16,col=palette[i])
    }
  }
  legend('bottomright',legend=as.character(c(seq(1:length(BAMMrunfiles)))),col=palette,lty=1,bty='n',lwd=2)
  mcmcout <- do.call("rbind", list.mcmc)
  plot(mcmcout$logLik ~ mcmcout$generation)
  plot(mcmcout$logLik ~ mcmcout$N_shifts)
  cat(paste('ESS Lik',effectiveSize(mcmcout$logLik),sep=': '))
  cat(paste('ESS Nshifts',effectiveSize(mcmcout$N_shifts),sep=': '))
  write.csv(mcmcout,file=paste(name,'_mcmcout_burninremoved.txt',sep=''))
  post_probs.div <- table(mcmcout$N_shifts) / nrow(mcmcout)
  names(post_probs.div)
  #plot the posterior distribution of Nrateshifts
  plot(names(post_probs.div),post_probs.div)
  
}

#this function combines a list of event data files (as a vector in BAMMeventfiles) with burnins (as a vector in burnin.list) and outputs a *name*_concatenate_event_data.txt
#this assumes that the samples in the event file have been taken every 10000 generation
concatenate_event_data<-function(BAMMeventfiles,burnin.list,name){
  #get ngenerations from files
  ngens<-sapply(BAMMeventfiles, function (x) {y<-system(sprintf('tail -n1 %s | cut -d "," -f1',x),intern=TRUE); y<-as.numeric(sub(y,pattern=x,replacement=''))})
  #apply burnin to each file
  ngens.new<-round_any(ngens*burnin.list,10000)
  csv.list<-list()
  for (i in 1:length(ngens)){
    csv.list[[i]]<-read.csv(BAMMeventfiles[i],header=TRUE)
    csv.list[[i]]<-csv.list[[i]][csv.list[[i]]$generation>=ngens.new[i],]
    csv.list[[i]]$generation<-csv.list[[i]]$generation-ngens.new[i]
    cat(paste(i,' read',sep=''),'\n')
  }
  #changing the number of generations in each file
  for (i in 1:length(csv.list)){
    if(i == 1){
      end<-csv.list[[i]][nrow(csv.list[[i]]),'generation']
    }
    else{
      csv.list[[i]]$generation<-csv.list[[i]]$generation+end+10000
      end<-csv.list[[i]][nrow(csv.list[[i]]),'generation']
    }
    cat(paste(i,' processed',sep=''),'\n')
  }
  #read all eventfiles into one
  concat.csv <- do.call("rbind", csv.list)
  write.csv(concat.csv,file=paste(name,'_concatenate_event_data.txt',sep=''),quote=F,row.names = FALSE)
  
}




#this plots prior vs posterior with lines instead of barplots
#warning: xlim for the graph is now set as (0,300): change inside the function if different

plotPriormod<-function (mcmc, expectedNumberOfShifts = 1, burnin = 0.15, priorCol = "light blue", postCol = "red", legendPos = "topright", maxNumberShifts=800, ...) 
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
  plot(c(0,0),xlim=c(0,maxNumberShifts),ylim=c(0,0.05),xlab='n shifts',type='n', xaxs='i',yaxs='i',xaxt='n',ylab='',yaxt='n')
  #barplot(prior[, 2], xlim=c(0,maxNumberShifts),names.arg = prior[, 1], ylim = c(0,0.2),type='n', border = NA, col = priorCol, 
  #        xlab = "n shifts", xaxs='i',yaxs='i',xaxt='n',yaxt='n', ...)
  lines(prior[, 2], xlim=c(0,maxNumberShifts),names.arg = prior[, 1], ylim = c(0,0.2), col = priorCol, 
        xlab = "n shifts", xaxs='i',yaxs='i',xaxt='n',yaxt='n', lwd=3)
  axis(1, at=seq(0,maxNumberShifts,by=25),cex.axis=.7)
  axis(2, at=seq(0,0.05,by=0.01),las=2,cex.axis=.7)
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
