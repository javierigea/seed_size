
library(BAMMtools)
library(coda)
library(plyr)

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
combine_multiple_BAMM_runs<-function(BAMMrunfiles,burnin.list){
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
  post_probs.div <- table(mcmcout$N_shifts) / nrow(mcmcout)
  names(post_probs.div)
  #plot the posterior distribution of Nrateshifts
  plot(names(post_probs.div),post_probs.div)
  plotPriormod(mcmcout,expectedNumberOfShifts = 500,burnin=0)
}




combine_multiple_BAMM_runs_output<-function(BAMMrunfiles,burnin.list,name){
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
  post_probs.div <- table(mcmcout$N_shifts) / nrow(mcmcout)
  names(post_probs.div)
  #plot the posterior distribution of Nrateshifts
  plot(names(post_probs.div),post_probs.div)
  first<-mcmcout$generation[1]
  mcmcout$generation<-mcmcout$generation-first
  #plotPriormod(mcmcout,expectedNumberOfShifts = 500,burnin=0)
  write.csv(mcmcout,file=paste('./',name,'_mcmc_out.txt',sep=''),quote=F,row.names=F)
}

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
  plot(c(0,0),xlim=c(0,800),ylim=c(0,0.05),xlab='n shifts',type='n', xaxs='i',yaxs='i',xaxt='n',ylab='',yaxt='n',border=NA)
  #barplot(prior[, 2], xlim=c(0,300),names.arg = prior[, 1], ylim = c(0,0.2),type='n', border = NA, col = priorCol, 
  #        xlab = "n shifts", xaxs='i',yaxs='i',xaxt='n',yaxt='n', ...)
  lines(prior[, 2], xlim=c(0,800),names.arg = prior[, 1], ylim = c(0,0.05), border = NA, col = priorCol, 
        xlab = "n shifts", xaxs='i',yaxs='i',xaxt='n',yaxt='n', lwd=3, ...)
  axis(1, at=seq(0,800,by=25),cex.axis=.7)
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


BAMMeventfiles<-c("/Users/javier/Desktop/hydrogen/seed_size/trait/event_data_trait_BAMM_Species_tree_500.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/trait/run500_event_data_trait_BAMM_Species_tree_500.txt","/Users/javier/Desktop/hydrogen/seed_size/trait/run500_2__event_data_trait_BAMM_Species_tree_500.txt","/Users/javier/Desktop/hydrogen/seed_size/trait/run500_3__event_data_trait_BAMM_Species_tree_500.txt")

burnin.list<-c(0.7,0.1,0.1,0.1)
#example: run plot_multiple_BAMM_runs first, select burnin and then run combine_multiple_BAMM_runs

###03April
#keep running 25 shifts
#keep running 250 shifts

######25 shifts
#this leaves the final 60 million generations
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_3_2_mcmc_out_div_BAMM_Species_tree_25.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_4_mcmc_out_div_BAMM_Species_tree_25.txt"),burnin=c(0.829,0))
concatenate_event_data(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_3_2_event_data_div_BAMM_Species_tree_25.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_4_event_data_div_BAMM_Species_tree_25.txt"),burnin=c(0.829,0),'25_1')

combine_multiple_BAMM_runs_output(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_3_mcmc_out_div_BAMM_Species_tree_25.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_4_mcmc_out_div_BAMM_Species_tree_25.txt"),c(0.5,0),name='25_div')
combine_multiple_BAMM_runs_output(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_3_mcmc_out_div_BAMM_Species_tree_50_laptop.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_4_mcmc_out_div_BAMM_Species_tree_50_laptop.txt"),c(0.5,0),name='50_div')
combine_multiple_BAMM_runs_output(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_3_mcmc_out_div_BAMM_Species_tree_100_laptop.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_4_mcmc_out_div_BAMM_Species_tree_100_laptop.txt"),c(0.5,0),name='100_div')
combine_multiple_BAMM_runs_output(c("/Users/javier/Desktop/hydrogen/seed_size/div/run3_mcmc_out_div_BAMM_Species_tree_250.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run4_mcmc_out_div_BAMM_Species_tree_250.txt"),c(0.5,0),name='250_div')


combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_3_mcmc_out_div_BAMM_Species_tree_50_laptop.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_4_mcmc_out_div_BAMM_Species_tree_50_laptop.txt"),burnin=c(0.25,0))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run3_mcmc_out_div_BAMM_Species_tree_250.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run4_mcmc_out_div_BAMM_Species_tree_250.txt"),c(0,0))

concatenate_event_data(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_3_mcmc_out_div_BAMM_Species_tree_50_laptop.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_4_mcmc_out_div_BAMM_Species_tree_50_laptop.txt"),burnin=c(0.2,0))
concatenate_event_data(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_3_event_data_div_BAMM_Species_tree_50_laptop.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_4_event_data_div_BAMM_Species_tree_50_laptop.txt"),c(0.2,0),'50_speciation')
#25 SHIFTS
plot_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/mcmc_out_div_BAMM_Species_tree_25_laptop.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/div/run_1_mcmc_out_div_BAMM_Species_tree_25.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_2_mcmc_out_div_BAMM_Species_tree_25.txt"))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_2_mcmc_out_div_BAMM_Species_tree_25.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_3_mcmc_out_div_BAMM_Species_tree_25.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_3_2_mcmc_out_div_BAMM_Species_tree_25.txt"),burnin=c(0.5,0,0))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_2_mcmc_out_div_BAMM_Species_tree_25.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_3_mcmc_out_div_BAMM_Species_tree_25.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_3_2_mcmc_out_div_BAMM_Species_tree_25.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_4_mcmc_out_div_BAMM_Species_tree_25.txt"),burnin=c(0.5,0,0,0))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_3_mcmc_out_div_BAMM_Species_tree_25.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_3_2_mcmc_out_div_BAMM_Species_tree_25.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_4_mcmc_out_div_BAMM_Species_tree_25.txt"),burnin=c(0,0,0))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_3_mcmc_out_div_BAMM_Species_tree_25.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_3_2_mcmc_out_div_BAMM_Species_tree_25.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_4_mcmc_out_div_BAMM_Species_tree_25.txt"),burnin=c(0,0,0))

combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_3_2_mcmc_out_div_BAMM_Species_tree_25.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_4_mcmc_out_div_BAMM_Species_tree_25.txt"),burnin=c(0.9,0))
#50 SHIFTS**CONVERGED
plot_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/mcmc_out_div_BAMM_Species_tree_50.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/div/run_1_mcmc_out_div_BAMM_Species_tree_50_laptop.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_2_mcmc_out_div_BAMM_Species_tree_50_laptop.txt"))
#this works
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_3_mcmc_out_div_BAMM_Species_tree_50_laptop.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_4_mcmc_out_div_BAMM_Species_tree_50_laptop.txt"),burnin=c(0.2,0))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_2_mcmc_out_div_BAMM_Species_tree_50_laptop.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_3_mcmc_out_div_BAMM_Species_tree_50_laptop.txt"),burnin=c(0.5,0))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_2_mcmc_out_div_BAMM_Species_tree_50_laptop.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_3_mcmc_out_div_BAMM_Species_tree_50_laptop.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_4_mcmc_out_div_BAMM_Species_tree_50_laptop.txt"),burnin=c(0.9,0,0))

####THIS IS BEST?
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_3_mcmc_out_div_BAMM_Species_tree_50_laptop.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_4_mcmc_out_div_BAMM_Species_tree_50_laptop.txt"),burnin=c(0.2,0))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_3_mcmc_out_div_BAMM_Species_tree_50_laptop.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_4_mcmc_out_div_BAMM_Species_tree_50_laptop.txt"),burnin=c(0.9,0))



#100 SHIFTS ****CONVERGED
plot_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/mcmc_out_div_BAMM_Species_tree_100_laptop.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/div/run_1_mcmc_out_div_BAMM_Species_tree_100_laptop.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_2_mcmc_out_div_BAMM_Species_tree_100_laptop.txt"))
#this works
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_2_mcmc_out_div_BAMM_Species_tree_100_laptop.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_3_mcmc_out_div_BAMM_Species_tree_100_laptop.txt"),burnin=c(0,0))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_2_mcmc_out_div_BAMM_Species_tree_100_laptop.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_3_mcmc_out_div_BAMM_Species_tree_100_laptop.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_4_mcmc_out_div_BAMM_Species_tree_100_laptop.txt"),burnin=c(0.9,0,0))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_3_mcmc_out_div_BAMM_Species_tree_100_laptop.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_4_mcmc_out_div_BAMM_Species_tree_100_laptop.txt"),burnin=c(0,0))

#250 SHIFTS ****CONVERGED
plot_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/mcmc_out_div_BAMM_Species_tree_250.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/div/run1_mcmc_out_div_BAMM_Species_tree_250.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run2_mcmc_out_div_BAMM_Species_tree_250.txt"))
#this is ok
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run1_mcmc_out_div_BAMM_Species_tree_250.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run2_mcmc_out_div_BAMM_Species_tree_250.txt"),c(0.8,0))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run2_mcmc_out_div_BAMM_Species_tree_250.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run3_mcmc_out_div_BAMM_Species_tree_250.txt"),c(0,0))
#this is also ok
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run1_mcmc_out_div_BAMM_Species_tree_250.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run2_mcmc_out_div_BAMM_Species_tree_250.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run3_mcmc_out_div_BAMM_Species_tree_250.txt"),c(0.8,0,0))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run3_mcmc_out_div_BAMM_Species_tree_250.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run4_mcmc_out_div_BAMM_Species_tree_250.txt"),c(0,0))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run2_mcmc_out_div_BAMM_Species_tree_250.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run3_mcmc_out_div_BAMM_Species_tree_250.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run4_mcmc_out_div_BAMM_Species_tree_250.txt"),c(0.5,0,0))


#500 SHIFTS
plot_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/mcmc_out_div_BAMM_Species_tree_500_mac.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/div/run_1_mcmc_out_div_BAMM_Species_tree_500.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_2_1st_mcmc_out_div_BAMM_Species_tree_500.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_2_2nd_mcmc_out_div_BAMM_Species_tree_500.txt"))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_2_1st_mcmc_out_div_BAMM_Species_tree_500.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_2_2nd_mcmc_out_div_BAMM_Species_tree_500.txt"),burnin=c(0.8,0))

plot_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/mcmc_out_trait_BAMM_Species_tree_750.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/run750_mcmc_out_trait_BAMM_Species_tree_750.txt"))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/mcmc_out_trait_BAMM_Species_tree_750.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/run750_mcmc_out_trait_BAMM_Species_tree_750.txt"),c(0.5,0.5))


plot_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/mcmc_out_trait_BAMM_Species_tree_750.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/run750_mcmc_out_trait_BAMM_Species_tree_750.txt"))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/mcmc_out_trait_BAMM_Species_tree_750.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/run750_mcmc_out_trait_BAMM_Species_tree_750.txt"),c(0.5,0.5))


plot_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/trait/mcmc_out_trait_BAMM_Species_tree_150.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/trait/run150_mcmc_out_trait_BAMM_Species_tree_150.txt"))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/trait/mcmc_out_trait_BAMM_Species_tree_150.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/trait/run150_mcmc_out_trait_BAMM_Species_tree_150.txt"),c(0.7,0.5))


plot_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/trait/mcmc_out_trait_BAMM_Species_tree_250.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/trait/run250_mcmc_out_trait_BAMM_Species_tree_250.txt"))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/trait/mcmc_out_trait_BAMM_Species_tree_250.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/trait/run250_mcmc_out_trait_BAMM_Species_tree_250.txt"),c(0.7,0.5))

plot_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/trait/mcmc_out_trait_BAMM_Species_tree_500.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/trait/run500_mcmc_out_trait_BAMM_Species_tree_500.txt","/Users/javier/Desktop/hydrogen/seed_size/trait/run500_2__mcmc_out_trait_BAMM_Species_tree_500.txt","/Users/javier/Desktop/hydrogen/seed_size/trait/run500_3__mcmc_out_trait_BAMM_Species_tree_500.txt"))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/trait/mcmc_out_trait_BAMM_Species_tree_500.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/trait/run500_mcmc_out_trait_BAMM_Species_tree_500.txt","/Users/javier/Desktop/hydrogen/seed_size/trait/run500_2__mcmc_out_trait_BAMM_Species_tree_500.txt","/Users/javier/Desktop/hydrogen/seed_size/trait/run500_3__mcmc_out_trait_BAMM_Species_tree_500.txt"),c(0.5,0.1,0.1,0.1))

plot_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/trait/mcmc_out_trait_BAMM_Species_tree_750.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/trait/run750_mcmc_out_trait_BAMM_Species_tree_750.txt"))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/trait/mcmc_out_trait_BAMM_Species_tree_750.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/trait/run750_mcmc_out_trait_BAMM_Species_tree_750.txt"),c(0.5,0.5))

plot_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/trait/mcmc_out_trait_BAMM_Species_tree_25.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/trait/run_25_mcmc_out_trait_BAMM_Species_tree_25.txt"))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/trait/mcmc_out_trait_BAMM_Species_tree_25.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/trait/run_25_mcmc_out_trait_BAMM_Species_tree_25.txt"),c(0.5,0.1))

plot_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/trait/mcmc_out_trait_BAMM_Species_tree_50.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/trait/run50_mcmc_out_trait_BAMM_Species_tree_50.txt"))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/trait/mcmc_out_trait_BAMM_Species_tree_50.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/trait/run50_mcmc_out_trait_BAMM_Species_tree_50.txt"),c(0.5,0.1))

plot_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/mcmc_out_div_BAMM_Species_tree_250_pc.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/div/mcmc_out_div_BAMM_Species_tree_250.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run1_mcmc_out_div_BAMM_Species_tree_250.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run2_mcmc_out_div_BAMM_Species_tree_250.txt"))

plot_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/mcmc_out_div_BAMM_Species_tree_50_laptop.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/div/run_1_mcmc_out_div_BAMM_Species_tree_50_laptop.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_2_mcmc_out_div_BAMM_Species_tree_50_laptop.txt"))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_1_mcmc_out_div_BAMM_Species_tree_50_laptop.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_2_mcmc_out_div_BAMM_Species_tree_50_laptop.txt"),c(0.5,0.1))

combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/mcmc_out_div_BAMM_Species_tree_250_pc.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/div/mcmc_out_div_BAMM_Species_tree_250.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run1_mcmc_out_div_BAMM_Species_tree_250.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run2_mcmc_out_div_BAMM_Species_tree_250.txt"),c(1,1,0.1,0.1))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run1_mcmc_out_div_BAMM_Species_tree_250.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run2_mcmc_out_div_BAMM_Species_tree_250.txt"),c(0.8,0.1))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run1_mcmc_out_div_BAMM_Species_tree_250.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run2_mcmc_out_div_BAMM_Species_tree_250.txt"),c(0.9,0.1))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/mcmc_out_div_BAMM_Species_tree_250_pc.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/div/run1_mcmc_out_div_BAMM_Species_tree_250.txt","/Users/javier/Desktop/hydrogen/seed_size/div/mcmc_out_div_BAMM_Species_tree_250.txt"),c(0.5,0.1,0.5))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/mcmc_out_div_BAMM_Species_tree_250_pc.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/div/run1_mcmc_out_div_BAMM_Species_tree_250_15M.txt"),c(0.5,0.1))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/trait/mcmc_out_trait_BAMM_Species_tree_50.txt" ,"/Users/javier/Desktop/hydrogen/seed_size/trait/run50_mcmc_out_trait_BAMM_Species_tree_50.txt"),c(0.5,0.1))

combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run1_mcmc_out_div_BAMM_Species_tree_250.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run2_mcmc_out_div_BAMM_Species_tree_250.txt"),c(0.8,0))
concatenate_event_data(c("/Users/javier/Desktop/hydrogen/seed_size/div/run1_event_data_div_BAMM_Species_tree_250.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run2_event_data_div_BAMM_Species_tree_250.txt"),c(0.8,0),'run1_2_250')


plot_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_1_mcmc_out_div_BAMM_Species_tree_25.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_2_mcmc_out_div_BAMM_Species_tree_25.txt"))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_1_mcmc_out_div_BAMM_Species_tree_25.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_2_mcmc_out_div_BAMM_Species_tree_25.txt"),c(0.9,0))

plot_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_1_mcmc_out_div_BAMM_Species_tree_500.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_2_mcmc_out_div_BAMM_Species_tree_500.txt"))
combine_multiple_BAMM_runs(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_1_mcmc_out_div_BAMM_Species_tree_500.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_2_mcmc_out_div_BAMM_Species_tree_500.txt"),c(0.99,0))

