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

####main BAMM analysis function: generates tables, runs strapp, etc
####it can also generate rate transformed trees (commented out now because it takes long to run)
#treefile: tree to analyse
#eventfile: event_data file
#traitdata:BAMMPhyndr_GenusSeed_trait_data.txt (no header;genus<\t>log10SeedWeight)
#burnin: burnin for STRAPP
#mode: 'diversification' or 'trait'
#name: name to append to analysis
analyse_BAMM_output<-function(treefile,eventfile,traitdata,burnin,mode,name){
  #load tree
  tree<-read.tree(treefile)
  #get event data
  print('Analysing eventfile')
  Div_edata <- getEventData(tree, eventfile, burnin = burnin, nsamples = 5000,type=mode)
  #number of shifts
  shift_probs <- summary(Div_edata)
  if (mode=='diversification'){
    #getTipRates
    rates<-getTipRates(Div_edata)
    div.rates<-getTipRates(Div_edata,returnNetDiv = T)
    mean.lambda.rates<-rates$lambda.avg
    mean.mu.rates<-rates$mu.avg
    mean.div.rates<-div.rates$netdiv.avg
    rates.df<-data.frame(names(mean.lambda.rates),mean.lambda.rates,mean.mu.rates,mean.div.rates)
    row.names(rates.df)<-NULL
    colnames(rates.df)[1]<-'names'
    trait_input<-read.table(traitdata,sep='\t')
    rates.df<-merge(rates.df,trait_input,by.x='names',by.y='V1')
    colnames(rates.df)[ncol(rates.df)]<-'log.Seed.Weight'    
    write.table(rates.df,file=paste('./output/tables/BAMM_diversification_rates_',name,'.txt',sep=''),quote=F,row.names=F,sep='\t')
  }
  if (mode=='trait'){
    #getTipRates
    rates<-getTipRates(Div_edata)
    mean.beta.rates<-rates$beta.avg
    rates.df<-data.frame(names(mean.beta.rates),mean.beta.rates)
    row.names(rates.df)<-NULL
    colnames(rates.df)[1]<-'names'
    #add SeedWeight info to table
    trait_input<-read.table(traitdata,sep='\t')
    rates.df<-merge(rates.df,trait_input,by.x='names',by.y='V1')
    colnames(rates.df)[ncol(rates.df)]<-'log.Seed.Weight'    
    write.table(rates.df,file=paste('./output/tables/BAMM_trait_rates_',name,'.txt',sep=''),quote=F,row.names=F,sep='\t')
  }
}  
#treefile: tree to analyse
#eventfile: event_data file
#burnin: burnin for event file
#mode:"trait" or "diversification"
#parameter:"speciation","extinction","ndr",'trait'
rate_scaled_trees<-function(treefile,eventfile,burnin,mode,parameter){
  #load tree
  tree<-read.tree(treefile)
  #get event data
  print('Analysing eventfile')
  Div_edata <- getEventData(tree, eventfile, burnin = burnin, nsamples = 5000,type=mode)
  #rate scaled trees
  print('Getting rate scaled tree')
  parameter.tree<-getMeanBranchLengthTree(Div_edata,rate=parameter)
  #write rate scaled trees
  file<-paste('./output/trees/',parameter,'_transformed','.tree',sep='')
  write.tree(parameter.tree$phy,file)
}



 # if (mode=='trait'){
    
##    
#  }
    #load trait data
##trait_input<-read.table(traitdata,sep='\t')
##trait_data<-trait_input$V2
##names(trait_data)<-trait_input$V1
###load trait data***CoefficientOfVariation
##trait_input.cov<-read.table(traitdatacov,sep='\t',header=T)
##trait_data.cov<-trait_input.cov$cov
##names(trait_data.cov)<-trait_input.sd$tip.label
###load trait rate evolution data
##traitrate_input<-read.table(traitratedata,sep='\t',header=T)
##trait_data.rate<-traitrate_input$mean.beta.rates
##names(trait_data.rate)<-traitrate_input$names
###add trait data to rates info,write table
##rates.trait.df<-merge(rates.df,trait_input,by.x='names',by.y='V1')
##colnames(rates.trait.df)[5]<-'log.Real.Seed.Weight'
##rates.trait.df.file<-paste('../output/tables/BAMMPhyndr_rates_traitdata_',name,'.txt',sep='')
##write.table(rates.trait.df,rates.trait.df.file,sep='\t',quote=F,row.names=F)
##plot(rates.trait.df$log.Real.Seed.Weight,log10(rates.trait.df$mean.lambda.rates))
###add trait data to rates info,write table
##rates.trait.rates.df<-merge(rates.trait.df,traitrate_input,by.x='names',by.y='names')
##rates.trait.rates.df$log.Real.Seed.Weight.y<-NULL
##colnames(rates.trait.rates.df)[5]<-'log10.Real.Seed.Weight'
##rates.trait.rates.df.file<-paste('../output/tables/BAMMPhyndr_rates_traitrates_data_',name,'.txt',sep='')
##write.table(rates.trait.rates.df,rates.trait.rates.df.file,sep='\t',quote=F,row.names=F)
###add stdev and cov data info, write table
##rates.trait.rates.variation.df<-merge(rates.trait.rates.df,trait_input.sd,by.x='names',by.y='tip.label')
##rates.trait.rates.variation.df<-merge(rates.trait.rates.variation.df,trait_input.cov,by.x='names',by.y='tip.label')
##rates.trait.rates.variation.df.file<-paste('../output/tables/BAMMPhyndr_rates_traitrates_variation_data_',name,'.txt',sep='')
##
##########rate scaled trees
###print('Getting lambda rate scaled tree')
###lambda.tree<-getMeanBranchLengthTree(Div_edata,rate='speciation')
###print('Getting mu rate scaled tree')
###mu.tree<-getMeanBranchLengthTree(Div_edata,rate = 'extinction')
###print('Getting div rate scaled trees')
###ndr.tree<-getMeanBranchLengthTree(Div_edata,rate  = 'ndr')
###write rate scaled trees
###lambdatreefile<-paste('../output/trees/lambda_transformed_',name,'.tree',sep='')
###mutreefile<-paste('../output/trees/mu_transformed_',name,'.tree',sep='')
###divtreefile<-paste('../output/trees/div_transformed_',name,'.tree',sep='')
###write.tree(lambda.tree$phy,lambdatreefile)
###write.tree(mu.tree$phy,mutreefile)
###write.tree(ndr.tree$phy,divtreefile)
##########strapp_seedmass
##print('Running strapp tests')
##strapp.seed.lambda<-traitDependentBAMM(Div_edata, traits = trait_data, 5000, rate = 'speciation', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
##strapp.seed.mu<-traitDependentBAMM(Div_edata, traits = trait_data, 5000, rate = 'extinction', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
##strapp.seed.div<-traitDependentBAMM(Div_edata, traits = trait_data, 5000, rate = 'net diversification', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
##strapp.seedrate.lambda<-traitDependentBAMM(Div_edata, traits = trait_data.rate, 5000, rate = 'speciation', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
##strapp.seedrate.mu<-traitDependentBAMM(Div_edata, traits = trait_data.rate, 5000, rate = 'extinction', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
##strapp.seedrate.div<-traitDependentBAMM(Div_edata, traits = trait_data.rate, 5000, rate = 'net diversification', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
##
##plot.strapp.histogram(strapp.seed.lambda,'lambda_seed')
##plot.strapp.histogram(strapp.seed.mu,'mu_seed')
##plot.strapp.histogram(strapp.seed.div,'div_seed')
##
###write table 
##strapp.seed.results<-data.frame(c('lambda','mu','div'),c(strapp.seed.lambda$estimate,strapp.seed.mu$estimate,strapp.seed.div$estimate),c(strapp.seed.lambda$p.value,strapp.seed.mu$p.value,strapp.seed.div$p.value))
##colnames(strapp.seed.results)<-c('parameter','estimate','pvalue')
##strappfile<-paste('../output/tables/Phyndr_strapp_',name,'.txt',sep='')
##write.table(strapp.seed.results,strappfile,sep='\t',quote=F,row.names=F)
##
####strapp.lambda.median
##strapp.lambda.diff<-abs(strapp.seed.lambda$obs.corr)-abs(strapp.seed.lambda$null)
##strapp.lambda.diff.df<-data.frame(strapp.lambda.diff,strapp.seed.lambda$gen)
##median.lambda<-median(strapp.lambda.diff.df$strapp.lambda.diff)
##median.strapp.lambda<-strapp.lambda.diff.df[which.min(abs(strapp.lambda.diff.df$strapp.lambda.diff - median.lambda)),]$strapp.seed.lambda.gen
###get the sample number where rho=median in the distribution
##median.strapp.lambda.eventdata <- cbind(Div_edata$edge, Div_edata$eventVectors[[median.strapp.lambda]])
##median.strapp.lambda.eventdata <- median.strapp.lambda.eventdata[median.strapp.lambda.eventdata[,2] <= length(tree$tip.label), ]
##median.strapp.lambda.eventdata <- data.frame(Div_edata$tip.label, median.strapp.lambda.eventdata[,3], Div_edata$tipLambda[[median.strapp.lambda]],Div_edata$tipMu[[median.strapp.lambda]],(Div_edata$tipLambda[[median.strapp.lambda]]-Div_edata$tipMu[[median.strapp.lambda]]))
##colnames(median.strapp.lambda.eventdata)<-c('tip.label','event','tip.lambda','tip.mu','tip.div')
##median.strapp.lambda.eventdata <- merge(median.strapp.lambda.eventdata,trait_input,by.x='tip.label',by.y='V1')
##colnames(median.strapp.lambda.eventdata)[6]<-'log10.Real.Seed.Weight'
##median.strapp.lambda.eventdata <- merge(median.strapp.lambda.eventdata,traitrate_input,by.x='tip.label',by.y='names')
##medianstrapplambdafile<-paste('../output/tables/BAMMPhyndr_medianstrappalambda_table_',name,'.txt',sep='')
##write.table(median.strapp.lambda.eventdata,medianstrapplambdafile,sep='\t',quote=F,row.names=F)
##
####strapp.mu.median
##strapp.mu.diff<-abs(strapp.seed.mu$obs.corr)-abs(strapp.seed.mu$null)
##strapp.mu.diff.df<-data.frame(strapp.mu.diff,strapp.seed.mu$gen)
##median.mu<-median(strapp.mu.diff.df$strapp.mu.diff)
##median.strapp.mu<-strapp.mu.diff.df[which.min(abs(strapp.mu.diff.df$strapp.mu.diff - median.mu)),]$strapp.seed.mu.gen
###get the sample number where rho=median in the distribution
##median.strapp.mu.eventdata <- cbind(Div_edata$edge, Div_edata$eventVectors[[median.strapp.mu]])
##median.strapp.mu.eventdata <- median.strapp.mu.eventdata[median.strapp.mu.eventdata[,2] <= length(tree$tip.label), ]
##median.strapp.mu.eventdata <- data.frame(Div_edata$tip.label, median.strapp.mu.eventdata[,3], Div_edata$tipLambda[[median.strapp.mu]],Div_edata$tipMu[[median.strapp.mu]],(Div_edata$tipLambda[[median.strapp.mu]]-Div_edata$tipMu[[median.strapp.mu]]))
##colnames(median.strapp.mu.eventdata)<-c('tip.label','event','tip.lambda','tip.mu','tip.div')
##median.strapp.mu.eventdata <- merge(median.strapp.mu.eventdata,trait_input,by.x='tip.label',by.y='V1')
##colnames(median.strapp.mu.eventdata)[6]<-'log10.Real.Seed.Weight'
##median.strapp.mu.eventdata <- merge(median.strapp.mu.eventdata,traitrate_input,by.x='tip.label',by.y='names')
##medianstrappmufile<-paste('../output/tables/BAMMPhyndr_medianstrappamu_table_',name,'.txt',sep='')
##write.table(median.strapp.mu.eventdata,medianstrappmufile,sep='\t',quote=F,row.names=F)
##
####strapp.div.median
##strapp.div.diff<-abs(strapp.seed.div$obs.corr)-abs(strapp.seed.div$null)
##strapp.div.diff.df<-data.frame(strapp.div.diff,strapp.seed.div$gen)
##median.div<-median(strapp.div.diff.df$strapp.div.diff)
##median.strapp.div<-strapp.div.diff.df[which.min(abs(strapp.div.diff.df$strapp.div.diff - median.div)),]$strapp.seed.div.gen
###get the sample number where rho=median in the distribution
##median.strapp.div.eventdata <- cbind(Div_edata$edge, Div_edata$eventVectors[[median.strapp.div]])
##median.strapp.div.eventdata <- median.strapp.div.eventdata[median.strapp.div.eventdata[,2] <= length(tree$tip.label), ]
##median.strapp.div.eventdata <- data.frame(Div_edata$tip.label, median.strapp.div.eventdata[,3], Div_edata$tipLambda[[median.strapp.div]],Div_edata$tipMu[[median.strapp.div]],(Div_edata$tipLambda[[median.strapp.div]]-Div_edata$tipMu[[median.strapp.div]]))
##colnames(median.strapp.div.eventdata)<-c('tip.label','event','tip.lambda','tip.mu','tip.div')
##median.strapp.div.eventdata <- merge(median.strapp.div.eventdata,trait_input,by.x='tip.label',by.y='V1')
##colnames(median.strapp.div.eventdata)[6]<-'log10.Real.Seed.Weight'
##median.strapp.div.eventdata <- merge(median.strapp.div.eventdata,traitrate_input,by.x='tip.label',by.y='names')
##medianstrappdivfile<-paste('../output/tables/BAMMPhyndr_medianstrappadiv_table_',name,'.txt',sep='')
##write.table(median.strapp.div.eventdata,medianstrappdivfile,sep='\t',quote=F,row.names=F)
##
##########strapp_seedmass****STDEV
##strapp.seed.stdev.lambda<-traitDependentBAMM(Div_edata, traits = trait_data.sd, 1000, rate = 'speciation', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
##strapp.seed.stdev.mu<-traitDependentBAMM(Div_edata, traits = trait_data.sd, 1000, rate = 'extinction', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
##strapp.seed.stdev.div<-traitDependentBAMM(Div_edata, traits = trait_data.sd, 1000, rate = 'net diversification', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
###write table 
##strapp.seed.st.dev.results<-data.frame(c('lambda','mu','div'),c(strapp.seed.stdev.lambda$estimate,strapp.seed.stdev.mu$estimate,strapp.seed.stdev.div$estimate),c(strapp.seed.stdev.lambda$p.value,strapp.seed.stdev.mu$p.value,strapp.seed.stdev.div$p.value))
##plot.strapp.histogram(strapp.seed.stdev.lambda,'lambda_sd_seed')
##plot.strapp.histogram(strapp.seed.stdev.mu,'mu_sd_seed')
##plot.strapp.histogram(strapp.seed.stdev.div,'div_sd_seed')
##
##colnames(strapp.seed.st.dev.results)<-c('parameter','estimate','pvalue')
##strappsdfile<-paste('../output/tables/Phyndr_strapp_sd_',name,'.txt',sep='')
##write.table(strapp.seed.st.dev.results,strappsdfile,sep='\t',quote=F,row.names=F)
##########strapp_seedmass****COV
##strapp.seed.cov.lambda<-traitDependentBAMM(Div_edata, traits = trait_data.cov, 5000, rate = 'speciation', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
##strapp.seed.cov.mu<-traitDependentBAMM(Div_edata, traits = trait_data.cov, 5000, rate = 'extinction', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
##strapp.seed.cov.div<-traitDependentBAMM(Div_edata, traits = trait_data.cov, 5000, rate = 'net diversification', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
###write table 
##plot.strapp.histogram(strapp.seed.cov.lambda,'lambda_cov_seed')
##plot.strapp.histogram(strapp.seed.cov.mu,'mu_cov_seed')
##plot.strapp.histogram(strapp.seed.cov.div,'div_cov_seed')
##strapp.seed.cov.results<-data.frame(c('lambda','mu','div'),c(strapp.seed.cov.lambda$estimate,strapp.seed.cov.mu$estimate,strapp.seed.cov.div$estimate),c(strapp.seed.cov.lambda$p.value,strapp.seed.cov.mu$p.value,strapp.seed.cov.div$p.value))
##colnames(strapp.seed.cov.results)<-c('parameter','estimate','pvalue')
##strappcovfile<-paste('../output/tables/Phyndr_strapp_cov_',name,'.txt',sep='')
##write.table(strapp.seed.cov.results,strappcovfile,sep='\t',quote=F,row.names=F)
  ##
  ####}
#####FROM THIS POINT ON IT'S ONLY COHORT RELATED STUFF. USEFUL BUT REMOVE FROM CLEAN VERSIONS  
###   
###   
###   #####best shift + get cohorts info
###   print('Getting Best Shift Configuration')
###   best <- getBestShiftConfiguration(Div_edata, expectedNumberOfShifts=nshifts)
###   #build dataframe with cohort adscription, rates, tiplabels and trait data
###   best.df<-cbind(best$edge, best$eventVectors[[1]])
###   best.df <- best.df[best.df[,2] <= length(tree$tip.label), ]
###   #add info on tip parameters
###   best.df <- data.frame(best$tip.label, best.df[,3], best$tipLambda[[1]],best$tipMu[[1]],(best$tipLambda[[1]]-best$tipMu[[1]]))
###   colnames(best.df)<-c('tip.label','event','tip.lambda','tip.mu','tip.div')
###   #add info on trait data to tips
###   best.df <- merge(best.df,trait_input,by.x='tip.label',by.y='V1')
###   colnames(best.df)[6]<-'log10.Real.Seed.Weight'
###   best.df <- merge(best.df,traitrate_input,by.x='tip.label',by.y='names')
###   bestfile<-paste('../output/tables/BAMMPhyndr_best_table_',name,'.txt',sep='')
###   write.table(best.df,bestfile,sep='\t',quote=F,row.names=F)
###   #get the correlation across the whole posterior of samples
###   nsamples<-length(Div_edata$eventVectors)
###   cor.posterior.lambda<-matrix(NA,nrow=nsamples,ncol=3)
###   cor.posterior.mu<-matrix(NA,nrow=nsamples,ncol=3)
###   cor.posterior.div<-matrix(NA,nrow=nsamples,ncol=3)
###   spearman.posterior.lambda<-matrix(NA,nrow=nsamples,ncol=2)
###   spearman.posterior.mu<-matrix(NA,nrow=nsamples,ncol=2)
###   spearman.posterior.div<-matrix(NA,nrow=nsamples,ncol=2)
###   print('Running correlations through posterior')
###   for (i in 1:nsamples){
###     cat(i, '\n')	
###     eventdata <- cbind(Div_edata$edge, Div_edata$eventVectors[[i]])
###     eventdata <- eventdata[eventdata[,2] <= length(tree$tip.label), ]
###     eventdata <- data.frame(Div_edata$tip.label, eventdata[,3], Div_edata$tipLambda[[i]],Div_edata$tipMu[[i]],(Div_edata$tipLambda[[i]]-Div_edata$tipMu[[i]]))
###     colnames(eventdata)<-c('tip.label','event','tip.lambda','tip.mu','tip.div')
###     eventdata <- merge(eventdata,trait_input,by.x='tip.label',by.y='V1')
###     colnames(eventdata)[6]<-'log10.Real.Seed.Weight'
###     eventdata$event<-as.factor(eventdata$event)
###     table.events<-table(eventdata$event)
###     table.events<-data.frame(names(table.events),unname(table.events))
###     table.events$Var1<-NULL
###     colnames(table.events)<-c('event','count')
###     eventdata.lambda<-aggregate(tip.lambda~event, data=eventdata, mean)
###     eventdata.lambda<-merge(eventdata.lambda,table.events,by.x='event',by.y='event')
###     eventdata.mu<-aggregate(tip.mu~event, data=eventdata, mean)
###     eventdata.mu<-merge(eventdata.mu,table.events,by.x='event',by.y='event')
###     eventdata.div<-aggregate(tip.div~event, data=eventdata, mean)
###     eventdata.div<-merge(eventdata.div,table.events,by.x='event',by.y='event')
###     eventdata.seed<-aggregate(log10.Real.Seed.Weight~event, data=eventdata, mean)
###     eventdata.seed<-merge(eventdata.seed,table.events,by.x='event',by.y='event')
###     
###     p.lambda<-lm(log10(eventdata.lambda$tip.lambda)~eventdata.seed$log10.Real.Seed.Weight)
###     p.mu<-lm(log10(eventdata.mu$tip.mu)~eventdata.seed$log10.Real.Seed.Weight)
###     p.div<-lm(log10(eventdata.div$tip.div)~eventdata.seed$log10.Real.Seed.Weight)
###     
###     spearman.lambda<-cor.test(log10(eventdata.lambda$tip.lambda),eventdata.seed$log10.Real.Seed.Weight,method='spearman')
###     spearman.mu<-cor.test(log10(eventdata.mu$tip.mu),eventdata.seed$log10.Real.Seed.Weight,method='spearman')
###     spearman.div<-cor.test(log10(eventdata.div$tip.div),eventdata.seed$log10.Real.Seed.Weight,method='spearman')
###     
###     cor.posterior.lambda[i,]<-c(i,p.lambda$coefficients[[1]],p.lambda$coefficients[[2]])
###     cor.posterior.mu[i,]<-c(i,p.mu$coefficients[[1]],p.mu$coefficients[[2]])
###     cor.posterior.div[i,]<-c(i,p.div$coefficients[[1]],p.div$coefficients[[2]])
###     
###     spearman.posterior.lambda[i,]<-c(i,as.numeric(spearman.lambda$estimate))
###     spearman.posterior.mu[i,]<-c(i,as.numeric(spearman.mu$estimate))
###     spearman.posterior.div[i,]<-c(i,as.numeric(spearman.div$estimate))
###     
###   }
###   spearman.posterior.lambda<-as.data.frame(spearman.posterior.lambda)
###   colnames(spearman.posterior.lambda)<-c('name','rho')
###   spearmanlambdafile<-paste('../output/tables/spearman_posterior_lambda_correlations_',name,'.txt',sep='')
###   write.table(spearman.posterior.lambda,spearmanlambdafile,sep='\t',quote=F,row.names=F)
###   spearman.posterior.mu<-as.data.frame(spearman.posterior.mu)
###   colnames(spearman.posterior.mu)<-c('name','rho')
###   spearmanmufile<-paste('../output/tables/spearman_posterior_mu_correlations_',name,'.txt',sep='')
###   write.table(spearman.posterior.mu,spearmanmufile,sep='\t',quote=F,row.names=F)
###   spearman.posterior.div<-as.data.frame(spearman.posterior.div)
###   colnames(spearman.posterior.div)<-c('name','rho')
###   spearmandivfile<-paste('../output/tables/spearman_posterior_div_correlations_',name,'.txt',sep='')
###   write.table(spearman.posterior.div,spearmandivfile,sep='\t',quote=F,row.names=F)
###   
###   ########cohort for sample where rho=median in the distribution
###   #get the sample number where rho=median in the distribution
###   median<-median(spearman.posterior.lambda$rho)
###   median.value<-spearman.posterior.lambda[spearman.posterior.lambda$rho==median,]$name
###   #get the sample number where rho=median in the distribution
###   median.eventdata <- cbind(Div_edata$edge, Div_edata$eventVectors[[median.value]])
###   median.eventdata <- median.eventdata[median.eventdata[,2] <= length(tree$tip.label), ]
###   median.eventdata <- data.frame(Div_edata$tip.label, median.eventdata[,3], Div_edata$tipLambda[[median.value]],Div_edata$tipMu[[median.value]],(Div_edata$tipLambda[[median.value]]-Div_edata$tipMu[[median.value]]))
###   colnames(median.eventdata)<-c('tip.label','event','tip.lambda','tip.mu','tip.div')
###   median.eventdata <- merge(median.eventdata,trait_input,by.x='tip.label',by.y='V1')
###   colnames(median.eventdata)[6]<-'log10.Real.Seed.Weight'
###   medianfile<-paste('../output/tables/BAMMPhyndr_median_table_',name,'.txt',sep='')
###   write.table(median.eventdata,medianfile,sep='\t',quote=F,row.names=F)
###   
###   ###########SEED MASS RATE ANALYSIS
###   ########strapp_seedmassrates
###   print('Running strapp rate tests')
###   strapp.seedrate.lambda<-traitDependentBAMM(Div_edata, traits = trait_data.rate, 1000, rate = 'speciation', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
###   strapp.seedrate.mu<-traitDependentBAMM(Div_edata, traits = trait_data.rate, 1000, rate = 'extinction', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
###   strapp.seedrate.div<-traitDependentBAMM(Div_edata, traits = trait_data.rate, 1000, rate = 'net diversification', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
###   plot.strapp.histogram(strapp.seedrate.lambda,'lambda_rate')
###   plot.strapp.histogram(strapp.seedrate.mu,'mu_rate')
###   plot.strapp.histogram(strapp.seedrate.div,'div_rate')
###   #write table 
###   strapp.seedrate.results<-data.frame(c('lambda','mu','div'),c(strapp.seedrate.lambda$estimate,strapp.seedrate.mu$estimate,strapp.seedrate.div$estimate),c(strapp.seedrate.lambda$p.value,strapp.seedrate.mu$p.value,strapp.seedrate.div$p.value))
###   colnames(strapp.seedrate.results)<-c('parameter','estimate','pvalue')
###   strappratefile<-paste('../output/tables/Phyndr_strapp_rate_',name,'.txt',sep='')
###   write.table(strapp.seedrate.results,strappratefile,sep='\t',quote=F,row.names=F)
###   
###}}
############MSC INSTEAD OF BEST SHIFT
######msc_tree <- maximumShiftCredibility(Div_edata)
######msceventdata <- cbind(Div_edata$edge, Div_edata$eventVectors[[msc_tree$sampleindex]])
######msceventdata <- msceventdata[msceventdata[,2] <= length(tree$tip.label), ]
######msceventdata <- data.frame(Div_edata$tip.label, msceventdata[,3], Div_edata$tipLambda[[msc_tree$sampleindex]],Div_edata$tipMu[[msc_tree$sampleindex]],(Div_edata$tipLambda[[msc_tree$sampleindex]]-Div_edata$tipMu[[msc_tree$sampleindex]]))
######colnames(msceventdata)<-c('tip.label','event','tip.lambda','tip.mu','tip.div')
######msceventdata <- merge(msceventdata,trait_input,by.x='tip.label',by.y='V1')
######colnames(msceventdata)[6]<-'log10.Real.Seed.Weight'
######msceventdata$event<-as.factor(msceventdata$event)
######table.events<-table(msceventdata$event)
######table.events<-data.frame(names(table.events),unname(table.events))
######table.events$Var1<-NULL
######colnames(table.events)<-c('event','count')
######msceventdata.lambda<-aggregate(tip.lambda~event, data=msceventdata, mean)
######msceventdata.lambda<-merge(msceventdata.lambda,table.events,by.x='event',by.y='event')
######msceventdata.mu<-aggregate(tip.mu~event, data=msceventdata, mean)
######msceventdata.mu<-merge(msceventdata.mu,table.events,by.x='event',by.y='event')
######msceventdata.div<-aggregate(tip.div~event, data=msceventdata, mean)
######msceventdata.div<-merge(msceventdata.div,table.events,by.x='event',by.y='event')
######msceventdata.seed<-aggregate(log10.Real.Seed.Weight~event, data=msceventdata, mean)
######msceventdata.seed<-merge(msceventdata.seed,table.events,by.x='event',by.y='event')
######spearman.msc.lambda<-cor.test(log10(msceventdata.lambda$tip.lambda),msceventdata.seed$log10.Real.Seed.Weight,method='spearman')
######spearman.msc.mu<-cor.test(log10(msceventdata.mu$tip.mu),msceventdata.seed$log10.Real.Seed.Weight,method='spearman')
######spearman.msc.div<-cor.test(log10(msceventdata.div$tip.div),msceventdata.seed$log10.Real.Seed.Weight,method='spearman')
