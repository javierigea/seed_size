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
library(ape)
library(geiger)
library(Taxonstand)
library(BAMMtools)
library(coda)
library(taxonlookup)
library(phytools)
library(diversitree)
library(phyndr)
library(gdata)
library(MonoPhy)
library(plyr)



annual_vs_perennial<-function(treefile,traitfile,cvaluetablefile,lifecycletablefile){
  tree<-read.tree(treefile)
  #load the trait dataset
  trait<-read.table(traitfile,header=T,sep='\t',stringsAsFactors = F)
  tree.trait<-drop.tip(tree,setdiff(tree$tip.label,trait$Species))
  #run pgls
  tree.trait$node.label<-NULL
  comp.data<-comparative.data(tree.trait,trait,names.col = 'Species')
  #not working with caper
  #pgls.Seed.LifeCycle<-pgls(log10.SeedWeight~cycle,data=comp.data, lambda="ML")
  pgls.seedLifeCycle<-gls(log10.SeedWeight ~ cycle, correlation = corPagel(1,phy = tree.trait,fixed=FALSE),data = trait, method = "ML")
  
  pgls.seedwoody<-gls(log10.SeedWeight ~ woodiness.state, correlation = corPagel(1,phy = tree.trait,fixed=FALSE),data = trait, method = "ML")
  summary(pgls.Seed.Annuality)
  ######plot seed mass in annuals and perennials
  seedCvaluelifecycle.strict$Life.Cycle<-NA
  seedCvaluelifecycle.strict[seedCvaluelifecycle.strict$Annuality==0,]$Life.Cycle<-'perennial'
  seedCvaluelifecycle.strict[seedCvaluelifecycle.strict$Annuality==1,]$Life.Cycle<-'annual'
  seedCvaluelifecycle.strict$Life.Cycle<-factor(seedCvaluelifecycle.strict$Life.Cycle)
  pdf(file='./output/plots/annual_vs_perennial_boxplot.pdf',paper='a4')
  boxplot(log10.Seed.Weight~as.factor(Life.Cycle),data=seedCvaluelifecycle.strict,horizontal=T, las=1,boxwex=0.3,cex=0.3,col=c('light grey','white'), main='Life Cycle and Seed Mass',xaxt='n',yaxt='n',xlab='Genus Seed Weight (g/seed)',notch=T)
  legend(x='topright', fill = c('light grey', 'white'), legend = c('annual','perennial'), bty = 'n',cex=.7)
  ticks <- seq(-6, 3, by=1)
  labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
  axis(1, at=ticks, labels=labels,cex.axis=.7)
  dev.off()
}



#seedtablefile BAMMPhyndr_GenusSeed_trait_data.txt
#cvaluetablefile KewCvaluedata_GenusLevel_taxonlookup.txt
#lifecycletablefile KewLifeCycle_GenusLevel.txt
run_seedlifecycleCvaluewoodiness_STRAPP<-function(treefile,eventfile,burnin,traitfile){
  tree<-read.tree(treefile)
  #get event data
  if (class(eventfile)!='bammdata'){
    Div_edata <- getEventData(tree, eventfile, burnin = burnin, nsamples = 5000,type='diversification')
  }else if (class(eventfile)=='bammdata'){
    Div_edata<-eventfile
  }
  #load traitdata
  trait<-read.table(traitfile,header=T,sep='\t',stringsAsFactors = F)
  #change lifecycle to categorical
  #lifecycle: perennial-0; annual-1
  trait[trait$cycle=='annual',]$cycle<-1
  trait[trait$cycle=='perennial',]$cycle<-0
  
  #extracting subtreeBAMM
  subtree<-subtreeBAMM(Div_edata,tips=trait$Species)
  #all four analysis
  all.Cvalue_data<-trait$C.value
  names(all.Cvalue_data)<-trait$Species
  all.seed_data<-trait$log10.SeedWeight
  names(all.seed_data)<-trait$Species
  all.LifeCycle_data<-trait$cycle
  names(all.LifeCycle_data)<-trait$Species
  all.LifeCycle_data<-as.numeric(trait$cycle)
  names(all.LifeCycle_data)<-trait$Species
  all.Woodiness_data<-trait$woodiness.state
  names(all.Woodiness_data)<-trait$Species
  
  #running STRAPP for all dataset
  #seed mass STRAPP
  strapp.all.seed.lambda<-traitDependentBAMM(subtree, traits = all.seed_data, 1000, rate = 'speciation', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.seed.mu<-traitDependentBAMM(subtree, traits = all.seed_data, 1000, rate = 'extinction', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.seed.div<-traitDependentBAMM(subtree, traits = all.seed_data, 1000, rate = 'net diversification', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.seed.lambda.diff<-abs(strapp.all.seed.lambda$obs.corr)-abs(strapp.all.seed.lambda$null)
  strapp.all.seed.mu.diff<-abs(strapp.all.seed.mu$obs.corr)-abs(strapp.all.seed.mu$null)
  strapp.all.seed.div.diff<-abs(strapp.all.seed.div$obs.corr)-abs(strapp.all.seed.div$null)
  #Cvalue STRAPP
  strapp.all.Cvalue.lambda<-traitDependentBAMM(subtree, traits = all.Cvalue_data, 1000, rate = 'speciation', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.Cvalue.mu<-traitDependentBAMM(subtree, traits = all.Cvalue_data, 1000, rate = 'extinction', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.Cvalue.div<-traitDependentBAMM(subtree, traits = all.Cvalue_data, 1000, rate = 'net diversification', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.Cvalue.lambda.diff<-abs(strapp.all.Cvalue.lambda$obs.corr)-abs(strapp.all.Cvalue.lambda$null)
  strapp.all.Cvalue.mu.diff<-abs(strapp.all.Cvalue.mu$obs.corr)-abs(strapp.all.Cvalue.mu$null)
  strapp.all.Cvalue.div.diff<-abs(strapp.all.Cvalue.div$obs.corr)-abs(strapp.all.Cvalue.div$null)
  #LifeCycle STRAPP
  strapp.all.LifeCycle.lambda<-traitDependentBAMM(subtree, traits = all.LifeCycle_data, 1000, rate = 'speciation', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.LifeCycle.mu<-traitDependentBAMM(subtree, traits = all.LifeCycle_data, 1000, rate = 'extinction', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.LifeCycle.div<-traitDependentBAMM(subtree, traits = all.LifeCycle_data, 1000, rate = 'net diversification', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.LifeCycle.lambda.diff<-abs(strapp.all.LifeCycle.lambda$obs.corr)-abs(strapp.all.LifeCycle.lambda$null)
  strapp.all.LifeCycle.mu.diff<-abs(strapp.all.LifeCycle.mu$obs.corr)-abs(strapp.all.LifeCycle.mu$null)
  strapp.all.LifeCycle.div.diff<-abs(strapp.all.LifeCycle.div$obs.corr)-abs(strapp.all.LifeCycle.div$null)
  #Woodiness STRAPP
  strapp.all.Woodiness.lambda<-traitDependentBAMM(subtree, traits = all.Woodiness_data, 1000, rate = 'speciation', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.Woodiness.mu<-traitDependentBAMM(subtree, traits = all.Woodiness_data, 1000, rate = 'extinction', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.Woodiness.div<-traitDependentBAMM(subtree, traits = all.Woodiness_data, 1000, rate = 'net diversification', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.Woodiness.lambda.diff<-abs(strapp.all.Woodiness.lambda$obs.corr)-abs(strapp.all.Woodiness.lambda$null)
  strapp.all.Woodiness.mu.diff<-abs(strapp.all.Woodiness.mu$obs.corr)-abs(strapp.all.Woodiness.mu$null)
  strapp.all.Woodiness.div.diff<-abs(strapp.all.Woodiness.div$obs.corr)-abs(strapp.all.Woodiness.div$null)
  
  #table for seed results
  seed.results<-data.frame(cbind(c(strapp.all.seed.lambda$estimate,strapp.all.seed.mu$estimate,strapp.all.seed.div$estimate),c(strapp.all.seed.lambda$p.value,strapp.all.seed.mu$p.value,strapp.all.seed.div$p.value)))
  row.names(seed.results)<-c('lambda','mu','div')
  colnames(seed.results)<-c('estimate','p-value')
  #table for Cvalue results
  Cvalue.results<-data.frame(cbind(c(strapp.all.Cvalue.lambda$estimate,strapp.all.Cvalue.mu$estimate,strapp.all.Cvalue.div$estimate),c(strapp.all.Cvalue.lambda$p.value,strapp.all.Cvalue.mu$p.value,strapp.all.Cvalue.div$p.value)))
  row.names(Cvalue.results)<-c('lambda','mu','div')
  colnames(Cvalue.results)<-c('estimate','p-value')
  #table for LifeCycle results
  LifeCycle.results<-data.frame(cbind(c(strapp.all.LifeCycle.lambda$estimate,strapp.all.LifeCycle.mu$estimate,strapp.all.LifeCycle.div$estimate),c(strapp.all.LifeCycle.lambda$p.value,strapp.all.LifeCycle.mu$p.value,strapp.all.LifeCycle.div$p.value)))
  row.names(LifeCycle.results)<-c('lambda','mu','div')
  colnames(LifeCycle.results)<-c('estimate','p-value')
  #table for Woodiness results
  Woodiness.results<-data.frame(cbind(c(strapp.all.Woodiness.lambda$estimate,strapp.all.Woodiness.mu$estimate,strapp.all.Woodiness.div$estimate),c(strapp.all.Woodiness.lambda$p.value,strapp.all.Woodiness.mu$p.value,strapp.all.Woodiness.div$p.value)))
  row.names(Woodiness.results)<-c('lambda','mu','div')
  colnames(Woodiness.results)<-c('estimate','p-value')
  
  #seedcvaluelifecyclewoodiness.df<-cbind(seed.results,Cvalue.results,LifeCycle.results,Woodiness.results)
  seedcvaluelifecyclewoodiness.df<-cbind(seed.results,Cvalue.results,Woodiness.results)
  #colnames(seedcvaluelifecyclewoodiness.df)<-c('seedmass.estimate','seedmass.pvalue','Cvalue.estimate','Cvalue.pvalue','LifeCycle.estimate','LifeCycle.pvalue','Woodiness.estimate','Woodiness.pvalue')
  colnames(seedcvaluelifecyclewoodiness.df)<-c('seedmass.estimate','seedmass.pvalue','Cvalue.estimate','Cvalue.pvalue','Woodiness.estimate','Woodiness.pvalue')
  seedcvaluelifecyclewoodiness.df$parameter<-row.names(seedcvaluelifecyclewoodiness.df)
  row.names(seedcvaluelifecyclewoodiness.df)<-NULL
  seedcvaluelifecyclewoodiness.df<-seedcvaluelifecyclewoodiness.df[c(ncol(seedcvaluelifecyclewoodiness.df),1:(ncol(seedcvaluelifecyclewoodiness.df)-1))]
  write.table(seedcvaluelifecyclewoodiness.df,'./output/tables/seedCvalueLifeCycle_STRAPP.txt',quote=F,row.names=F,sep='\t')
  #write plot for lambda
  pdf(file='./output/plots/lambda_3trait_strapp_density.pdf',paper='a4')
  plot(density(strapp.all.seed.lambda.diff),ylim=c(0,10),xlim=c(-0.4,0.3),main=paste('SeedMass,C.value,Life Cycle (n=',length(all.Cvalue_data),')',sep=''),ylab='',xlab='absolute.diff in rho(obs-null)',col='blue')
  lines(density(strapp.all.Cvalue.lambda.diff),col='red')
  lines(density(strapp.all.LifeCycle.lambda.diff),col='dark green')
  lines(density(strapp.all.Woodiness.lambda.diff),col='purple')
  abline(v=0,lty=2)
  abline(v=mean(strapp.all.seed.lambda.diff),col='blue',lty=2)
  abline(v=mean(strapp.all.Cvalue.lambda.diff),col='red',lty=2)
  abline(v=mean(strapp.all.LifeCycle.lambda.diff),col='dark green',lty=2)
  abline(v=mean(strapp.all.Woodiness.lambda.diff),col='purple',lty=2)
  legend('topright',lty=1,col=c('red','blue','dark green','purple'),c('C.value','SeedMass','Life Cycle','Woodiness'),cex=.5,bty='n')
  dev.off()
  #write plot for mu
  pdf(file='./output/plots/mu_3trait_strapp_density.pdf',paper='a4')
  plot(density(strapp.all.seed.mu.diff),ylim=c(0,10),xlim=c(-0.4,0.3),main=paste('SeedMass,C.value,Life Cycle (n=',length(all.Cvalue_data),')',sep=''),ylab='',xlab='absolute.diff in rho(obs-null)',col='blue')
  lines(density(strapp.all.Cvalue.mu.diff),col='red')
  lines(density(strapp.all.LifeCycle.mu.diff),col='dark green')
  lines(density(strapp.all.Woodiness.mu.diff),col='purple')
  abline(v=0,lty=2)
  abline(v=mean(strapp.all.seed.mu.diff),col='blue',lty=2)
  abline(v=mean(strapp.all.Cvalue.mu.diff),col='red',lty=2)
  abline(v=mean(strapp.all.LifeCycle.mu.diff),col='dark green',lty=2)
  abline(v=mean(strapp.all.Woodiness.mu.diff),col='purple',lty=2)
  legend('topright',lty=1,col=c('red','blue','dark green','purple'),c('C.value','SeedMass','Life Cycle','Woodiness'),cex=.5,bty='n')
  dev.off()
  #write plot for div
  pdf(file='./output/plots/div_3trait_strapp_density.pdf',paper='a4')
  plot(density(strapp.all.seed.div.diff),ylim=c(0,10),xlim=c(-0.4,0.3),main=paste('SeedMass,C.value,Life Cycle (n=',length(all.Cvalue_data),')',sep=''),ylab='',xlab='absolute.diff in rho(obs-null)',col='blue')
  lines(density(strapp.all.Cvalue.div.diff),col='red')
  lines(density(strapp.all.LifeCycle.div.diff),col='dark green')
  lines(density(strapp.all.Woodiness.div.diff),col='purple')
  abline(v=0,lty=2)
  abline(v=mean(strapp.all.seed.div.diff),col='blue',lty=2)
  abline(v=mean(strapp.all.Cvalue.div.diff),col='red',lty=2)
  abline(v=mean(strapp.all.LifeCycle.div.diff),col='dark green',lty=2)
  abline(v=mean(strapp.all.Woodiness.div.diff),col='purple',lty=2)
  legend('topright',lty=1,col=c('red','blue','dark green','purple'),c('C.value','SeedMass','Life Cycle','Woodiness'),cex=.5,bty='n')
  dev.off()
}
