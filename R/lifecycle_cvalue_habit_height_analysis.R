library(ape)
library(Taxonstand)
library(taxonlookup)
library(phytools)
library(diversitree)
library(geiger)
library(BAMMtools)
library(nlme)
library(phyndr)
library(MonoPhy)
library(picante)
library(RPANDA)
library(parallel)
library(PBD)
library(caper)
library(MuMIn)

################this function runs the STRAPP analysis in a dataset with seed,Cvalue,lifecycle and woodiness data (kewsid_lifecycle_woodiness_TPL.txt)
#traitfile = kewsid_lifecycle_woodiness_TPL.txt (produced by sourcing this very same .R)
run_seedlifecycleCvaluewoodiness_STRAPP<-function(treefile,eventfile,burnin,traitfile){
  tree<-read.tree(treefile)
  #get event data
  if (class(eventfile)!='bammdata'){
    Div_edata <- getEventData(tree, eventfile, burnin = burnin, nsamples = 1000,type='diversification')
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
  seedcvaluelifecyclewoodiness.df<-cbind(seed.results,Cvalue.results,LifeCycle.results,Woodiness.results)
  #colnames(seedcvaluelifecyclewoodiness.df)<-c('seedmass.estimate','seedmass.pvalue','Cvalue.estimate','Cvalue.pvalue','LifeCycle.estimate','LifeCycle.pvalue','Woodiness.estimate','Woodiness.pvalue')
  colnames(seedcvaluelifecyclewoodiness.df)<-c('seedmass.estimate','seedmass.pvalue','Cvalue.estimate','Cvalue.pvalue','LifeCycle.estimate','LifeCycle.pvalue','Woodiness.estimate','Woodiness.pvalue')
  seedcvaluelifecyclewoodiness.df$parameter<-row.names(seedcvaluelifecyclewoodiness.df)
  row.names(seedcvaluelifecyclewoodiness.df)<-NULL
  seedcvaluelifecyclewoodiness.df<-seedcvaluelifecyclewoodiness.df[c(ncol(seedcvaluelifecyclewoodiness.df),1:(ncol(seedcvaluelifecyclewoodiness.df)-1))]
  write.table(seedcvaluelifecyclewoodiness.df,'./output/tables/seedCvalueLifeCycle_STRAPP.txt',quote=F,row.names=F,sep='\t')
  #write plot for lambda
  pdf(file='./output/plots/lambda_3trait_strapp_density.pdf',paper='a4')
  plot(density(strapp.all.seed.lambda.diff),ylim=c(0,10),xlim=c(-0.4,0.3),main=paste('SeedMass,C.value,Life Cycle (n=',length(all.Cvalue_data),')',sep=''),ylab='',xlab='absolute.diff in rho(obs-null)',col='cyan')
  lines(density(strapp.all.Cvalue.lambda.diff),col='red')
  lines(density(strapp.all.LifeCycle.lambda.diff),col='dark green')
  lines(density(strapp.all.Woodiness.lambda.diff),col='purple')
  abline(v=0,lty=2)
  abline(v=mean(strapp.all.seed.lambda.diff),col='cyan',lty=2)
  abline(v=mean(strapp.all.Cvalue.lambda.diff),col='red',lty=2)
  abline(v=mean(strapp.all.LifeCycle.lambda.diff),col='dark green',lty=2)
  abline(v=mean(strapp.all.Woodiness.lambda.diff),col='purple',lty=2)
  legend('topright',lty=1,col=c('red','cyan','dark green','purple'),c('C.value','SeedMass','Life Cycle','Woodiness'),cex=.5,bty='n')
  dev.off()
  #write plot for mu
  pdf(file='./output/plots/mu_3trait_strapp_density.pdf',paper='a4')
  plot(density(strapp.all.seed.mu.diff),ylim=c(0,10),xlim=c(-0.4,0.3),main=paste('SeedMass,C.value,Life Cycle (n=',length(all.Cvalue_data),')',sep=''),ylab='',xlab='absolute.diff in rho(obs-null)',col='cyan')
  lines(density(strapp.all.Cvalue.mu.diff),col='red')
  lines(density(strapp.all.LifeCycle.mu.diff),col='dark green')
  lines(density(strapp.all.Woodiness.mu.diff),col='purple')
  abline(v=0,lty=2)
  abline(v=mean(strapp.all.seed.mu.diff),col='cyan',lty=2)
  abline(v=mean(strapp.all.Cvalue.mu.diff),col='red',lty=2)
  abline(v=mean(strapp.all.LifeCycle.mu.diff),col='dark green',lty=2)
  abline(v=mean(strapp.all.Woodiness.mu.diff),col='purple',lty=2)
  legend('topright',lty=1,col=c('red','cyan','dark green','purple'),c('C.value','SeedMass','Life Cycle','Woodiness'),cex=.5,bty='n')
  dev.off()
  #write plot for div
  pdf(file='./output/plots/div_3trait_strapp_density.pdf',paper='a4')
  plot(density(strapp.all.seed.div.diff),ylim=c(0,10),xlim=c(-0.4,0.3),main=paste('SeedMass,C.value,Life Cycle (n=',length(all.Cvalue_data),')',sep=''),ylab='',xlab='absolute.diff in rho(obs-null)',col='cyan')
  lines(density(strapp.all.Cvalue.div.diff),col='red')
  lines(density(strapp.all.LifeCycle.div.diff),col='dark green')
  lines(density(strapp.all.Woodiness.div.diff),col='purple')
  abline(v=0,lty=2)
  abline(v=mean(strapp.all.seed.div.diff),col='cyan',lty=2)
  abline(v=mean(strapp.all.Cvalue.div.diff),col='red',lty=2)
  abline(v=mean(strapp.all.LifeCycle.div.diff),col='dark green',lty=2)
  abline(v=mean(strapp.all.Woodiness.div.diff),col='purple',lty=2)
  legend('topright',lty=1,col=c('red','cyan','dark green','purple'),c('C.value','SeedMass','Life Cycle','Woodiness'),cex=.5,bty='n')
  dev.off()
}

run_seedlifecycleCvaluewoodiness_seedrate_STRAPP<-function(treefile,eventfile,burnin,traitfile,seedratetable){
  tree<-read.tree(treefile)
  #get event data
  if (class(eventfile)!='bammdata'){
    Div_edata <- getEventData(tree, eventfile, burnin = burnin, nsamples = 1000,type='diversification')
  }else if (class(eventfile)=='bammdata'){
    Div_edata<-eventfile
  }
  #load traitdata
  trait<-read.table(traitfile,header=T,sep='\t',stringsAsFactors = F)
  #change lifecycle to categorical
  #lifecycle: perennial-0; annual-1
  trait[trait$cycle=='annual',]$cycle<-1
  trait[trait$cycle=='perennial',]$cycle<-0
  #read seedratetable and add to trait data
  seedratetable<-read.table(seedratetable,header=T,sep='\t',stringsAsFactors = F)
  seedratetable<-seedratetable[,c('names','mean.beta.rates')]
  colnames(seedratetable)[1]<-c('Species')
  trait<-merge(trait,seedratetable)
  #extracting subtreeBAMM
  subtree<-subtreeBAMM(Div_edata,tips=trait$Species)
  #all five analysis
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
  all.seedrate_data<-trait$mean.beta.rates
  names(all.seedrate_data)<-trait$Species
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
  #seedrate STRAPP
  strapp.all.seedrate.lambda<-traitDependentBAMM(subtree, traits = all.seedrate_data, 1000, rate = 'speciation', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.seedrate.mu<-traitDependentBAMM(subtree, traits = all.seedrate_data, 1000, rate = 'extinction', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.seedrate.div<-traitDependentBAMM(subtree, traits = all.seedrate_data, 1000, rate = 'net diversification', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.seedrate.lambda.diff<-abs(strapp.all.seedrate.lambda$obs.corr)-abs(strapp.all.seedrate.lambda$null)
  strapp.all.seedrate.mu.diff<-abs(strapp.all.seedrate.mu$obs.corr)-abs(strapp.all.seedrate.mu$null)
  strapp.all.seedrate.div.diff<-abs(strapp.all.seedrate.div$obs.corr)-abs(strapp.all.seedrate.div$null)
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
  #table for seedrate results
  seedrate.results<-data.frame(cbind(c(strapp.all.seedrate.lambda$estimate,strapp.all.seedrate.mu$estimate,strapp.all.seedrate.div$estimate),c(strapp.all.seedrate.lambda$p.value,strapp.all.seedrate.mu$p.value,strapp.all.seedrate.div$p.value)))
  row.names(seedrate.results)<-c('lambda','mu','div')
  colnames(seedrate.results)<-c('estimate','p-value')
  #seedcvaluelifecyclewoodiness.df<-cbind(seed.results,Cvalue.results,LifeCycle.results,Woodiness.results)
  seedcvaluelifecyclewoodiness.df<-cbind(seed.results,Cvalue.results,LifeCycle.results,Woodiness.results,seedrate.results)
  #colnames(seedcvaluelifecyclewoodiness.df)<-c('seedmass.estimate','seedmass.pvalue','Cvalue.estimate','Cvalue.pvalue','LifeCycle.estimate','LifeCycle.pvalue','Woodiness.estimate','Woodiness.pvalue')
  colnames(seedcvaluelifecyclewoodiness.df)<-c('seedmass.estimate','seedmass.pvalue','Cvalue.estimate','Cvalue.pvalue','LifeCycle.estimate','LifeCycle.pvalue','Woodiness.estimate','Woodiness.pvalue','seedrate.estimate','seedrate.pvalue')
  seedcvaluelifecyclewoodiness.df$parameter<-row.names(seedcvaluelifecyclewoodiness.df)
  row.names(seedcvaluelifecyclewoodiness.df)<-NULL
  seedcvaluelifecyclewoodiness.df<-seedcvaluelifecyclewoodiness.df[c(ncol(seedcvaluelifecyclewoodiness.df),1:(ncol(seedcvaluelifecyclewoodiness.df)-1))]
  write.table(seedcvaluelifecyclewoodiness.df,'./output/tables/seedCvalueLifeCycleSeedRate_STRAPP.txt',quote=F,row.names=F,sep='\t')
  #write plot for lambda
  pdf(file='./output/plots/lambda_4trait_strapp_density.pdf',paper='a4')
  plot(density(strapp.all.seed.lambda.diff),ylim=c(0,10),xlim=c(-0.4,0.65),main=paste('SeedMass,C.value,Life Cycle, SeedRate (n=',length(all.Cvalue_data),')',sep=''),ylab='',xlab='absolute.diff in rho(obs-null)',col='cyan')
  lines(density(strapp.all.Cvalue.lambda.diff),col='red')
  lines(density(strapp.all.LifeCycle.lambda.diff),col='dark green')
  lines(density(strapp.all.Woodiness.lambda.diff),col='purple')
  lines(density(strapp.all.seedrate.lambda.diff),col='dark blue')
  abline(v=0,lty=2)
  abline(v=mean(strapp.all.seed.lambda.diff),col='cyan',lty=2)
  abline(v=mean(strapp.all.Cvalue.lambda.diff),col='red',lty=2)
  abline(v=mean(strapp.all.LifeCycle.lambda.diff),col='dark green',lty=2)
  abline(v=mean(strapp.all.Woodiness.lambda.diff),col='purple',lty=2)
  abline(v=mean(strapp.all.seedrate.lambda.diff),col='dark blue',lty=2)
  legend('topright',lty=1,col=c('red','cyan','dark green','purple','dark blue'),c('C.value','SeedMass','Life Cycle','Woodiness', 'SeedRate'),cex=.5,bty='n')
  dev.off()
  #write plot for mu
  pdf(file='./output/plots/mu_4trait_strapp_density.pdf',paper='a4')
  plot(density(strapp.all.seed.mu.diff),ylim=c(0,10),xlim=c(-0.4,0.65),main=paste('SeedMass,C.value,Life Cycle,SeedRate (n=',length(all.Cvalue_data),')',sep=''),ylab='',xlab='absolute.diff in rho(obs-null)',col='cyan')
  lines(density(strapp.all.Cvalue.mu.diff),col='red')
  lines(density(strapp.all.LifeCycle.mu.diff),col='dark green')
  lines(density(strapp.all.Woodiness.mu.diff),col='purple')
  lines(density(strapp.all.seedrate.mu.diff),col='dark blue')
  abline(v=0,lty=2)
  abline(v=mean(strapp.all.seed.mu.diff),col='cyan',lty=2)
  abline(v=mean(strapp.all.Cvalue.mu.diff),col='red',lty=2)
  abline(v=mean(strapp.all.LifeCycle.mu.diff),col='dark green',lty=2)
  abline(v=mean(strapp.all.Woodiness.mu.diff),col='purple',lty=2)
  abline(v=mean(strapp.all.seedrate.mu.diff),col='dark blue',lty=2)
  legend('topright',lty=1,col=c('red','cyan','dark green','purple','dark blue'),c('C.value','SeedMass','Life Cycle','Woodiness','SeedRate'),cex=.5,bty='n')
  dev.off()
  #write plot for div
  pdf(file='./output/plots/div_4trait_strapp_density.pdf',paper='a4')
  plot(density(strapp.all.seed.div.diff),ylim=c(0,10),xlim=c(-0.4,0.65),main=paste('SeedMass,C.value,Life Cycle,SeedRate (n=',length(all.Cvalue_data),')',sep=''),ylab='',xlab='absolute.diff in rho(obs-null)',col='cyan')
  lines(density(strapp.all.Cvalue.div.diff),col='red')
  lines(density(strapp.all.LifeCycle.div.diff),col='dark green')
  lines(density(strapp.all.Woodiness.div.diff),col='purple')
  lines(density(strapp.all.seedrate.div.diff),col='dark blue')
  abline(v=0,lty=2)
  abline(v=mean(strapp.all.seed.div.diff),col='cyan',lty=2)
  abline(v=mean(strapp.all.Cvalue.div.diff),col='red',lty=2)
  abline(v=mean(strapp.all.LifeCycle.div.diff),col='dark green',lty=2)
  abline(v=mean(strapp.all.Woodiness.div.diff),col='purple',lty=2)
  abline(v=mean(strapp.all.seedrate.div.diff),col='dark blue',lty=2)
  legend('topright',lty=1,col=c('red','cyan','dark green','purple','dark blue'),c('C.value','SeedMass','Life Cycle','Woodiness','SeedRate'),cex=.5,bty='n')
  dev.off()
}

run_seedlifecycleCvaluewoodinessheight_phenorate_STRAPP<-function(treefile,eventfile,burnin,traitfile,seedratetable,Cvalueratetable,heightratetable){
  tree<-read.tree(treefile)
  #get event data
  if (class(eventfile)!='bammdata'){
    Div_edata <- getEventData(tree, eventfile, burnin = burnin, nsamples = 1000,type='diversification')
  }else if (class(eventfile)=='bammdata'){
    Div_edata<-eventfile
  }
  #load traitdata
  trait<-read.table(traitfile,header=T,sep='\t',stringsAsFactors = F)
  #change lifecycle to categorical
  #lifecycle: perennial-0; annual-1
  trait[trait$cycle=='annual',]$cycle<-1
  trait[trait$cycle=='perennial',]$cycle<-0
  #read seedratetable and add to trait data
  seedratetable<-read.table(seedratetable,header=T,sep='\t',stringsAsFactors = F)
  seedratetable<-seedratetable[,c('names','mean.beta.rates')]
  colnames(seedratetable)[1]<-c('Species')
  trait<-merge(trait,seedratetable)
  #read cvalueratetable and add to trait data
  Cvalueratetable<-read.table(Cvalueratetable,header=T,sep='\t',stringsAsFactors = F)
  Cvalueratetable<-Cvalueratetable[,c('names','mean.beta.rates')]
  colnames(Cvalueratetable)<-c('Species','mean.Cvalue.beta.rates')
  trait<-merge(trait,Cvalueratetable)
  #read cvalueratetable and add to trait data
  heightratetable<-read.table(heightratetable,header=T,sep='\t',stringsAsFactors = F)
  heightratetable<-heightratetable[,c('names','mean.beta.rates')]
  colnames(heightratetable)<-c('Species','mean.height.beta.rates')
  trait<-merge(trait,heightratetable)
  write.table(trait,file='./output/tables/seedlifecycleCvaluewoodinessheight_phenorates.txt',quote=F,sep='\t',row.names=F)
  #extracting subtreeBAMM
  subtree<-subtreeBAMM(Div_edata,tips=trait$Species)
  #all five analysis
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
  all.Height_data<-trait$log10.height
  names(all.Height_data)<-trait$Species
  
  all.seedrate_data<-trait$mean.beta.rates
  names(all.seedrate_data)<-trait$Species
  all.Cvaluerate_data<-trait$mean.Cvalue.beta.rates
  names(all.Cvaluerate_data)<-trait$Species
  all.heightrate_data<-trait$mean.height.beta.rates
  names(all.heightrate_data)<-trait$Species
  
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
  #Height STRAPP
  strapp.all.Height.lambda<-traitDependentBAMM(subtree, traits = all.Height_data, 1000, rate = 'speciation', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.Height.mu<-traitDependentBAMM(subtree, traits = all.Height_data, 1000, rate = 'extinction', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.Height.div<-traitDependentBAMM(subtree, traits = all.Height_data, 1000, rate = 'net diversification', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.Height.lambda.diff<-abs(strapp.all.Height.lambda$obs.corr)-abs(strapp.all.Height.lambda$null)
  strapp.all.Height.mu.diff<-abs(strapp.all.Height.mu$obs.corr)-abs(strapp.all.Height.mu$null)
  strapp.all.Height.div.diff<-abs(strapp.all.Height.div$obs.corr)-abs(strapp.all.Height.div$null)
  #seedrate STRAPP
  strapp.all.seedrate.lambda<-traitDependentBAMM(subtree, traits = all.seedrate_data, 1000, rate = 'speciation', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.seedrate.mu<-traitDependentBAMM(subtree, traits = all.seedrate_data, 1000, rate = 'extinction', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.seedrate.div<-traitDependentBAMM(subtree, traits = all.seedrate_data, 1000, rate = 'net diversification', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.seedrate.lambda.diff<-abs(strapp.all.seedrate.lambda$obs.corr)-abs(strapp.all.seedrate.lambda$null)
  strapp.all.seedrate.mu.diff<-abs(strapp.all.seedrate.mu$obs.corr)-abs(strapp.all.seedrate.mu$null)
  strapp.all.seedrate.div.diff<-abs(strapp.all.seedrate.div$obs.corr)-abs(strapp.all.seedrate.div$null)
  #Cvaluerate STRAPP
  strapp.all.Cvaluerate.lambda<-traitDependentBAMM(subtree, traits = all.Cvaluerate_data, 1000, rate = 'speciation', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.Cvaluerate.mu<-traitDependentBAMM(subtree, traits = all.Cvaluerate_data, 1000, rate = 'extinction', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.Cvaluerate.div<-traitDependentBAMM(subtree, traits = all.Cvaluerate_data, 1000, rate = 'net diversification', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.Cvaluerate.lambda.diff<-abs(strapp.all.Cvaluerate.lambda$obs.corr)-abs(strapp.all.Cvaluerate.lambda$null)
  strapp.all.Cvaluerate.mu.diff<-abs(strapp.all.Cvaluerate.mu$obs.corr)-abs(strapp.all.Cvaluerate.mu$null)
  strapp.all.Cvaluerate.div.diff<-abs(strapp.all.Cvaluerate.div$obs.corr)-abs(strapp.all.Cvaluerate.div$null)
  #heightrate STRAPP
  strapp.all.heightrate.lambda<-traitDependentBAMM(subtree, traits = all.heightrate_data, 1000, rate = 'speciation', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.heightrate.mu<-traitDependentBAMM(subtree, traits = all.heightrate_data, 1000, rate = 'extinction', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.heightrate.div<-traitDependentBAMM(subtree, traits = all.heightrate_data, 1000, rate = 'net diversification', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.heightrate.lambda.diff<-abs(strapp.all.heightrate.lambda$obs.corr)-abs(strapp.all.heightrate.lambda$null)
  strapp.all.heightrate.mu.diff<-abs(strapp.all.heightrate.mu$obs.corr)-abs(strapp.all.heightrate.mu$null)
  strapp.all.heightrate.div.diff<-abs(strapp.all.heightrate.div$obs.corr)-abs(strapp.all.heightrate.div$null)
  
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
  #table for Height results
  Height.results<-data.frame(cbind(c(strapp.all.Height.lambda$estimate,strapp.all.Height.mu$estimate,strapp.all.Height.div$estimate),c(strapp.all.Height.lambda$p.value,strapp.all.Height.mu$p.value,strapp.all.Height.div$p.value)))
  row.names(Height.results)<-c('lambda','mu','div')
  colnames(Height.results)<-c('estimate','p-value')
  #table for seedrate results
  seedrate.results<-data.frame(cbind(c(strapp.all.seedrate.lambda$estimate,strapp.all.seedrate.mu$estimate,strapp.all.seedrate.div$estimate),c(strapp.all.seedrate.lambda$p.value,strapp.all.seedrate.mu$p.value,strapp.all.seedrate.div$p.value)))
  row.names(seedrate.results)<-c('lambda','mu','div')
  colnames(seedrate.results)<-c('estimate','p-value')
  #table for Cvaluerate results
  Cvaluerate.results<-data.frame(cbind(c(strapp.all.Cvaluerate.lambda$estimate,strapp.all.Cvaluerate.mu$estimate,strapp.all.Cvaluerate.div$estimate),c(strapp.all.Cvaluerate.lambda$p.value,strapp.all.Cvaluerate.mu$p.value,strapp.all.Cvaluerate.div$p.value)))
  row.names(Cvaluerate.results)<-c('lambda','mu','div')
  colnames(Cvaluerate.results)<-c('estimate','p-value')
  #table for height results
  heightrate.results<-data.frame(cbind(c(strapp.all.heightrate.lambda$estimate,strapp.all.heightrate.mu$estimate,strapp.all.heightrate.div$estimate),c(strapp.all.heightrate.lambda$p.value,strapp.all.heightrate.mu$p.value,strapp.all.heightrate.div$p.value)))
  row.names(heightrate.results)<-c('lambda','mu','div')
  colnames(heightrate.results)<-c('estimate','p-value')
  #seedcvaluelifecyclewoodiness.df<-cbind(seed.results,Cvalue.results,LifeCycle.results,Woodiness.results)
  seedcvaluelifecyclewoodinessheight.df<-cbind(seed.results,Cvalue.results,LifeCycle.results,Woodiness.results,Height.results,seedrate.results,Cvaluerate.results,heightrate.results)
  #colnames(seedcvaluelifecyclewoodiness.df)<-c('seedmass.estimate','seedmass.pvalue','Cvalue.estimate','Cvalue.pvalue','LifeCycle.estimate','LifeCycle.pvalue','Woodiness.estimate','Woodiness.pvalue')
  colnames(seedcvaluelifecyclewoodinessheight.df)<-c('seedmass.estimate','seedmass.pvalue','Cvalue.estimate','Cvalue.pvalue','LifeCycle.estimate','LifeCycle.pvalue','Woodiness.estimate','Woodiness.pvalue','height.estimate','height.pvalue','seedrate.estimate','seedrate.pvalue','Cvaluerate.estimate','Cvaluerate.pvalue','heightrate.estimate','heightrate.pvalue')
  seedcvaluelifecyclewoodinessheight.df$parameter<-row.names(seedcvaluelifecyclewoodinessheight.df)
  row.names(seedcvaluelifecyclewoodinessheight.df)<-NULL
  seedcvaluelifecyclewoodinessheight.df<-seedcvaluelifecyclewoodinessheight.df[c(ncol(seedcvaluelifecyclewoodinessheight.df),1:(ncol(seedcvaluelifecyclewoodinessheight.df)-1))]
  write.table(seedcvaluelifecyclewoodinessheight.df,'./output/tables/seedCvalueLifeCycleSeedRate_STRAPP.txt',quote=F,row.names=F,sep='\t')
  #write plot for lambda
  colours<-rainbow(n=8)
  pdf(file='./output/plots/lambda_alltraitrates_strapp_density.pdf',paper='a4')
  plot(density(strapp.all.seed.lambda.diff),ylim=c(0,10),xlim=c(-0.4,0.65),main=paste('SeedMass,C.value,Life Cycle, SeedRate (n=',length(all.Cvalue_data),')',sep=''),ylab='',xlab='absolute.diff in rho(obs-null)',col=colours[1])
  lines(density(strapp.all.Cvalue.lambda.diff),col=colours[2])
  lines(density(strapp.all.LifeCycle.lambda.diff),col=colours[3])
  lines(density(strapp.all.Woodiness.lambda.diff),col=colours[4])
  lines(density(strapp.all.Height.lambda.diff),col=colours[5])
  lines(density(strapp.all.seedrate.lambda.diff),col=colours[6])
  lines(density(strapp.all.Cvaluerate.lambda.diff),col=colours[7])
  lines(density(strapp.all.heightrate.lambda.diff),col=colours[8])
  
  abline(v=0,lty=2)
  abline(v=mean(strapp.all.seed.lambda.diff),col=colours[1],lty=2)
  abline(v=mean(strapp.all.Cvalue.lambda.diff),col=colours[2],lty=2)
  abline(v=mean(strapp.all.LifeCycle.lambda.diff),col=colours[3],lty=2)
  abline(v=mean(strapp.all.Woodiness.lambda.diff),col=colours[4],lty=2)
  abline(v=mean(strapp.all.Height.lambda.diff),col=colours[5],lty=2)
  abline(v=mean(strapp.all.seedrate.lambda.diff),col=colours[6],lty=2)
  abline(v=mean(strapp.all.Cvaluerate.lambda.diff),col=colours[7],lty=2)
  abline(v=mean(strapp.all.heightrate.lambda.diff),col=colours[8],lty=2)
  
  legend('topright',lty=1,col=colours,c('seed','C.value','Life Cycle','Woodiness','Height','SeedRate','CvalueRate','HeightRate'),cex=.5,bty='n')
  dev.off()
  #write plot for mu
  colours<-rainbow(n=8)
  pdf(file='./output/plots/mu_alltraitrates_strapp_density.pdf',paper='a4')
  plot(density(strapp.all.seed.mu.diff),ylim=c(0,10),xlim=c(-0.4,0.65),main=paste('SeedMass,C.value,Life Cycle, SeedRate (n=',length(all.Cvalue_data),')',sep=''),ylab='',xlab='absolute.diff in rho(obs-null)',col=colours[1])
  lines(density(strapp.all.Cvalue.mu.diff),col=colours[2])
  lines(density(strapp.all.LifeCycle.mu.diff),col=colours[3])
  lines(density(strapp.all.Woodiness.mu.diff),col=colours[4])
  lines(density(strapp.all.Height.mu.diff),col=colours[5])
  lines(density(strapp.all.seedrate.mu.diff),col=colours[6])
  lines(density(strapp.all.Cvaluerate.mu.diff),col=colours[7])
  lines(density(strapp.all.heightrate.mu.diff),col=colours[8])
  
  abline(v=0,lty=2)
  abline(v=mean(strapp.all.seed.mu.diff),col=colours[1],lty=2)
  abline(v=mean(strapp.all.Cvalue.mu.diff),col=colours[2],lty=2)
  abline(v=mean(strapp.all.LifeCycle.mu.diff),col=colours[3],lty=2)
  abline(v=mean(strapp.all.Woodiness.mu.diff),col=colours[4],lty=2)
  abline(v=mean(strapp.all.Height.mu.diff),col=colours[5],lty=2)
  abline(v=mean(strapp.all.seedrate.mu.diff),col=colours[6],lty=2)
  abline(v=mean(strapp.all.Cvaluerate.mu.diff),col=colours[7],lty=2)
  abline(v=mean(strapp.all.heightrate.mu.diff),col=colours[8],lty=2)
  
  legend('topright',lty=1,col=colours,c('seed','C.value','Life Cycle','Woodiness','Height','SeedRate','CvalueRate','HeightRate'),cex=.5,bty='n')
  dev.off()
  #write plot for div
  pdf(file='./output/plots/div_alltraitrates_strapp_density.pdf',paper='a4')
  plot(density(strapp.all.seed.div.diff),ylim=c(0,10),xlim=c(-0.4,0.65),main=paste('SeedMass,C.value,Life Cycle, SeedRate (n=',length(all.Cvalue_data),')',sep=''),ylab='',xlab='absolute.diff in rho(obs-null)',col=colours[1])
  lines(density(strapp.all.Cvalue.div.diff),col=colours[2])
  lines(density(strapp.all.LifeCycle.div.diff),col=colours[3])
  lines(density(strapp.all.Woodiness.div.diff),col=colours[4])
  lines(density(strapp.all.Height.div.diff),col=colours[5])
  lines(density(strapp.all.seedrate.div.diff),col=colours[6])
  lines(density(strapp.all.Cvaluerate.div.diff),col=colours[7])
  lines(density(strapp.all.heightrate.div.diff),col=colours[8])
  
  abline(v=0,lty=2)
  abline(v=mean(strapp.all.seed.div.diff),col=colours[1],lty=2)
  abline(v=mean(strapp.all.Cvalue.div.diff),col=colours[2],lty=2)
  abline(v=mean(strapp.all.LifeCycle.div.diff),col=colours[3],lty=2)
  abline(v=mean(strapp.all.Woodiness.div.diff),col=colours[4],lty=2)
  abline(v=mean(strapp.all.Height.div.diff),col=colours[5],lty=2)
  abline(v=mean(strapp.all.seedrate.div.diff),col=colours[6],lty=2)
  abline(v=mean(strapp.all.Cvaluerate.div.diff),col=colours[7],lty=2)
  abline(v=mean(strapp.all.heightrate.div.diff),col=colours[8],lty=2)
  
  legend('topright',lty=1,col=colours,c('seed','C.value','Life Cycle','Woodiness','Height','SeedRate','CvalueRate','HeightRate'),cex=.5,bty='n')
  dev.off()
}

######################################################
#function to collate tree + seed + lifecycle + Cvalue + woodiness datasets
#load the KewSID + Qian tree dataset
prepare_lifecycle_Cvalue_woodiness_dataset<-function(){
  kewsid.TPL<-read.table('./output/tables/BAMM_Species_Seeddata.txt',sep='\t',header=F)
  colnames(kewsid.TPL)<-c('Merged','log10.SeedWeight')
  tree<-read.tree('./output/trees/BAMM_Species_tree.tree')
  #############LIFE CYCLE
  lifecycle.TPL<-read.table('./output/tables/lifecycle_dataset_TPL.txt',header=T,sep='\t')
  annuals<-read.csv('./raw_data/annual_CvalueDB.txt', header = T, sep = '\t')
  perennials<-read.csv('./raw_data/perennial_CvalueDB.txt',header = T, sep = '\t')
  #annuals
  annuals$Merged<-paste(annuals$Genus, annuals$Species, sep = '_')
  perennials$Merged<-paste(perennials$Genus, perennials$Species, sep = '_')
  lifecycle.TPL$Old.Species<-paste(lifecycle.TPL$Genus,lifecycle.TPL$Species,sep='_')
  annuals<-annuals[,c(6,3)]
  annuals.TPL<-merge(annuals,lifecycle.TPL,by.x='Merged',by.y='Old.Species')
  colnames(annuals.TPL)[c(1,2)]<-c('Old.Species','C.value')
  annuals.TPL$Merged<-paste(annuals.TPL$New.Genus,annuals.TPL$New.Species,sep='_')
  annuals.TPL<-annuals.TPL[!duplicated(annuals.TPL$Merged),]
  annuals.TPL<-annuals.TPL[,c('Merged','C.value')]
  annuals.TPL$cycle<-'annual'
  #perennials
  perennials$Merged<-paste(perennials$Genus, perennials$Species, sep = '_')
  perennials$Merged<-paste(perennials$Genus, perennials$Species, sep = '_')
  lifecycle.TPL$Old.Species<-paste(lifecycle.TPL$Genus,lifecycle.TPL$Species,sep='_')
  perennials<-perennials[,c(6,3)]
  perennials.TPL<-merge(perennials,lifecycle.TPL,by.x='Merged',by.y='Old.Species')
  colnames(perennials.TPL)[c(1,2)]<-c('Old.Species','C.value')
  perennials.TPL$Merged<-paste(perennials.TPL$New.Genus,perennials.TPL$New.Species,sep='_')
  perennials.TPL<-perennials.TPL[!duplicated(perennials.TPL$Merged),]
  perennials.TPL<-perennials.TPL[,c('Merged','C.value')]
  perennials.TPL$cycle<-'perennial'
  lifecycle.TPL<-rbind(annuals.TPL,perennials.TPL)
  #sort out duplicates
  lifecycle.TPL<-lifecycle.TPL[!duplicated(lifecycle.TPL$Merged),]
  lifecycle.TPL<-unique(lifecycle.TPL)
  
  #############WOODINESS
  woodiness<-read.csv('./raw_data/GlobalWoodinessDatabase.csv',stringsAsFactors = F)
  #just using the Zanne et al coding here
  woodiness$woody.state<-NA
  woodiness[woodiness$woodiness=='H',]$woody.state<-0
  woodiness[woodiness$woodiness=='W',]$woody.state<-1
  #dropping the "variable species"
  woodiness<-na.omit(woodiness)
  woodiness$Merged<-sub(woodiness$gs,pattern=' ',replacement='_')
  woodiness<-woodiness[,c("Merged","woody.state")]
  woodiness.TPL<-read.table('./output/tables/woodiness_dataset_TPL.txt',header=T,sep='\t')
  woodiness.TPL$Old.Species<-paste(woodiness.TPL$Genus,woodiness.TPL$Species,sep='_')
  woodiness.TPL<-merge(woodiness,woodiness.TPL,by.x='Merged',by.y='Old.Species')
  colnames(woodiness.TPL)[c(1,2)]<-c('Old.Species','woodiness.state')
  woodiness.TPL$Merged<-paste(woodiness.TPL$New.Genus,woodiness.TPL$New.Species,sep='_')
  woodiness.TPL<-woodiness.TPL[!duplicated(woodiness.TPL$Merged),]
  woodiness.TPL<-woodiness.TPL[,c('Merged','woodiness.state')]
  
  
  
  #merge all datasets together
  kewsid.TPL.lifecycle<-merge(kewsid.TPL,lifecycle.TPL)
  kewsid.TPL.lifecycle.woodiness<-merge(kewsid.TPL.lifecycle,woodiness.TPL)
  colnames(kewsid.TPL.lifecycle.woodiness)[1]<-c('Species')
  write.table(kewsid.TPL.lifecycle.woodiness,'./output/tables/kewsid_lifecycle_woodiness_TPL.txt',sep='\t',quote=F,row.names=F)
  
}

prepare_lifecycle_Cvalue_woodiness_height_dataset<-function(){
  kewsid.TPL<-read.table('./output/tables/BAMM_Species_Seeddata.txt',sep='\t',header=F)
  colnames(kewsid.TPL)<-c('Merged','log10.SeedWeight')
  tree<-read.tree('./output/trees/BAMM_Species_tree.tree')
  #############LIFE CYCLE
  lifecycle.TPL<-read.table('./output/tables/lifecycle_dataset_TPL.txt',header=T,sep='\t')
  annuals<-read.csv('./raw_data/annual_CvalueDB.txt', header = T, sep = '\t')
  perennials<-read.csv('./raw_data/perennial_CvalueDB.txt',header = T, sep = '\t')
  #annuals
  annuals$Merged<-paste(annuals$Genus, annuals$Species, sep = '_')
  perennials$Merged<-paste(perennials$Genus, perennials$Species, sep = '_')
  lifecycle.TPL$Old.Species<-paste(lifecycle.TPL$Genus,lifecycle.TPL$Species,sep='_')
  annuals<-annuals[,c(6,3)]
  annuals.TPL<-merge(annuals,lifecycle.TPL,by.x='Merged',by.y='Old.Species')
  colnames(annuals.TPL)[c(1,2)]<-c('Old.Species','C.value')
  annuals.TPL$Merged<-paste(annuals.TPL$New.Genus,annuals.TPL$New.Species,sep='_')
  annuals.TPL<-annuals.TPL[!duplicated(annuals.TPL$Merged),]
  annuals.TPL<-annuals.TPL[,c('Merged','C.value')]
  annuals.TPL$cycle<-'annual'
  #perennials
  perennials$Merged<-paste(perennials$Genus, perennials$Species, sep = '_')
  perennials$Merged<-paste(perennials$Genus, perennials$Species, sep = '_')
  lifecycle.TPL$Old.Species<-paste(lifecycle.TPL$Genus,lifecycle.TPL$Species,sep='_')
  perennials<-perennials[,c(6,3)]
  perennials.TPL<-merge(perennials,lifecycle.TPL,by.x='Merged',by.y='Old.Species')
  colnames(perennials.TPL)[c(1,2)]<-c('Old.Species','C.value')
  perennials.TPL$Merged<-paste(perennials.TPL$New.Genus,perennials.TPL$New.Species,sep='_')
  perennials.TPL<-perennials.TPL[!duplicated(perennials.TPL$Merged),]
  perennials.TPL<-perennials.TPL[,c('Merged','C.value')]
  perennials.TPL$cycle<-'perennial'
  lifecycle.TPL<-rbind(annuals.TPL,perennials.TPL)
  #sort out duplicates
  lifecycle.TPL<-lifecycle.TPL[!duplicated(lifecycle.TPL$Merged),]
  lifecycle.TPL<-unique(lifecycle.TPL)
  
  #############WOODINESS
  woodiness<-read.csv('./raw_data/GlobalWoodinessDatabase.csv',stringsAsFactors = F)
  #just using the Zanne et al coding here
  woodiness$woody.state<-NA
  woodiness[woodiness$woodiness=='H',]$woody.state<-0
  woodiness[woodiness$woodiness=='W',]$woody.state<-1
  #dropping the "variable species"
  woodiness<-na.omit(woodiness)
  woodiness$Merged<-sub(woodiness$gs,pattern=' ',replacement='_')
  woodiness<-woodiness[,c("Merged","woody.state")]
  woodiness.TPL<-read.table('./output/tables/woodiness_dataset_TPL.txt',header=T,sep='\t')
  woodiness.TPL$Old.Species<-paste(woodiness.TPL$Genus,woodiness.TPL$Species,sep='_')
  woodiness.TPL<-merge(woodiness,woodiness.TPL,by.x='Merged',by.y='Old.Species')
  colnames(woodiness.TPL)[c(1,2)]<-c('Old.Species','woodiness.state')
  woodiness.TPL$Merged<-paste(woodiness.TPL$New.Genus,woodiness.TPL$New.Species,sep='_')
  woodiness.TPL<-woodiness.TPL[!duplicated(woodiness.TPL$Merged),]
  woodiness.TPL<-woodiness.TPL[,c('Merged','woodiness.state')]
  #height
  height.TPL<-read.table(file='./output/tables/BAMM_Species_tree_height_table.txt',sep='\t',header=F)
  colnames(height.TPL)<-c('Merged','log10.height')
  #merge all datasets together
  kewsid.TPL.lifecycle<-merge(kewsid.TPL,lifecycle.TPL)
  kewsid.TPL.lifecycle.woodiness<-merge(kewsid.TPL.lifecycle,woodiness.TPL)
  kewsid.TPL.lifecycle.woodiness.height<-merge(kewsid.TPL.lifecycle.woodiness,height.TPL)
  colnames(kewsid.TPL.lifecycle.woodiness.height)[1]<-c('Species')
  write.table(kewsid.TPL.lifecycle.woodiness.height,'./output/tables/kewsid_lifecycle_woodiness_height_TPL.txt',sep='\t',quote=F,row.names=F)
  
}
####this runs a phylogenetic anova of differences in seed weight of annuals vs perennials
lifecycle_seed_phylanova<-function(treefile,traitfile){
  tree<-read.tree(treefile)
  #load traitdata
  trait<-read.table(traitfile,header=T,sep='\t',stringsAsFactors = F)
  tree.trait<-drop.tip(tree,setdiff(tree$tip.label,trait$Species))
  cycle<-trait$cycle
  names(cycle)<-trait$Species
  seed<-trait$log10.SeedWeight
  names(seed)<-trait$Species
  phylANOVAseedcycle<-phylANOVA(tree.trait,cycle,seed)
  return(phylANOVAseedcycle)
}

####this runs a phylogenetic anova of differences in seed weight of woody vs herbaceous
woody_seed_phylanova<-function(treefile,traitfile){
  tree<-read.tree(treefile)
  #load traitdata
  trait<-read.table(traitfile,header=T,sep='\t',stringsAsFactors = F)
  tree.trait<-drop.tip(tree,setdiff(tree$tip.label,trait$Species))
  woody<-trait$woodiness.state
  names(woody)<-trait$Species
  seed<-trait$log10.SeedWeight
  names(seed)<-trait$Species
  phylANOVAseedwoody<-phylANOVA(tree.trait,woody,seed)
  return(phylANOVAseedwoody)
}

annual_vs_perennial<-function(traitfile){
  #load the trait dataset
  trait<-read.table(traitfile,header=T,sep='\t',stringsAsFactors = F)
  ######plot seed mass in annuals and perennials
  trait$cycle<-factor(trait$cycle)
  pdf(file='./output/plots/annual_vs_perennial_boxplot.pdf',paper='a4')
  boxplot(log10.SeedWeight~as.factor(cycle),data=trait,horizontal=T, las=1,boxwex=0.3,cex=0.3,col=c('light grey','white'), main='Life Cycle and Seed Mass',xaxt='n',yaxt='n',xlab='Species Seed Weight',notch=T)
  legend(x='topright', fill = c('light grey', 'white'), legend = c('annual','perennial'), bty = 'n',cex=.7)
  ticks <- seq(-6, 3, by=1)
  labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
  axis(1, at=ticks, labels=labels,cex.axis=.7)
  dev.off()
}
