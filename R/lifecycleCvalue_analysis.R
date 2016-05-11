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


##########build the Cvalue dataset, write output to KewCvaluedata_GenusLevel_taxonlookup.txt ( to be used later)
#tablefile: 'genome_size_allKew.txt'
build_Cvalue_dataset<-function(tablefile){
  Cvalue<-read.table('./raw_data/genome_size_allKew.txt', sep="\t", header = T)
  Cvalue$Merged<-paste(Cvalue$Genus, Cvalue$Species, sep = '_')
  Cvalue<-Cvalue[!duplicated(Cvalue$Merged),]
  colnames(Cvalue)[3]<-'C.value'
  Genera.Cvalue<-aggregate(C.value~Genus, data=Cvalue, mean)
  colnames(Genera.Cvalue)<-c('genus','C.value')
  Cvalue<-Genera.Cvalue
  write.table(Cvalue, file = './output/tables/KewCvalue_GenusLevel.txt', quote = FALSE, sep = '\t', row.names = FALSE)
}  


build_lifecycle_dataset<-function(annualsfile,perennialsfile){
  annuals<-read.csv(annualsfile, header = T, sep = '\t')
  perennials<-read.csv(perennialsfile,header = T, sep = '\t')
  annuals$Merged<-paste(annuals$Genus, annuals$Species, sep = '_')
  perennials$Merged<-paste(perennials$Genus, perennials$Species, sep = '_')
  annuals<-annuals[!duplicated(annuals$Merged),]
  annuals$cycle<-'annual'
  perennials<-perennials[!duplicated(perennials$Merged),]
  perennials$cycle<-'perennial'
  all<-rbind(annuals,perennials)
  all$Original.Reference<-NULL
  all$Paper<-NULL
  colnames(all)[3]<-'C.value'
  all<-all[!(all$Species==""),]
  #adding "annuality" info (calculating proportion of annual species for each genera)
  cyclecount<-count(all, c("Genus", "cycle"))
  cyclecount<-cyclecount[order(cyclecount$Genus),]
  cyclecount$Genus<-as.character(cyclecount$Genus)
  annuality<-function(x){
    annual.count<-subset(cyclecount,cyclecount$Genus==x & cyclecount$cycle=='annual')
    perennial.count<-subset(cyclecount,cyclecount$Genus==x & cyclecount$cycle=='perennial')
    annual.count<-annual.count$freq
    perennial.count<-perennial.count$freq
    if (length(annual.count)==0){annual.count<-0}
    if (length(perennial.count)==0){perennial.count<-0}
    annuality<-annual.count/(annual.count+perennial.count)
    res<-data.frame(x,annuality)
    return(res)
  }
  list<-lapply(unique(cyclecount$Genus),function(x) annuality(x))
  annual <- data.frame(matrix(unlist(list), nrow=length(list), byrow=T))
  annual$X1<-unique(cyclecount$Genus)
  colnames(annual)<-c('Genus','Annuality')
  write.table(annual,file='./output/tables/KewLifeCycle_GenusLevel.txt',sep='\t',row.names=F,quote=F)
}
annual_vs_perennial<-function(treefile,seedtablefile,cvaluetablefile,lifecycletablefile){
  tree<-read.tree(treefile)
  #load the seed mass dataset (n=4105, run through BAMM), NO HEADER
  seed<-read.table(seedtablefile,header=F,sep='\t')
  colnames(seed)<-c('genus','log10.Seed.Weight')
  #load the Cvalue dataset
  Cvalue<-read.table(cvaluetablefile,header=T,sep='\t')
  colnames(Cvalue)<-c('genus','C.value')
  annual<-read.table(lifecycletablefile,header=T,sep='\t')
  #merge seedmass and Cvalue
  seedCvalue<-merge(seed,Cvalue,by.x='genus',by.y='genus')
  #merge seedmass and lifecycle
  seedlifecycle<-merge(seed,annual,by.x='genus',by.y='Genus')
  #merge all three
  seedCvaluelifecycle<-merge(seedCvalue,seedlifecycle,by.x='genus',by.y='genus')
  seedCvaluelifecycle.strict<-seedCvaluelifecycle
  seedCvaluelifecycle.strict<-subset(seedCvaluelifecycle.strict,seedCvaluelifecycle.strict$Annuality==1 | seedCvaluelifecycle.strict$Annuality==0)
  colnames(seedCvaluelifecycle.strict)[2]<-'log10.Seed.Weight'
  seedCvaluelifecycle.strict$log10.Seed.Weight.y<-NULL
  #run pgls
  tree$node.label<-NULL
  comp.data<-comparative.data(tree,seedCvaluelifecycle.strict,names.col = 'genus')
  pgls.Seed.Annuality<-pgls(log10.Seed.Weight~Annuality, data=comp.data, lambda="ML")
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
run_seedlifecycleCvalue_STRAPP<-function(treefile,eventfile,burnin,seedtablefile,cvaluetablefile,lifecycletablefile){
  tree<-read.tree(treefile)
  #get event data
  print('Analysing eventfile')
  Div_edata <- getEventData(tree, eventfile, burnin = burnin, nsamples = 5000,type='diversification')
  #load the seed mass dataset (n=4105, run through BAMM), NO HEADER
  seed<-read.table(seedtablefile,header=F,sep='\t')
  colnames(seed)<-c('genus','log10.Seed.Weight')
  #load the Cvalue dataset
  Cvalue<-read.table(cvaluetablefile,header=T,sep='\t')
  colnames(Cvalue)<-c('genus','C.value')
  annual<-read.table(lifecycletablefile,header=T,sep='\t')
  #merge seedmass and Cvalue
  seedCvalue<-merge(seed,Cvalue,by.x='genus',by.y='genus')
  #merge seedmass and lifecycle
  seedlifecycle<-merge(seed,annual,by.x='genus',by.y='Genus')
  #merge all three
  seedCvaluelifecycle<-merge(seedCvalue,seedlifecycle,by.x='genus',by.y='genus')
  seedCvaluelifecycle.strict<-seedCvaluelifecycle
  seedCvaluelifecycle.strict<-subset(seedCvaluelifecycle.strict,seedCvaluelifecycle.strict$Annuality==1 | seedCvaluelifecycle.strict$Annuality==0)
  colnames(seedCvaluelifecycle.strict)[2]<-'log10.Seed.Weight'
  seedCvaluelifecycle.strict$log10.Seed.Weight.y<-NULL
  #all three analysis
  all.Cvalue_data<-seedCvaluelifecycle.strict$C.value
  names(all.Cvalue_data)<-seedCvaluelifecycle.strict$genus
  all.seed_data<-seedCvaluelifecycle.strict$log10.Seed.Weight
  names(all.seed_data)<-seedCvaluelifecycle.strict$genus
  all.LifeCycle_data<-seedCvaluelifecycle.strict$Annuality
  names(all.LifeCycle_data)<-seedCvaluelifecycle.strict$genus
  #running STRAPP for all dataset
  #seed mass STRAPP
  strapp.all.seed.lambda<-traitDependentBAMM(Div_edata, traits = all.seed_data, 1000, rate = 'speciation', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.seed.mu<-traitDependentBAMM(Div_edata, traits = all.seed_data, 1000, rate = 'extinction', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.seed.div<-traitDependentBAMM(Div_edata, traits = all.seed_data, 1000, rate = 'net diversification', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.seed.lambda.diff<-abs(strapp.all.seed.lambda$obs.corr)-abs(strapp.all.seed.lambda$null)
  strapp.all.seed.mu.diff<-abs(strapp.all.seed.mu$obs.corr)-abs(strapp.all.seed.mu$null)
  strapp.all.seed.div.diff<-abs(strapp.all.seed.div$obs.corr)-abs(strapp.all.seed.div$null)
  #Cvalue STRAPP
  strapp.all.Cvalue.lambda<-traitDependentBAMM(Div_edata, traits = all.Cvalue_data, 1000, rate = 'speciation', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.Cvalue.mu<-traitDependentBAMM(Div_edata, traits = all.Cvalue_data, 1000, rate = 'extinction', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.Cvalue.div<-traitDependentBAMM(Div_edata, traits = all.Cvalue_data, 1000, rate = 'net diversification', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.Cvalue.lambda.diff<-abs(strapp.all.Cvalue.lambda$obs.corr)-abs(strapp.all.Cvalue.lambda$null)
  strapp.all.Cvalue.mu.diff<-abs(strapp.all.Cvalue.mu$obs.corr)-abs(strapp.all.Cvalue.mu$null)
  strapp.all.Cvalue.div.diff<-abs(strapp.all.Cvalue.div$obs.corr)-abs(strapp.all.Cvalue.div$null)
  #LifeCycle.Strict STRAPP
  strapp.all.LifeCycle.strict.lambda<-traitDependentBAMM(Div_edata, traits = all.LifeCycle_data, 1000, rate = 'speciation', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.LifeCycle.strict.mu<-traitDependentBAMM(Div_edata, traits = all.LifeCycle_data, 1000, rate = 'extinction', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.LifeCycle.strict.div<-traitDependentBAMM(Div_edata, traits = all.LifeCycle_data, 1000, rate = 'net diversification', return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  strapp.all.LifeCycle.strict.lambda.diff<-abs(strapp.all.LifeCycle.strict.lambda$obs.corr)-abs(strapp.all.LifeCycle.strict.lambda$null)
  strapp.all.LifeCycle.strict.mu.diff<-abs(strapp.all.LifeCycle.strict.mu$obs.corr)-abs(strapp.all.LifeCycle.strict.mu$null)
  strapp.all.LifeCycle.strict.div.diff<-abs(strapp.all.LifeCycle.strict.div$obs.corr)-abs(strapp.all.LifeCycle.strict.div$null)
  #table for seed results
  seed.results<-data.frame(cbind(c(strapp.all.seed.lambda$estimate,strapp.all.seed.mu$estimate,strapp.all.seed.div$estimate),c(strapp.all.seed.lambda$p.value,strapp.all.seed.mu$p.value,strapp.all.seed.div$p.value)))
  row.names(seed.results)<-c('lambda','mu','div')
  colnames(seed.results)<-c('estimate','p-value')
  #table for Cvalue results
  Cvalue.results<-data.frame(cbind(c(strapp.all.Cvalue.lambda$estimate,strapp.all.Cvalue.mu$estimate,strapp.all.Cvalue.div$estimate),c(strapp.all.Cvalue.lambda$p.value,strapp.all.Cvalue.mu$p.value,strapp.all.Cvalue.div$p.value)))
  row.names(Cvalue.results)<-c('lambda','mu','div')
  colnames(Cvalue.results)<-c('estimate','p-value')
  #table for LifeCycle results
  LifeCycle.results<-data.frame(cbind(c(strapp.all.LifeCycle.strict.lambda$estimate,strapp.all.LifeCycle.strict.mu$estimate,strapp.all.LifeCycle.strict.div$estimate),c(strapp.all.LifeCycle.strict.lambda$p.value,strapp.all.LifeCycle.strict.mu$p.value,strapp.all.LifeCycle.strict.div$p.value)))
  row.names(LifeCycle.results)<-c('lambda','mu','div')
  colnames(LifeCycle.results)<-c('estimate','p-value')
  seedcvaluelifecycle.df<-cbind(seed.results,Cvalue.results,LifeCycle.results)
  colnames(seedcvaluelifecycle.df)<-c('seedmass.estimate','seedmass.pvalue','Cvalue.estimate','Cvalue.pvalue','LifeCycle.estimate','LifeCycle.pvalue')
  seedcvaluelifecycle.df$parameter<-row.names(seedcvaluelifecycle.df)
  row.names(seedcvaluelifecycle.df)<-NULL
  seedcvaluelifecycle.df<-seedcvaluelifecycle.df[c(ncol(seedcvaluelifecycle.df),1:(ncol(seedcvaluelifecycle.df)-1))]
  write.table(seedcvaluelifecycle.df,'./output/tables/seedCvalueLifeCycle_STRAPP.txt',quote=F,row.names=F,sep='\t')
  #write plot for lambda
  pdf(file='./output/plots/lambda_3trait_strapp_density.pdf',paper='a4')
  plot(density(strapp.all.seed.lambda.diff),ylim=c(0,10),xlim=c(-0.4,0.3),main=paste('SeedMass,C.value,Life Cycle (n=',length(all.Cvalue_data),')',sep=''),ylab='',xlab='absolute.diff in rho(obs-null)',col='blue')
  lines(density(strapp.all.Cvalue.lambda.diff),col='red')
  lines(density(strapp.all.LifeCycle.strict.lambda.diff),col='dark green')
  abline(v=0,lty=2)
  abline(v=mean(strapp.all.seed.lambda.diff),col='blue',lty=2)
  abline(v=mean(strapp.all.Cvalue.lambda.diff),col='red',lty=2)
  abline(v=mean(strapp.all.LifeCycle.strict.lambda.diff),col='dark green',lty=2)
  legend('topright',lty=1,col=c('red','blue','dark green'),c('C.value','SeedMass','Life Cycle'),cex=.5,bty='n')
  dev.off()
  #write plot for mu
  pdf(file='./output/plots/mu_3trait_strapp_density.pdf',paper='a4')
  plot(density(strapp.all.seed.mu.diff),ylim=c(0,10),xlim=c(-0.4,0.3),main=paste('SeedMass,C.value,Life Cycle (n=',length(all.Cvalue_data),')',sep=''),ylab='',xlab='absolute.diff in rho(obs-null)',col='blue')
  lines(density(strapp.all.Cvalue.mu.diff),col='red')
  lines(density(strapp.all.LifeCycle.strict.mu.diff),col='dark green')
  abline(v=0,lty=2)
  abline(v=mean(strapp.all.seed.mu.diff),col='blue',lty=2)
  abline(v=mean(strapp.all.Cvalue.mu.diff),col='red',lty=2)
  abline(v=mean(strapp.all.LifeCycle.strict.mu.diff),col='dark green',lty=2)
  legend('topright',lty=1,col=c('red','blue','dark green'),c('C.value','SeedMass','Life Cycle'),cex=.5,bty='n')
  dev.off()
  #write plot for div
  pdf(file='./output/plots/div_3trait_strapp_density.pdf',paper='a4')
  plot(density(strapp.all.seed.div.diff),ylim=c(0,10),xlim=c(-0.4,0.3),main=paste('SeedMass,C.value,Life Cycle (n=',length(all.Cvalue_data),')',sep=''),ylab='',xlab='absolute.diff in rho(obs-null)',col='blue')
  lines(density(strapp.all.Cvalue.div.diff),col='red')
  lines(density(strapp.all.LifeCycle.strict.div.diff),col='dark green')
  abline(v=0,lty=2)
  abline(v=mean(strapp.all.seed.div.diff),col='blue',lty=2)
  abline(v=mean(strapp.all.Cvalue.div.diff),col='red',lty=2)
  abline(v=mean(strapp.all.LifeCycle.strict.div.diff),col='dark green',lty=2)
  legend('topright',lty=1,col=c('red','blue','dark green'),c('C.value','SeedMass','Life Cycle'),cex=.5,bty='n')
  dev.off()
}
  