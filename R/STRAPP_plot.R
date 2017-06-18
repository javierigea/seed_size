library(BAMMtools)
library(coda)
library(ape)
library(mvtnorm)
library(lattice)

#runs STRAPP and outputs density and tables
#traitdata: 'BAMM_trait_rates_seed.txt'
#name:'seed' or 'seedrate'
STRAPP_plot<-function(treefile,eventfile,traitdata,burnin,name){
  #load tree
  tree<-read.tree(treefile)
  #get event data
  print('Analysing eventfile')
  if (class(eventfile)!='bammdata'){
    Div_edata <- getEventData(tree, eventfile, burnin = burnin, nsamples = 1000,type='diversification')
  }else if (class(eventfile)=='bammdata'){
    Div_edata<-eventfile
  }
    trait.input<-read.table(traitdata,header=T,sep='\t')
  if(name=='seed'){
    #read trait input
    trait_data<-trait.input$log.Seed.Weight
    names(trait_data)<-trait.input$names
  }
  if(name=='seedrate'){
    #read trait input
    trait.input<-read.table(traitdata,header=T,sep='\t')
    trait_data<-trait.input$mean.beta.rates
    names(trait_data)<-trait.input$names
  }
  if(name=='seedcov'){
    #read trait input
    trait.input<-read.table(traitdata,header=T,sep='\t')
    trait_data<-trait.input$cov
    names(trait_data)<-trait.input$tip.label
  }
  print('Running strapp tests')
  strapp.lambda<-traitDependentBAMM(Div_edata, traits = trait_data, 5000, rate = 'speciation', return.full = TRUE, method = 'spearman', logrates = TRUE, two.tailed = TRUE) 
  strapp.mu<-traitDependentBAMM(Div_edata, traits = trait_data, 5000, rate = 'extinction', return.full = TRUE, method = 'spearman', logrates = TRUE, two.tailed = TRUE) 
  strapp.div<-traitDependentBAMM(Div_edata, traits = trait_data, 5000, rate = 'net diversification', return.full = TRUE, method = 'spearman', logrates = TRUE, two.tailed = TRUE) 
  print('Plotting strapp histograms')
  plot_strapp_density(strapp.lambda, paste('lambda_',name,sep=''))
  plot_strapp_density(strapp.mu,paste('mu_',name,sep=''))
  plot_strapp_density(strapp.div,paste('div_',name,sep=''))
  print('Outputting table to file')
  strapp.results<-data.frame(c('lambda','mu','div'),c(strapp.lambda$estimate,strapp.mu$estimate,strapp.div$estimate),c(strapp.lambda$p.value,strapp.mu$p.value,strapp.div$p.value))
  colnames(strapp.results)<-c('parameter','estimate','pvalue')
  strappfile<-paste('./output/tables/strappresults_',name,'.txt',sep='')
  write.table(strapp.results,strappfile,sep='\t',quote=F,row.names=F)
}

  
