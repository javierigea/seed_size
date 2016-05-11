library(BAMMtools)
library(coda)
library(ape)
library(mvtnorm)
library(lattice)

####function to generate BAMM rates tables
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


#############function to generate trees with branches scaled by BAMM parameter (lambda, mu, r)
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

#####function to generate spearman correlations in cohorts across stapp samples (takes some time), outputs to a table
#traitdata: 'BAMM_trait_rates_seed.txt' (for seed or seedrate) or BAMMPhyndr_GenusSeed_cov_trait_data.txt (for seedcov)
#name='seed', 'seedrate' or 'seedcov'
#parameter='speciation','extinction','net diversification'

generate_correlations_strapp<-function(treefile,eventfile,traitdata,burnin,name,parameter){
  #load tree
  tree<-read.tree(treefile)
  #get event data
  print('Analysing eventfile')
  Div_edata <- getEventData(tree, eventfile, burnin = burnin, nsamples = 5000,type='diversification')
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
    trait_data<-log10(trait.input$cov)
    names(trait_data)<-trait.input$tip.label
  }
  print('Running strapp tests')
  strapp<-traitDependentBAMM(Div_edata, traits = trait_data, 1000, rate = parameter, return.full = TRUE, method = 'spearman', logrates = FALSE, two.tailed = TRUE) 
  #build an empty matrix to store correlations from posterior
  trait_data<-data.frame(names(trait_data),unname(trait_data))
  colnames(trait_data)<-c('V1','V2')
  spearman.regressions.df<-matrix(nrow=length(strapp$gen),ncol=2)
  #for all strapp permutations
  for (i in 1:length(strapp$gen)){
    cat(i, '\n')	 
    #run through posterior
    gen<-strapp$gen[i]
    #get parameters from eventobject
    strapp.eventdata <- cbind(Div_edata$edge, Div_edata$eventVectors[[gen]])
    strapp.eventdata <- strapp.eventdata[strapp.eventdata[,2] <= length(tree$tip.label), ]
    strapp.eventdata <- data.frame(Div_edata$tip.label, strapp.eventdata[,3], Div_edata$tipLambda[[gen]],Div_edata$tipMu[[gen]],(Div_edata$tipLambda[[gen]]-Div_edata$tipMu[[gen]]))
    colnames(strapp.eventdata)<-c('tip.label','event','tip.lambda','tip.mu','tip.div')
    #add trait data to parameters
    strapp.eventdata<-merge(strapp.eventdata,trait_data,by.x='tip.label','V1')
    if(name=='seedcov'){
      strapp.eventdata<-na.omit(strapp.eventdata)  
    }
    #to plot spearman correlations
    if(name=='seedrate'){
      strapp.eventdata$V2<-log10(strapp.eventdata$V2)
    }
    #this vector determines the number of correlations to test (we need to get the value where the correlation equals 0)
    #see http://stats.stackexchange.com/questions/64938/if-linear-regression-is-related-to-pearsons-correlation-are-there-any-regressi/110112#110112
    x<-seq(-10,10,by=0.01)
    beta<-vector('numeric')
    if(parameter=='speciation'){
      for (a in 1:length(x)){
        beta[a]<-cor(log10(strapp.eventdata$tip.lambda)-(x[a]*strapp.eventdata$V2),strapp.eventdata$V2,method='s')  
      }
    }
    if(parameter=='extinction'){
      for (a in 1:length(x)){
        beta[a]<-cor(log10(strapp.eventdata$tip.mu)-(x[a]*strapp.eventdata$V2),strapp.eventdata$V2,method='spearman')  
      }
    }
    if(parameter=='net diversification'){
      #this removes any negative values from diversification (they're only 0.3% of the rows approx)
      strapp.eventdata<-strapp.eventdata[strapp.eventdata$tip.div>0,]
      for (a in 1:length(x)){
        beta[a]<-cor(log10(strapp.eventdata$tip.div)-(x[a]*strapp.eventdata$V2),strapp.eventdata$V2,method='spearman')  
      }
    }
    
    df<-data.frame(x,beta)
    #get the value of beta that's closest to 0
    if(name=='seedcov'){
      df<-na.omit(df)
      #get the value of beta that's closest to 0
      value<-which.min(abs(df$beta - 0))
    }
    else {
      value<-which.min(abs(beta - 0))  
    }
    slope<-df[df$beta==df$beta[value],]$x
    if (parameter=='speciation'){
      model<-lm(log10(strapp.eventdata$tip.lambda)-(slope*strapp.eventdata$V2)~0)  
    }
    if (parameter=='extinction'){
      model<-lm(log10(strapp.eventdata$tip.mu)-(slope*strapp.eventdata$V2)~0)  
    }
    if (parameter=='net diversification'){
      model<-lm(log10(strapp.eventdata$tip.div)-(slope*strapp.eventdata$V2)~0)  
    }
    intercept<-median(model$residuals)
    spearman.regressions.df[i,1]<-intercept
    spearman.regressions.df[i,2]<-slope 
  }    
  spearman.regressions.df<-data.frame(spearman.regressions.df)
  if (parameter=='speciation'){
    output.parameter<-'lambda'
  }
  if (parameter=='extinction'){
    output.parameter<-'mu'
  }
  if (parameter=='net diversification'){
    output.parameter<-'div'
  }
  correlationfile<-paste('./output/tables/strapp_all_',output.parameter,'_',name,'.txt',sep='')
  write.table(spearman.regressions.df,correlationfile,row.names=F,sep='\t',quote=F)
###pdffile<-paste('./output/plots/strapp_all_',output.parameter,'_',name,'.pdf',sep='')
###pdf(pdffile,paper='a4')
###if(parameter=='speciation'){
###  plot(strapp.eventdata$V2,log10(strapp.eventdata$tip.lambda),type='n',xlab=name,ylab=output.parameter)
###  col=rgb(214,142,148,maxColorValue = 255)
###}
###if(parameter=='extinction'){
###  plot(strapp.eventdata$V2,log10(strapp.eventdata$tip.mu),type='n',xlab=name,ylab=output.parameter)
###  col=rgb(150,192,231,maxColorValue = 255)
###}
###if(parameter=='net diversification'){
###  plot(strapp.eventdata$V2,log10(strapp.eventdata$tip.div),type='n',xlab=name,ylab=output.parameter)
###  col=rgb(216,198,222,maxColorValue = 255)
###}
###apply(spearman.regressions.df,1,function(x) abline(x[1],x[2],col=col))
###median.slope<-which.min(abs(median(spearman.regressions.df$X2)-spearman.regressions.df$X2))
###abline(spearman.regressions.df$X1[median.slope],spearman.regressions.df$X2[median.slope])  
  ###dev.off()
}
