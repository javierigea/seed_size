library(BAMMtools)
library(coda)
library(ape)
library(mvtnorm)
library(colorRamps)
library(Hmisc)
library(plotrix)
library(scales)
library(gplots)
library(caper)


#############plot phylogeny + "heatmap" variables
#treefile:genus level tree
#tablefile:BAMMdivtrait_rates_seed.txt
#name: to attach to file
plot_phylogeny_heatmap<-function(treefile,tablefile,name){
  tree<-read.tree(treefile)
  hc <- as.hclust(tree) 
  dend <- as.dendrogram(hc)
  #tree$tip.label<-as.character(tree$tip.label)
  table<-read.table(tablefile,sep='\t',header=T)
  table<-table[match(tree$tip.label,table$names),]
  #check the order is correct
  identical(as.character(table$names),as.character(tree$tip.label))
  #plot the tree+bar
  #pdf(pdffile,paper='a4')
  #par(fig=c(0,0.3,0,0.9))
  #plot(tree,show.tip.label = F,edge.width=.2)
  palette <- colorRampPalette(c("blue", "white", "red"))(n = 300)
  matrix<-data.frame(table$log.Seed.Weight,log10(table$mean.beta.rates),log10(table$mean.lambda.rates),log10(table$mean.mu.rates),log10(table$mean.div.rates))
  matrix<-data.matrix(matrix)
  colnames(matrix)<-c('Log10.Seed.Weight','mean.beta.rates','mean.lambda.rates','mean.mu.rates','mean.div.rates')
  row.names(matrix)<-table$names
  #par(fig=c(0.515,1,0,0.9))
  pdffile<-paste('./output/plots/',name,'.pdf',sep='')
  pdf(pdffile,paper='a4')
  heatmap.2(matrix,Rowv=dend,Colv = NA,col = palette,dendrogram='row', scale="column",labRow=FALSE,cexCol=0.5,labCol=c('seed mass','seed mass evolution','lambda','mu','r'),key=TRUE,trace='none',lwid=c(1,0.5),srtCol=45,keysize = 1,denscol = NA,key.title=NA,key.xlab='Z-score',cex=.7,key.ylab=NA)
  dev.off()
  
}

######plot number of rate shifts in BAMM
#mcmcout<-mcmc file from BAMM
#burnin
#name: to attach to output file
plot_number_shifts_BAMM<-function(mcmcout,burnin,name){
  mcmcout.div <- read.csv(mcmcout, header=T)
  plot(mcmcout.div$logLik ~ mcmcout.div$generation)
  #remove initial 10% (burnin)
  burnstart.div <- floor(0.2 * nrow(mcmcout.div))
  postburn.div <- mcmcout.div[burnstart.div:nrow(mcmcout.div), ]
  #compute posterior probabilities of models
  post_probs.div <- table(postburn.div$N_shifts) / nrow(postburn.div)
  names(post_probs.div)
  #plot the posterior distribution of Nrateshifts
  df<-as.data.frame(post_probs.div)
  df$Var1<-as.numeric(as.character(df$Var1))
  pdf(paste('./output/plots/rateshifts_',name,'.pdf',sep=''),paper='a4')
  plot(df$Var1,df$Freq,xlim=c(0,220),yaxs='i',xaxs='i',ylim=c(0,0.05),type='l',ylab='Posterior probability',xlab='number of rate shifts',las=1, main=name)
  dev.off()
}


###############plot density of strapp correlations
##strapp:strapp object
##name:name of parameter (lambda,mu,div)
plot_strapp_density<-function(strapp,name){
  pvalue<-strapp$p.value
  estimate<-strapp$estimate
  strapp.diff<-abs(strapp$obs.corr)-abs(strapp$null)
  title<-paste('rho = ',estimate,'(',pvalue,')',sep='')
  pdffile<-paste('./output/plots/',name,'_strapp_density.pdf',sep='')
  pdf(pdffile,paper='a4') 
  d<-density(strapp.diff)
  plot(d,xlim=c(-0.2,0.6),xaxs='i',yaxs='i',main=title,yaxt='n',ylab='',bty='n',xlab='rho[observed-null]')
  abline(v=0,lty=2)
  dev.off()
}

##############plot strapp correlations across posterior from file
##correlationfile: strapp_all_*.txt file
##parametername: 'lambda','mu' or 'div'
##name: 'seed' or 'seedrate' or 'seedcov'
#tablefile: BAMMdivtrait_rates_seed.txt (for seed and seedrate) or BAMMPhyndr_GenusSeed_cov_trait_data.txt (for cov), just to get limits of plot
plot_strapp_samples_from_file<-function(correlationfile,parameter,name,tablefile){
  spearman.regressions.df<-read.table(correlationfile,header=T,sep='\t')
  table<-read.table(tablefile,header=T,sep='\t')
  pdffile<-paste('./output/plots/strapp_all_main_',parameter,'_',name,'.pdf',sep='')
  pdf(pdffile,paper='a4')
  if (name=='seedrate'){
    min.y<-min(c(table$mean.lambda.rates,table$mean.mu.rates))
    max.y<-max(c(table$mean.lambda.rates,table$mean.mu.rates))
    min.x<-min(c(table$mean.beta.rates))
    max.x<-max(c(table$mean.beta.rates))
    if(parameter=='lambda'){
      plot(c(log10(min.x),log10(max.x)),c(log10(min.y),log10(max.y)),xlab='rate of seed mass change',ylab='speciation rate (lineages myr-1)',xaxt='n',yaxt='n',type='n')
      ticks.x<-seq(-3, 0, by = 1)
      axis(1, at = ticks.x, labels=10^(ticks.x), las=1)
      ticks.y<-seq(-3, 1, by = 1)
      axis(2, at = ticks.y, labels=10^(ticks.y), las=1)
      col.lines<-rgb(214,142,148,maxColorValue = 255)
    }
    if(parameter=='mu'){
      plot(c(log10(min.x),log10(max.x)),c(log10(min.y),log10(max.y)),xlab='rate of seed mass change',ylab='extinction rate (lineages myr-1)',xaxt='n',yaxt='n',type='n')
      ticks.x<-seq(-3, 0, by = 1)
      axis(1, at = ticks.x, labels=10^(ticks.x), las=1)
      ticks.y<-seq(-3, 1, by = 1)
      axis(2, at = ticks.y, labels=10^(ticks.y), las=1)
      col.lines<-rgb(150,192,231,maxColorValue = 255)
    }
    if(parameter=='div'){
      plot(c(log10(min.x),log10(max.x)),c(log10(min.y),log10(max.y)),xlab='rate of seed mass change',ylab='diversification rate (lineages myr-1)',xaxt='n',yaxt='n',type='n')
      ticks.x<-seq(-3, 0, by = 1)
      axis(1, at = ticks.x, labels=10^(ticks.x), las=1)
      ticks.y<-seq(-3, 1, by = 1)
      axis(2, at = ticks.y, labels=10^(ticks.y), las=1)
      col.lines<-rgb(216,198,222,maxColorValue = 255)
    }
  }
  if (name=='seed'){
    min.y<-min(c(table$mean.lambda.rates,table$mean.mu.rates))
    max.y<-max(c(table$mean.lambda.rates,table$mean.mu.rates))
    min.x<-min(c(table$log.Seed.Weight))
    max.x<-max(c(table$log.Seed.Weight))
    if(parameter=='lambda'){
      plot(c(min.x,max.x),c(log10(min.y),log10(max.y)),xlab='seed mass (g)',ylab='speciation rate (lineages myr-1)',xaxt='n',yaxt='n',type='n')
      ticks.x<-seq(-6, 4, by = 2)
      axis(1, at = ticks.x, labels=10^(ticks.x), las=1)
      ticks.y<-seq(-3, 1, by = 1)
      axis(2, at = ticks.y, labels=10^(ticks.y), las=1)
      col.lines<-rgb(214,142,148,maxColorValue = 255)
    }
    if(parameter=='mu'){
      plot(c(min.x,max.x),c(log10(min.y),log10(max.y)),xlab='seed mass (g)',ylab='extinction rate (lineages myr-1)',xaxt='n',yaxt='n',type='n')
      ticks.x<-seq(-6, 4, by = 2)
      axis(1, at = ticks.x, labels=10^(ticks.x), las=1)
      ticks.y<-seq(-3, 1, by = 1)
      axis(2, at = ticks.y, labels=10^(ticks.y), las=1)      
      col.lines<-rgb(150,192,231,maxColorValue = 255)
    }
    if(parameter=='div'){
      plot(c(min.x,max.x),c(log10(min.y),log10(max.y)),xlab='seed mass (g)',ylab='diversification rate (lineages myr-1)',xaxt='n',yaxt='n',type='n')
      ticks.x<-seq(-6, 4, by = 2)
      axis(1, at = ticks.x, labels=10^(ticks.x), las=1)
      ticks.y<-seq(-3, 1, by = 1)
      axis(2, at = ticks.y, labels=10^(ticks.y), las=1)      
      col.lines<-rgb(216,198,222,maxColorValue = 255)
    }
  
    
  }
  if (name=='seedcov'){
    min.y<-min(c(table$mean.lambda.rates,table$mean.mu.rates))
    max.y<-max(c(table$mean.lambda.rates,table$mean.mu.rates))
    min.x<-log10(0.1)
    max.x<-log10(1000)
    if(parameter=='lambda'){
      plot(c(min.x,max.x),c(log10(min.y),log10(max.y)),xlab='seed mass (g)',ylab='speciation rate (lineages myr-1)',xaxt='n',yaxt='n',type='n')
      ticks.x<-seq(-1, 3, by = 1)
      axis(1, at = ticks.x, labels=10^(ticks.x), las=1)
      ticks.y<-seq(-3, 1, by = 1)
      axis(2, at = ticks.y, labels=10^(ticks.y), las=1)
      col.lines<-rgb(214,142,148,maxColorValue = 255)
    }
    if(parameter=='mu'){
      plot(c(min.x,max.x),c(log10(min.y),log10(max.y)),xlab='seed mass (g)',ylab='extinction rate (lineages myr-1)',xaxt='n',yaxt='n',type='n')
      ticks.x<-seq(-1, 3, by = 1)
      axis(1, at = ticks.x, labels=10^(ticks.x), las=1)
      ticks.y<-seq(-3, 1, by = 1)
      axis(2, at = ticks.y, labels=10^(ticks.y), las=1)      
      col.lines<-rgb(150,192,231,maxColorValue = 255)
    }
    if(parameter=='div'){
      plot(c(min.x,max.x),c(log10(min.y),log10(max.y)),xlab='seed mass (g)',ylab='diversification rate (lineages myr-1)',xaxt='n',yaxt='n',type='n')
      ticks.x<-seq(-1, 3, by = 1)
      axis(1, at = ticks.x, labels=10^(ticks.x), las=1)
      ticks.y<-seq(-3, 1, by = 1)
      axis(2, at = ticks.y, labels=10^(ticks.y), las=1)      
      col.lines<-rgb(216,198,222,maxColorValue = 255)
    }
    
    
  }  
  
  
  apply(spearman.regressions.df,1,function(x) abline(x[1],x[2],col=col.lines))
  median.slope<-which.min(abs(median(spearman.regressions.df$X2)-spearman.regressions.df$X2))
  abline(spearman.regressions.df$X1[median.slope],spearman.regressions.df$X2[median.slope])  
  dev.off()
}
###plot genus mean seed mass vs coefficient of variation
#traitdatacov:'BAMMPhyndr_GenusSeed_cov_trait_data.txt'
plot_mean_cov_seed_mass<-function(traitdatacov){
  trait_input.cov<-read.table(traitdatacov,sep='\t',header=T)
  trait_input.cov<-na.omit(trait_input.cov)
  trait_input.cov$mean<-log10(trait_input.cov$mean)
  pdf('./output/plots/seed_mean_vs_cov.pdf',paper='a4') 
  plot(trait_input.cov$mean,log10(trait_input.cov$cov),las=1,xlab='seed mass(g.)',ylab='seed mass coefficient of variation',xaxt='n',yaxt='n')
  ticks.x<-seq(-6, 4, by = 2)
  axis(1, at = ticks.x, labels=10^(ticks.x), las=1)
  ticks.y<-seq(-1, 3, by = 1)
  axis(2, at = ticks.y, labels=10^(ticks.y), las=1)      
  dev.off()
}


####plot family level tree with family level incomplete sampling bar graph
#treefile: genus level tree
#tablefile: Phyndr_BAMM_Genus_SamplingFractions.txt (table with Genus\tFamily\tSamplingFraction)
plot_famtree_sampling<-function(treefile,tablefile){
  pdffile<-('./output/plots/familylevel_sampling.pdf')
  tree<-read.tree(treefile)
  tree$tip.label<-as.character(tree$tip.label)
  table<-read.table(tablefile,sep='\t',header=T)
  table<-table[match(tree$tip.label,table$Genus),]
  #check the order is correct
  identical(as.character(table$Genus),as.character(tree$tip.label))
  #get one tip per family
  tree$tip.label<-as.character(paste(table$family,table$Genus,sep='_'))
  family<-unique(sapply(strsplit(as.character(tree$tip.label),'_'),function(x) x[1]))
  family2<-sapply(family, function(x) x<- paste(x,'_',sep=''))
  ii<-sapply(family,function(x,y) grep(x,y,fixed=TRUE)[1],y=tree$tip.label)
  famtree<-drop.tip(tree,setdiff(tree$tip.label,tree$tip.label[ii]))
  #remove genus name from tips
  famtree$tip.label<-sapply(strsplit(famtree$tip.label,"_"),function(x) x[1])
  samplingdf<-data.frame(table$family,table$Sampling_Fraction)
  samplingdf<-unique(samplingdf)
  sampling<-samplingdf$table.Sampling_Fraction
  names(sampling)<-samplingdf$table.family
  #plot the tree+bar
  pdf(pdffile,paper='a4')
  par(fig=c(0,0.7,0,0.9))
  plot(famtree,show.tip.label = F,edge.width=.2)
  axisPhylo()
  mtext("(a)", adj = 0.05, padj = -3, cex = 0.9)
  mtext("Millions of Years", side=1, adj = 0.5, padj = 4)
  par(fig=c(0.515,1,0,0.9), new=TRUE)
  barplot(sampling[famtree$tip.label],horiz=TRUE,width=1,space=0,ylim=c(1,length(famtree$tip.label))-0.5,names="",col = 'red', axes = FALSE, xlab = 'Percentage of \ngenera sampled')
  axis(2, at=c(-10,310), labels = c('',''), hadj = -0.5 , las = 2, tick = TRUE, tcl=0)
  axis(4, at=c(-10,310), labels = c('',''), hadj = 1, las = 2, tick =TRUE, tcl=0)
  abline(v=0.5, lty = 2, col='gray80')
  mtext("0%", adj = -0.02, padj = -0.1, cex = 1)
  mtext("50%", adj = 0.53, padj = -0.1, cex = 1, col='gray55')
  mtext("100%", adj = 1.05, padj = -0.1, cex = 1)
  mtext("(b)", adj = 0.05, padj = -3, cex = 0.9)
  dev.off()
}
#####plot type I error results (from 'typeIerror_sims_*.txt')
#simfile: typeIerror_sims_*.txt
plot_typeIerror<-function(simfile){
  sims<-read.table(simfile)
  error<-sum(sims$V1<=0.05)/length(sims$V1)
  h<-hist(sims$V1)
  h$density<-h$counts/sum(h$counts)
  pdf(file='./output/plots/typeIerror_strapp.pdf',paper='a4')
  plot(h,freq=FALSE,ylim=c(0,1),xaxs='i',yaxs='i',xlab='STRAPP test p-value',ylab='Frequency',las=1,main=paste('p-values of STRAPP with neutral traits (Type I error=',error,')',sep = ''))
  dev.off()
}

####plot rate-scaled tree with branches coloured by rates of trait evolution
#parameter tree: lambda, mu or div  tree (the tree file that will determine branch lengths)
#colourtee: rate of evolution tree (the tree file that will determine branch colours)
#name: name to append to pdf
plot_colour_ratetree<-function(parametertree,colourtree,name){
  #read speciation rate and trait rate trees (WITH BRANCHES TRANSFORMED FOR RATES)
  parametertree<-read.tree(parametertree)
  colourtree<-read.tree(colourtree)
  #get a vector with the percentiles
  quantiles.trait<-quantile(colourtree$edge.length, c(0,.25,.50,.75,1))
  #create an empty vector to store the colours of the branches
  branchcolours <- vector(mode="character", length=length(colourtree$edge.length))
  #creat a colour palette
  colours<-c("#009933","#ccff99","#ff9933","#ff0000")
  #loop throught the trait rate branches and assign colours to each bin 
  for (i in 1:length(colourtree$edge.length)){
    #from 0 to 0.25
    if (colourtree$edge.length[i] <= quantiles.trait[2]){
      branchcolours[i]<-colours[1]
    }
    #from 0.25 to 0.5
    if ((colourtree$edge.length[i] > quantiles.trait[2]) &&  (colourtree$edge.length[i] <= quantiles.trait[3])){
      branchcolours[i]<-colours[2]
    }
    #from 0.5 to 0.75
    if ((colourtree$edge.length[i] > quantiles.trait[3]) &&  (colourtree$edge.length[i] < quantiles.trait[4])){
      branchcolours[i]<-colours[3]
    } 
    #from 0.75 to 1
    if (colourtree$edge.length[i] >= quantiles.trait[4]) {
      branchcolours[i]<-colours[4]
    } 
  }
  #plot the tree with the colours
  pdffile<-paste('./output/plots/',name,'_trait_paintedtree.pdf',sep='')
  pdf(pdffile,paper='a4')
  plot(parametertree,edge.color=branchcolours,show.tip.label = F,edge.width = .3)
  add.scale.bar(1,1,length = 0.5,cex=.5)
  legend(title= 'Relative rates of seed mass evolution', x='right',inset = 0.065, legend = c('First Quartile \n<0.024','\nSecond Quartile \n(0.24-0.034)', '\nThird Quartile \n(0.034-0.095)', '\nFourth Quartile \n(>0.095)'),bty='n',fill = c(colours[1],colours[2], colours[3],colours[4]),cex=.5)
  dev.off()
  #########to get clade names for labelling
  ######table<-read.table('./output/tables/table_Angiosperm_clades.txt',sep='\t')
  ###Rosids<-table[table$V1=='Rosids',]$genus
  ###Rosids<-subset(table,table$V1=='Rosids')
  ###Rosids<-Rosids$genus
  ###Asterids<-table[table$V1=='Asterids',]$genus
  ###Asterids<-subset(table,table$V1=='Asterids')
  ###Asterids<-Asterids$genus
  ###Monocots<-table[table$V1=='Monocots',]$genus
  ###Monocots<-subset(table,table$V1=='Monocots')
  ###Monocots<-Monocots$genus
  ###Magnoliids<-table[table$V1=='Magnoliids',]$genus
  ###Magnoliids<-subset(table,table$V1=='Magnoliids')
  ###Magnoliids<-Magnoliids$genus
  ###nodesRosids<-which(parametertree$tip.label %in% Rosids)
  ###nodesAsterids<-which(parametertree$tip.label %in% Asterids)
  ###nodesMonocots<-which(parametertree$tip.label %in% Monocots)
  ###nodesMagnoliids<-which(parametertree$tip.label %in% Magnoliids)
  ###plot(parametertree,edge.color=branchcolours,show.tip.label = F,edge.width = .3)
  ###nodelabels(tip=nodesRosids,cex=.3)
  ###plot(parametertree,edge.color=branchcolours,show.tip.label = F,edge.width = .3)
  ###tiplabels(tip=nodesAsterids,cex=.3)
  ###plot(parametertree,edge.color=branchcolours,show.tip.label = F,edge.width = .3)
  ###tiplabels(tip=nodesMonocots,cex=.3)
  ###plot(parametertree,edge.color=branchcolours,show.tip.label = F,edge.width = .3)
  ###tiplabels(tip=nodesMagnoliids,cex=.3)
}