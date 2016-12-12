library(picante)
library(BAMMtools)
library(RPANDA)
library(geiger)
library(taxonlookup)
library(ape)
library(parallel)
library(phytools)
library(PBD)
library(caper)
library(MuMIn)
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#from Rabosky 2016 "Challenges in estimation" in the supplementary info
###   This equation estimates the probability that a sample of k taxa from a clade of n total taxa
###       includes the root node, under a Yule process.
###
###   The equation is taken from:
###
###   Sanderson, M. J. 1996. How many taxa must be sampled to identify the root node of a large clade?
###                    Systematic Biology 45: 168 - 173
crownCaptureProbability <- function(n, k){
  p1 <- 2*(n-k)/((n-1)*(k+1));
  return(1-p1);
}

########################get phylogroups function by Victor Soria-Carrasco
get.phylogroups<-function(tree, minage, maxage, mincladesize, ncores){
  root.age<-max(node.depth.edgelength(tree))
  nodes<-seq(1,tree$Nnode)+length(tree$tip.label)
  nodes.ages<-as.numeric(unlist(mclapply(nodes, 
                                         function(x) root.age-nodeheight(tree,x), mc.cores=ncores)))
  nodes.desc<-as.numeric(unlist(mclapply(nodes, 
                                         function(x) length(getDescendants(tree,x)), mc.cores=ncores)))
  nodes<-data.frame(node=nodes,age=nodes.ages,ndesc=nodes.desc)
  chosen.nodes<-nodes[(nodes$age>=minage & nodes$age<=maxage &
                         nodes$ndesc>mincladesize),]
  
  # sort nodes by clade size and remove the descendendants
  # in order to keep the most inclusive clades only
  chosen.nodes.sort<-chosen.nodes[order(-chosen.nodes[,3]),]
  for(i in 1:nrow(chosen.nodes.sort)){
    if (!is.na(chosen.nodes.sort[i,1])){
      desc<-getDescendants(tree, chosen.nodes.sort[i,1])
      chosen.nodes.sort<-chosen.nodes.sort[!(chosen.nodes.sort[,1] %in% desc),]
    }
  }
  
  # extract clades
  phylogroups<-list()
  for (i in 1:nrow(chosen.nodes.sort)){
    n<-chosen.nodes.sort$node[i]
    phylogroups[[i]]<-extract.clade(tree,n)
  }
  class(phylogroups)<-"multiPhylo"
  
  return(list(nodes=chosen.nodes.sort, phylogroups=phylogroups))
  #return(phylogroups)
}
################################
################################
############################################################################
#runs RPANDA on clade with a set of models, performs model selection and gets the parameters of the best model
#it also runs fitcontinuous for the clade and seed size data
run_RPANDAextended_fitContinuous_phylogroup_select_size<-function(phylogroup,table){
  #cat(length(phylogroup$tip.label),'\n')
  genus.table<-table[match(phylogroup$tip.label,table$tip),]
  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
  genus.table.genera<-unique(genus.table.genera)
  genus.counts<-table(as.character(genus.table$genus))
  genus.counts<-data.frame(names(genus.counts),unname(genus.counts))
  genus.counts<-genus.counts[,c(1,3)]
  colnames(genus.counts)<-c('genus','n.of.tips')
  genus.table.genera<-merge(genus.table.genera,genus.counts)
  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
  genus.table$phylogroup.mean.seed.size<-mean(genus.table$Log10.Seed.Weight)
  genus.table$phylogroup.sd.seed.size<-sd(genus.table$Log10.Seed.Weight)
  genus.table$phylogroup.sampling.fraction<-sampling.fraction
  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
  if(length(unique(genus.table$family))==1){
    #cat(as.character(genus.table$family[1]),'\n')
  }
  tot_time<-max(node.age(phylogroup)$ages)
  f.lamb.cst<-function(t,y){y[1]}
  f.lamb.exp<-function(t,y){y[1]*exp(y[2]*t)}
  f.mu.cst.0<-function(t,y){0}
  f.mu.cst<-function(t,y){y[1]}
  f.mu.exp<-function(t,y){y[1]*exp(y[2]*t);y[1]*exp(y[2]*t)}
  lamb_par_init.exp<-c(0.5,0.01)
  mu_par_init.exp<-c(0.05,-0.01)
  lamb_par_init.cst<-c(0.5)
  mu_par_init.cst<-c(0.05)
  mu_par_init.0<-c()
  res.lambda.cst.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par_init.cst,mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE)
  if(res.lambda.cst.mu.0$lamb_par[1]<0){
    cat('negative','\n')
    res.lambda.cst.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par = c(5),mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE)
  }
  res.lambda.exp.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par_init.exp,mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE)
  if(res.lambda.exp.mu.0$lamb_par[1]<0){
    cat('negative','\n')
    res.lambda.exp.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par=c(5,0.05),mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE)
  }
  res.lambda.cst.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par_init.cst,mu_par_init.cst,f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE)
  if(res.lambda.cst.mu.cst$mu_par[1]<0){
    cat('negative','\n')
    res.lambda.cst.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par=c(5),mu_par=c(0.01),f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE)
  }  
  res.lambda.cst.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par_init.cst,mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE)
  if(res.lambda.cst.mu.exp$mu_par[1]<0){
    cat('negative','\n')
    res.lambda.cst.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par = c(5),mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE)
  }  
  res.lambda.exp.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.exp,lamb_par_init.exp,mu_par_init.exp,f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE)
  if(res.lambda.exp.mu.exp$mu_par[1]<0){
    cat('negative','\n')
    res.lambda.exp.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.exp,lamb_par=c(5,0.05),mu_par=c(0.05,-0.01),f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE)
  }  
  res.lambda.exp.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par_init.exp,mu_par_init.cst,f=sampling.fraction,expo.lamb=TRUE,cst.mu=TRUE)
  if(res.lambda.exp.mu.cst$mu_par[1]<0){
    cat('negative','\n')
    res.lambda.exp.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par=c(5,0.05),mu_par=c(0.01),f=sampling.fraction,expo.lamb=TRUE,cst.mu=TRUE)
  }  
  aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.exp.mu.0$aicc,res.lambda.exp.mu.exp$aicc,res.lambda.exp.mu.cst$aicc,res.lambda.cst.mu.exp$aicc,res.lambda.cst.mu.cst$aicc)
  names(aicc.vector)<-c('lambda.cst.mu0','lambda.exp.mu.0','lambda.exp.mu.exp','lambda.exp.mu.cst','lambda.cst.mu.exp','lambda.cst.mu.cst')
  aicc.weights<-sort(Weights(aicc.vector),decreasing=T)
  if (names(aicc.weights)[1]=='lambda.cst.mu0'){
    genus.table$lambda.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
    genus.table$mu.rpanda1<-0
    genus.table$ndr.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
    genus.table$phylogroup.size<-nrow(genus.table)
  }
  else if (names(aicc.weights)[1]=='lambda.exp.mu.0'){
    genus.table$lambda.rpanda1<-res.lambda.exp.mu.0$lamb_par[1]
    genus.table$mu.rpanda1<-0
    genus.table$ndr.rpanda1<-res.lambda.exp.mu.0$lamb_par[1]
    genus.table$phylogroup.size<-nrow(genus.table)
  }
  else if (names(aicc.weights)[1]=='lambda.exp.mu.exp'){
    genus.table$lambda.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]
    genus.table$mu.rpanda1<-res.lambda.exp.mu.exp$mu_par[1]
    genus.table$ndr.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]-res.lambda.exp.mu.exp$mu_par[1]
    genus.table$phylogroup.size<-nrow(genus.table)
  }
  else if (names(aicc.weights)[1]=='lambda.exp.mu.cst'){
    genus.table$lambda.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]
    genus.table$mu.rpanda1<-res.lambda.exp.mu.cst$mu_par[1]
    genus.table$ndr.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]-res.lambda.exp.mu.cst$mu_par[1]
    genus.table$phylogroup.size<-nrow(genus.table)
  }
  else if (names(aicc.weights)[1]=='lambda.cst.mu.exp'){
    genus.table$lambda.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]
    genus.table$mu.rpanda1<-res.lambda.cst.mu.exp$mu_par[1]
    genus.table$ndr.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]-res.lambda.cst.mu.exp$mu_par[1]
    genus.table$phylogroup.size<-nrow(genus.table)
  }
  else if (names(aicc.weights)[1]=='lambda.cst.mu.cst'){
    genus.table$lambda.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]
    genus.table$mu.rpanda1<-res.lambda.cst.mu.cst$mu_par[1]
    genus.table$ndr.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]-res.lambda.cst.mu.cst$mu_par[1]
    genus.table$phylogroup.size<-nrow(genus.table)
  }
  if (genus.table$mu.rpanda1[1]<0){
    cat('*****NEGATIVE EXTINCTION RATE','\n')
  }
  genus.table$aicc<-names(aicc.weights)[1]
  trait<-genus.table$Log10.Seed.Weight
  names(trait)<-genus.table$tip
  res.trait1<-fitContinuous(phylogroup,trait,model="BM")
  genus.table$sigsq<-res.trait1$opt$sigsq
  genus.table$z0<-res.trait1$opt$z0
  return(genus.table)
}
############################################################################
############################################################################
#calculates lambda with Magallon & Sanderson (with three extinction fractions 0, 0.5, 0.9) for a clade
#(for net diversification remove (1-0.5) and (1-0.9))
#it also runs fitcontinuous for the clade and seed size data
run_MSlambda_fitContinuous_phylogroup_select_size<-function(phylogroup,table){
  #cat(length(phylogroup$tip.label),'\n')
  genus.table<-table[match(phylogroup$tip.label,table$tip),]
  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
  genus.table.genera<-unique(genus.table.genera)
  genus.counts<-table(as.character(genus.table$genus))
  genus.counts<-data.frame(names(genus.counts),unname(genus.counts))
  genus.counts<-genus.counts[,c(1,3)]
  colnames(genus.counts)<-c('genus','n.of.tips')
  genus.table.genera<-merge(genus.table.genera,genus.counts)
  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
  genus.table$phylogroup.mean.seed.size<-mean(genus.table$Log10.Seed.Weight)
  genus.table$phylogroup.sd.seed.size<-sd(genus.table$Log10.Seed.Weight)
  genus.table$phylogroup.sampling.fraction<-sampling.fraction
  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
  if(length(unique(genus.table$family))==1){
    #cat(as.character(genus.table$family[1]),'\n')
  }
  tot_time<-max(node.age(phylogroup)$ages)
  genus.table$lambda.ms0<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0)
  genus.table$lambda.ms05<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.5)/(1-0.5)
  genus.table$lambda.ms09<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.9)/(1-0.9)
  genus.table$phylogroup.size<-nrow(genus.table)
  trait<-genus.table$Log10.Seed.Weight
  names(trait)<-genus.table$tip
  res.trait1<-fitContinuous(phylogroup,trait,model="BM")
  genus.table$sigsq<-res.trait1$opt$sigsq
  genus.table$z0<-res.trait1$opt$z0
  return(genus.table)
}
##########################################################
##########################################################


##########################################################
#####wrapper to run pgls analysis with RPANDA
run_clades_RPANDA_pgls<-function(minage,maxage,mincladesize,ncores,sampling,table){
  colnames(table)[3]<-'tip'
  cat('getting clades','\n')
  phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=mincladesize,ncores=ncores)
  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
  length(unique(phylogroups.species))
  cat('measuring diversification/trait evolution in clades','\n')
  phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDAextended_fitContinuous_phylogroup_select_size(x,table))
  df.phylogroups <- do.call("rbind", phylogroups.RPANDA)
  #dropping tips that are not in the phylogroups
  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroups.species))
  #then keep just one tip per phylogroup (the first)
  phylogroup.clades.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label[1]))
  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,phylogroup.clades.species))
  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
  colnames(names.phylogroups.tree)<-'Genus'
  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Genus',by.y='tip',all.x=TRUE)
  #hist(df.0to5_3$phylogroup.sd.seed.size,xlab='sd.log10.Seed.Weight',main='phylogroup.0to5.sd.seed.size',xaxs='i',breaks=20)
  #adding a step to select clades where seed size doesn't vary too much
  #df.0to5_3<-df.0to5_3[df.0to5_3$phylogroup.sd.seed.size<0.5,]
  #another filter to check that size of clade is correct
  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize,]
  #add another filter to select well sampled clades
  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
  #hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
  #hist(df.phylogroups$phylogroup.mean.seed.size,xlab='mean.log10.Seed.Weight',main='phylogroup.0to5.mean.seed.size',xaxs='i')
  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Genus))
  phylogroups.data<-data.frame(df.phylogroups$Genus,df.phylogroups$phylogroup.mean.seed.size,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$sigsq)
  colnames(phylogroups.data)<-c('Genus','phylogroup.mean.seed.size','lambda.rpanda','mu.rpanda','ndr.rpanda','sigsq')
  phylogroups.tree$node.label<-NULL
  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Genus',vcv=TRUE)
  #for mean family seed weight
  cat('running pgls','\n')
  pgls.seedlambda<-pgls(log10(lambda.rpanda)~phylogroup.mean.seed.size, data = phylogroups.tree.object, lambda='ML')
  summary(pgls.seedlambda)
  #pgls.seedmu<-pgls(log10(mu.rpanda)~phylogroup.mean.seed.size, data = phylogroups.tree.object, lambda='ML')
  #summary(pgls.seedmu)
  pgls.seedndr<-pgls(log10(ndr.rpanda)~phylogroup.mean.seed.size, data = phylogroups.tree.object, lambda='ML')
  summary(pgls.seedndr)
  pgls.lambdasigsq<-pgls(log10(lambda.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML')
  summary(pgls.lambdasigsq)
  #pgls.musigsq<-pgls(log10(mu.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML')
  #summary(pgls.musigsq)
  pgls.ndrsigsq<-pgls(log10(ndr.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML')
  summary(pgls.ndrsigsq)
  pdf(paste('./output/plots/clade_analyses/clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size_seedsize.pdf',sep=''),paper='a4')
  par(mfrow=c(2,2))
  plot(log10(phylogroups.tree.object$data$lambda.rpanda)~phylogroups.tree.object$data$phylogroup.mean.seed.size,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.5,1.1),ylim=c(-3,1))
  abline(pgls.seedlambda)
  #plot(log10(phylogroups.tree.object$data$mu.rpanda)~phylogroups.tree.object$data$phylogroup.mean.seed.size,xlab='log10.Seed.Weight',ylab='log10.mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedmu)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedmu)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
  #abline(pgls.seedmu)
  plot(log10(phylogroups.tree.object$data$ndr.rpanda)~phylogroups.tree.object$data$phylogroup.mean.seed.size,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.5,1.1),ylim=c(-3,1))
  abline(pgls.seedndr)
  dev.off()
  pdf(paste('./output/plots/clade_analyses/clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size_seedsizerate.pdf',sep=''),paper='a4')
  par(mfrow=c(2,2))
  plot(log10(phylogroups.tree.object$data$lambda.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.lambdasigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.lambdasigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3,0.5),ylim=c(-3,1))
  abline(pgls.lambdasigsq)
  #plot(log10(phylogroups.tree.object$data$mu.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.musigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.musigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
  #abline(pgls.musigsq)
  plot(log10(phylogroups.tree.object$data$ndr.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.ndrsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.ndrsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3,0.5),ylim=c(-3,1))
  abline(pgls.ndrsigsq)
  dev.off()
  pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
  #pgls.seedmuresults<-c('pgls.seedmu',round(summary(pgls.seedmu)$coefficients[2,1],8),round(summary(pgls.seedmu)$coefficients[2,4],8))
  pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
  pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.lambdasigsq)$coefficients[2,1],8),round(summary(pgls.lambdasigsq)$coefficients[2,4],8))
  #pgls.seedratemuresults<-c('pgls.seedratemu',round(summary(pgls.musigsq)$coefficients[2,1],8),round(summary(pgls.musigsq)$coefficients[2,4],8))
  pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.ndrsigsq)$coefficients[2,1],8),round(summary(pgls.ndrsigsq)$coefficients[2,4],8))
  results.df.lambda<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
  #results.df.mu<-as.data.frame(rbind(pgls.seedmuresults,pgls.seedratemuresults))
  results.df.ndr<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
  colnames(results.df.lambda)<-c('analysis','slope','pvalue')
  #colnames(results.df.mu)<-c('analysis','slope','pvalue')
  colnames(results.df.ndr)<-c('analysis','slope','pvalue')
  rownames(results.df.lambda)<-NULL
  #rownames(results.df.mu)<-NULL
  rownames(results.df.ndr)<-NULL
  write.table(results.df.lambda,file=paste('./output/plots/clade_analyses/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_lambda_results.txt',sep=''),quote=F,row.names=F,sep='\t')
  write.table(results.df.ndr,file=paste('./output/plots/clade_analyses/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_ndr_results.txt',sep=''),quote=F,row.names=F,sep='\t')
}
##################################################

##################################################
###wrapper to run clade analysis (with Magallon & Sanderson estimator) for congeneric species only
run_clades_MSlambda_congenerics_pgls<-function(minage,maxage,mincladesize,ncores,sampling,table){
  colnames(table)[3]<-'tip'
  cat('getting clades','\n')
  phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=mincladesize,ncores=ncores)
  #check for phylogroups that only contain congeneric species
  phylogroups.congenerics.positions<-unlist(lapply(phylogroups$phylogroups,function(x){genera<-length(unique(sapply(strsplit(as.character(x$tip.label),'_'),function(x) x[1]))==1)}))
  phylogroups.congenerics<-phylogroups$phylogroups[which(phylogroups.congenerics.positions==1)]
  phylogroups.congenerics.species<-unlist(lapply(phylogroups.congenerics,function(x)x$tip.label))
  length(unique(phylogroups.congenerics.species))
  cat('measuring diversification/trait evolution in clades','\n')
  phylogroups.MS<-lapply(phylogroups.congenerics,function(x) run_MS_fitContinuous_phylogroup_select_size(x,table))
  df.phylogroups <- do.call("rbind", phylogroups.MS)
  #dropping tips that are not in the phylogroups
  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroups.congenerics.species))
  #then keep just one tip per phylogroup (the first)
  phylogroup.clades.species<-unlist(lapply(phylogroups.congenerics,function(x)x$tip.label[1]))
  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,phylogroup.clades.species))
  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
  colnames(names.phylogroups.tree)<-'Genus'
  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Genus',by.y='tip',all.x=TRUE)
  pdf(paste('./output/plots/clade_analyses/MS_congenerics_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
  par(mfrow=c(2,2))
  #another filter to check that size of clade is correct
  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize,]
  #add another filter to select well sampled clades
  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
  hist(df.phylogroups$phylogroup.mean.seed.size,xlab='mean.log10.Seed.Weight',main='phylogroup.0to5.mean.seed.size',xaxs='i')
  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Genus))
  phylogroups.data<-data.frame(df.phylogroups$Genus,df.phylogroups$phylogroup.mean.seed.size,df.phylogroups$ms0,df.phylogroups$ms05,df.phylogroups$ms09,df.phylogroups$sigsq)
  colnames(phylogroups.data)<-c('Genus','phylogroup.mean.seed.size','ms0','ms05','ms09','sigsq')
  phylogroups.tree$node.label<-NULL
  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Genus',vcv=TRUE)
  #for mean family seed weight
  cat('running pgls','\n')
  pgls.seedlambda<-pgls(log10(ms09)~phylogroup.mean.seed.size, data = phylogroups.tree.object, lambda='ML')
  summary(pgls.seedlambda)
  pgls.seedsigsq<-pgls(log10(ms09)~log10(sigsq), data = phylogroups.tree.object, lambda='ML')
  summary(pgls.seedsigsq)
  plot(log10(phylogroups.tree.object$data$ms09)~phylogroups.tree.object$data$phylogroup.mean.seed.size,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,0.5))
  abline(pgls.seedlambda)
  plot(log10(phylogroups.tree.object$data$ms09)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,0.5))
  abline(pgls.seedsigsq)
  dev.off()
  pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
  pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
  colnames(results.df)<-c('analysis','slope','pvalue')
  rownames(results.df)<-NULL
  write.table(results.df,file=paste(minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MScongenerics_results.txt',sep=''),quote=F,row.names=F,sep='\t')
}

##################################################

##################################################
###wrapper to run clade analysis (with Magallon & Sanderson estimator)
run_clades_MSlambda_pgls<-function(minage,maxage,mincladesize,ncores,sampling,table){
  colnames(table)[3]<-'tip'
  cat('getting clades','\n')
  phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=mincladesize,ncores=ncores)
  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
  length(unique(phylogroups.species))
  cat('measuring diversification/trait evolution','\n')
  phylogroups.MS<-lapply(phylogroups$phylogroups,function(x) run_MSlambda_fitContinuous_phylogroup_select_size(x,table))
  df.phylogroups <- do.call("rbind", phylogroups.MS)
  #dropping tips that are not in the phylogroups
  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroups.species))
  #then keep just one tip per phylogroup (the first)
  phylogroup.clades.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label[1]))
  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,phylogroup.clades.species))
  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
  colnames(names.phylogroups.tree)<-'Genus'
  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Genus',by.y='tip',all.x=TRUE)
  pdf(paste('./output/plots/clade_analyses/MSlambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
  par(mfrow=c(2,2))
  #another filter to check that size of clade is correct
  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize,]
  #add another filter to select well sampled clades
  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
  hist(df.phylogroups$phylogroup.mean.seed.size,xlab='mean.log10.Seed.Weight',main='phylogroup.0to5.mean.seed.size',xaxs='i')
  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Genus))
  phylogroups.data<-data.frame(df.phylogroups$Genus,df.phylogroups$phylogroup.mean.seed.size,df.phylogroups$lambda.ms0,df.phylogroups$lambda.ms05,df.phylogroups$lambda.ms09,df.phylogroups$sigsq)
  colnames(phylogroups.data)<-c('Genus','phylogroup.mean.seed.size','lambda.ms0','lambda.ms05','lambda.ms09','sigsq')
  phylogroups.tree$node.label<-NULL
  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Genus',vcv=TRUE)
  #for mean family seed weight
  cat('running pgls','\n')
  pgls.seedlambda<-pgls(log10(lambda.ms05)~phylogroup.mean.seed.size, data = phylogroups.tree.object, lambda='ML')
  summary(pgls.seedlambda)
  pgls.seedsigsq<-pgls(log10(lambda.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML')
  summary(pgls.seedsigsq)
  plot(log10(phylogroups.tree.object$data$lambda.ms05)~phylogroups.tree.object$data$phylogroup.mean.seed.size,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
  abline(pgls.seedlambda)
  plot(log10(phylogroups.tree.object$data$lambda.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
  abline(pgls.seedsigsq)
  dev.off()
  pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
  pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
  colnames(results.df)<-c('analysis','slope','pvalue')
  rownames(results.df)<-NULL
  write.table(results.df,file=paste('./output/plots/clade_analyses/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSlambda_results.txt',sep=''),quote=F,row.names=F,sep='\t')
}

######################################################
cat('preparing dataset to run clade based analyses','\n')
#Qian and Jin created an updated version of Zanne tree - 'Appendix S3-v2' or 'QianTree.txt'
#reference doi: 10.1093/jpe/rtv047
Qian<-read.tree('./raw_data/QianTree.txt')
Qian_TPL<-read.table('./output/tables/Qian_TPL.txt',header=T,sep='\t',stringsAsFactors = F)
##fix a mistake with Ranunculus
Qian_TPL[Qian_TPL$Merged=='Ranunculus_maclovianus',]$Merged<-'Carex_fuscula'
#run through taxonlookup
Qian_TPLunique<-Qian_TPL[!duplicated(Qian_TPL$Merged),]
Qian_Lookup<-lookup_table(as.character(Qian_TPLunique$Merged),by_species = TRUE)
#get table with angiosperms only
Qian_Lookup_Angios<-subset(Qian_Lookup,Qian_Lookup$group=='Angiosperms')
Qian_Lookup_Angios$Fullspecies<-row.names(Qian_Lookup_Angios)
row.names(Qian_Lookup_Angios)<-NULL
Qian_Lookup_Angios_TPL<-merge(Qian_Lookup_Angios,Qian_TPLunique,by.x='Fullspecies',by.y='Merged')
colnames(Qian_Lookup_Angios_TPL)[1]<-'New_Species'
Qian_Lookup_Angios_TPL$Old_Species<-paste(Qian_Lookup_Angios_TPL$Genus,Qian_Lookup_Angios_TPL$Species,sep='_')
Qian_Lookup_Angios_TPL<-Qian_Lookup_Angios_TPL[,c(c(1:5),23)]
#drop tips in Qian tree (not angios, duplicated, hybrids)
remove_tips<-setdiff(Qian$tip.label,Qian_Lookup_Angios_TPL$Old_Species)
Qian_dropped<-drop.tip(Qian,remove_tips)
Qiandf<-data.frame(Qian_dropped$tip.label)
colnames(Qiandf)<-'Tipname'
#replace tip names in Qian Tree with TPL+taxonlookup alternative
Qian_Lookup_Angios_TPL<-merge(Qiandf,Qian_Lookup_Angios_TPL,by.x='Tipname',by.y='Old_Species')
Qian_Lookup_Angios_TPL$Tipname<-as.character(Qian_Lookup_Angios_TPL$Tipname)
Qian_Lookup_Angios_TPL<-Qian_Lookup_Angios_TPL[match(Qian_dropped$tip.label,Qian_Lookup_Angios_TPL$Tipname),]
Qian_dropped$tip.label<-Qian_Lookup_Angios_TPL$New_Species
#overlap seed and tree datasets
Kew.TaxLookUp<-read.table('./output/tables/KewSIDdata_taxonlookup.txt',header=T,quote='',sep='\t')
Seed_Data<-Kew.TaxLookUp$Log10.Seed.Weight
names(Seed_Data)<-Kew.TaxLookUp$Species
name.check.Seed<-name.check(Qian_dropped,Seed_Data)
Qian_dropped.Seed<-drop.tip(Qian_dropped,name.check.Seed$tree_not_data)
Kew.TaxLookUp.Qiandropped<-merge(data.frame(Qian_dropped.Seed$tip.label),Kew.TaxLookUp,by.x='Qian_dropped.Seed.tip.label',by.y='Species')

tree<-Qian_dropped.Seed
#get total info from TaxonLookUp
plant_lookup_version_current()
PlantLookup<-plant_lookup(include_counts = TRUE)
#get total number of genera in each family
PlantLookup.2<-aggregate(genus~family, data=PlantLookup, length)
PlantLookup<-merge(PlantLookup, PlantLookup.2, by.x = 'family', by.y = 'family', all.x = TRUE)
colnames(PlantLookup)[3]<-'Genus'
colnames(PlantLookup)[6]<-'number.of.genera'
#run TaxonLookup on taxa from the tree
tree.table<-lookup_table(tree$tip.label)
QianPhyndrGeneracount<-aggregate(genus~family, data=tree.table, length)
QianPhyndrSeednooutilers_GenusTree.table<-merge(tree.table,QianPhyndrGeneracount,by.x='family',by.y='family')
colnames(QianPhyndrSeednooutilers_GenusTree.table)[2]<-'Genus'
colnames(QianPhyndrSeednooutilers_GenusTree.table)[5]<-'number.of.genera'
#merge tree info with full taxonlookup info
Phyndr_Genus_Incomplete<-merge(QianPhyndrSeednooutilers_GenusTree.table, PlantLookup, by.x = 'Genus', by.y = 'Genus')
Phyndr_Genus_Incomplete$family.y<-NULL
Phyndr_Genus_Incomplete$order.y<-NULL
Phyndr_Genus_Incomplete$group.y<-NULL
colnames(Phyndr_Genus_Incomplete)[c(2,3,4,5,6,7)]<-c('family','order','group','n.genera.tree','n.species','n.total.genera')
Kew.TaxLookUp.Qiandropped<-merge(Kew.TaxLookUp.Qiandropped,Phyndr_Genus_Incomplete,by.x='genus',by.y='Genus')
Kew.TaxLookUp.Qiandropped$family.y<-NULL
Kew.TaxLookUp.Qiandropped$order.y<-NULL
Kew.TaxLookUp.Qiandropped$group.y<-NULL
colnames(Kew.TaxLookUp.Qiandropped)[c(3,4,5,9)]<-c('family','order','group','n.total.species.genus')
freqs.genera<-table(as.character(Kew.TaxLookUp.Qiandropped$genus))
freqs.genera<-data.frame(names(freqs.genera),unname(freqs.genera))
freqs.genera<-freqs.genera[,-2]
colnames(freqs.genera)<-c('Genus','n.species.tree.genus')
Kew.TaxLookUp.Qiandropped<-merge(Kew.TaxLookUp.Qiandropped,freqs.genera,by.x='genus',by.y='Genus')

genus.mean.seedmass<-aggregate(Log10.Seed.Weight~genus, data=Kew.TaxLookUp.Qiandropped, mean)
colnames(genus.mean.seedmass)[2]<-'genus.mean.Log10.Seed.Weight'
genus.sd.seedmass<-aggregate(Log10.Seed.Weight~genus, data=Kew.TaxLookUp.Qiandropped, sd)
colnames(genus.sd.seedmass)[2]<-'genus.sd.Log10.Seed.Weight'
family.mean.seedmass<-aggregate(Log10.Seed.Weight~family, data=Kew.TaxLookUp.Qiandropped, mean)
colnames(family.mean.seedmass)[2]<-'family.mean.Log10.Seed.Weight'
family.sd.seedmass<-aggregate(Log10.Seed.Weight~family, data=Kew.TaxLookUp.Qiandropped, sd)
colnames(family.sd.seedmass)[2]<-'family.sd.Log10.Seed.Weight'
Kew.TaxLookUp.Qiandropped<-merge(Kew.TaxLookUp.Qiandropped,genus.mean.seedmass,by.x='genus',by.y='genus')
Kew.TaxLookUp.Qiandropped<-merge(Kew.TaxLookUp.Qiandropped,genus.sd.seedmass,by.x='genus',by.y='genus')
Kew.TaxLookUp.Qiandropped<-merge(Kew.TaxLookUp.Qiandropped,family.mean.seedmass,by.x='family',by.y='family')
Kew.TaxLookUp.Qiandropped<-merge(Kew.TaxLookUp.Qiandropped,family.sd.seedmass,by.x='family',by.y='family')
family.species.count<-aggregate(number.of.species~family,PlantLookup,sum)
colnames(family.species.count)[2]<-c('n.total.species.family')
Kew.TaxLookUp.Qiandropped<-merge(Kew.TaxLookUp.Qiandropped,family.species.count)
Kew.TaxLookUp.Qiandropped.genera<-unique(Kew.TaxLookUp.Qiandropped[c(1,2,11,16)])
family.species<-aggregate(n.species.tree.genus~family,Kew.TaxLookUp.Qiandropped.genera,sum)
colnames(family.species)[2]<-'n.species.tree.family'
Kew.TaxLookUp.Qiandropped<-merge(Kew.TaxLookUp.Qiandropped,family.species)
Kew.TaxLookUp.Qiandropped$family.genus.sampling.fraction<-Kew.TaxLookUp.Qiandropped$n.genera.tree/Kew.TaxLookUp.Qiandropped$n.total.genera
Kew.TaxLookUp.Qiandropped$genus.species.sampling.fraction<-Kew.TaxLookUp.Qiandropped$n.species.tree.genus/Kew.TaxLookUp.Qiandropped$n.total.species.genus
Kew.TaxLookUp.Qiandropped$family.species.sampling.fraction<-Kew.TaxLookUp.Qiandropped$n.species.tree.family/Kew.TaxLookUp.Qiandropped$n.total.species.family
Kew.TaxLookUp.Qiandropped[Kew.TaxLookUp.Qiandropped$genus.species.sampling.fraction>1,]$genus.species.sampling.fraction<-1
Kew.TaxLookUp.Qiandropped$genus.CV.Log10.Seed.Weight<-(Kew.TaxLookUp.Qiandropped$genus.sd.Log10.Seed.Weight/Kew.TaxLookUp.Qiandropped$genus.mean.Log10.Seed.Weight)*100
Kew.TaxLookUp.Qiandropped$family.CV.Log10.Seed.Weight<-(Kew.TaxLookUp.Qiandropped$family.sd.Log10.Seed.Weight/Kew.TaxLookUp.Qiandropped$family.mean.Log10.Seed.Weight)*100


