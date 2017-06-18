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
###################################################################
###################################################################
###################################################################
###################################################################
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

######get phylogroups function by Victor Soria-Carrasco
######modified from extract_phylogroups at https://bitbucket.org/visoca/compdivrate/
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


######this runs RPANDA across phylogroups
####run_RPANDA_BM_phylogroup_select_size<-function(phylogroup,table,sampling,mincladesize){
####  genus.table<-table[match(phylogroup$tip.label,table$tip),]
####  cat(genus.table$New_Species[1],'\n')
####  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
####  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
####  genus.table.genera<-unique(genus.table.genera)
####  genus.counts<-as.data.frame(table(as.character(genus.table$genus)))
####  colnames(genus.counts)<-c('genus','n.of.tips')
####  genus.table.genera<-merge(genus.table.genera,genus.counts)
####  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
####  genus.table$phylogroup.sampling.fraction<-sampling.fraction
####  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
####  tot_time<-max(node.age(phylogroup)$ages)
####  genus.table$phylogroup.size<-nrow(genus.table)
####  genus.table$phylogroup.name<-phylogroup$tip.label[1]
####  #skip for small or undersampled clades
####  if(genus.table$phylogroup.size[1]<mincladesize){
####    genus.table$lambda.rpanda1<-NA
####    genus.table$mu.rpanda1<-NA
####    genus.table$ndr.rpanda1<-NA
####    genus.table$sigsq<-NA
####    genus.table$z0<-NA
####    genus.table$aicc<-NA
####    return(genus.table)
####  }
####  if(genus.table$phylogroup.sampling.fraction[1]<sampling){
####    genus.table$lambda.rpanda1<-NA
####    genus.table$mu.rpanda1<-NA
####    genus.table$ndr.rpanda1<-NA
####    genus.table$sigsq<-NA
####    genus.table$z0<-NA
####    genus.table$aicc<-NA
####    return(genus.table)
####  }
####  f.lamb.cst<-function(t,y){y[1]}
####  f.lamb.exp<-function(t,y){y[1]*exp(y[2]*t)}
####  f.mu.cst.0<-function(t,y){0}
####  f.mu.cst<-function(t,y){y[1]}
####  f.mu.exp<-function(t,y){y[1]*exp(y[2]*t);y[1]*exp(y[2]*t)}
####  lamb_par_init.exp<-c(0.5,0.01)
####  mu_par_init.exp<-c(0.05,-0.01)
####  lamb_par_init.cst<-c(0.5)
####  mu_par_init.cst<-c(0.05)
####  mu_par_init.0<-c()
####  res.lambda.cst.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par_init.cst,mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE)
####  if(res.lambda.cst.mu.0$lamb_par[1]<0){
####    cat('negative','\n')
####    res.lambda.cst.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par = c(5),mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE)
####  }
####  res.lambda.exp.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par_init.exp,mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE)
####  if(res.lambda.exp.mu.0$lamb_par[1]<0){
####    cat('negative','\n')
####    res.lambda.exp.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par=c(5,0.05),mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE)
####  }
####  res.lambda.cst.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par_init.cst,mu_par_init.cst,f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE)
####  if(res.lambda.cst.mu.cst$mu_par[1]<0){
####    cat('negative','\n')
####    res.lambda.cst.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par=c(5),mu_par=c(0.01),f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE)
####  }  
####  #res.lambda.cst.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par_init.cst,mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE)
####  #if(res.lambda.cst.mu.exp$mu_par[1]<0){
####  #  cat('negative','\n')
####  #  res.lambda.cst.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par = c(5),mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE)
####  #}  
####  #res.lambda.exp.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.exp,lamb_par_init.exp,mu_par_init.exp,f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE)
####  #if(res.lambda.exp.mu.exp$mu_par[1]<0){
####  #  res.lambda.exp.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.exp,lamb_par=c(5,0.05),mu_par=c(0.05,-0.01),f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE)
####  ##  cat('negative','\n')
####  #}  
####  res.lambda.exp.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par_init.exp,mu_par_init.cst,f=sampling.fraction,expo.lamb=TRUE,cst.mu=TRUE)
####  if(res.lambda.exp.mu.cst$mu_par[1]<0){
####    cat('negative','\n')
####    res.lambda.exp.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par=c(5,0.05),mu_par=c(0.01),f=sampling.fraction,expo.lamb=TRUE,cst.mu=TRUE)
####  }  
####  #aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.exp.mu.0$aicc,res.lambda.exp.mu.exp$aicc,res.lambda.exp.mu.cst$aicc,res.lambda.cst.mu.exp$aicc,res.lambda.cst.mu.cst$aicc)
####  #names(aicc.vector)<-c('lambda.cst.mu0','lambda.exp.mu.0','lambda.exp.mu.exp','lambda.exp.mu.cst','lambda.cst.mu.exp','lambda.cst.mu.cst')
####  aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.exp.mu.0$aicc,res.lambda.exp.mu.cst$aicc,res.lambda.cst.mu.cst$aicc)
####  names(aicc.vector)<-c('lambda.cst.mu0','lambda.exp.mu.0','lambda.exp.mu.cst','lambda.cst.mu.cst')
####  aicc.weights<-sort(Weights(aicc.vector),decreasing=T)
####  if (names(aicc.weights)[1]=='lambda.cst.mu0'){
####    genus.table$lambda.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
####    genus.table$mu.rpanda1<-0
####    genus.table$ndr.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.exp.mu.0'){
####    genus.table$lambda.rpanda1<-res.lambda.exp.mu.0$lamb_par[1]
####    genus.table$mu.rpanda1<-0
####    genus.table$ndr.rpanda1<-res.lambda.exp.mu.0$lamb_par[1]
####  }
####  #else if (names(aicc.weights)[1]=='lambda.exp.mu.exp'){
####  #  genus.table$lambda.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]
####  #  genus.table$mu.rpanda1<-res.lambda.exp.mu.exp$mu_par[1]
####  #  genus.table$ndr.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]-res.lambda.exp.mu.exp$mu_par[1]
####  #}
####  else if (names(aicc.weights)[1]=='lambda.exp.mu.cst'){
####    genus.table$lambda.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.exp.mu.cst$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]-res.lambda.exp.mu.cst$mu_par[1]
####  }
####  #else if (names(aicc.weights)[1]=='lambda.cst.mu.exp'){
####  #  genus.table$lambda.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]
####  #  genus.table$mu.rpanda1<-res.lambda.cst.mu.exp$mu_par[1]
####  #  genus.table$ndr.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]-res.lambda.cst.mu.exp$mu_par[1]
####  #}
####  else if (names(aicc.weights)[1]=='lambda.cst.mu.cst'){
####    genus.table$lambda.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.cst.mu.cst$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]-res.lambda.cst.mu.cst$mu_par[1]
####  }
####  if (genus.table$mu.rpanda1[1]<0){
####    cat('*****NEGATIVE EXTINCTION RATE','\n')
####  }
####  genus.table$aicc<-names(aicc.weights)[1]
####  trait<-as.numeric(genus.table$Log10.Seed.Weight)
####  names(trait)<-genus.table$tip
####  trait<-trait[!is.na(trait)]
####  if(length((trait))>1){
####    res.trait1<-fitContinuous(phylogroup,trait,model="BM")
####    genus.table$sigsq<-res.trait1$opt$sigsq
####    genus.table$z0<-res.trait1$opt$z0
####  }else{
####    genus.table$sigsq<-NA
####    genus.table$z0<-NA
####  }
####  
####  return(genus.table)
####}
##########this runs RPANDA across phylogroups
####run_RPANDA_EB_phylogroup_select_size<-function(phylogroup,table,sampling,mincladesize){
####  genus.table<-table[match(phylogroup$tip.label,table$tip),]
####  cat(genus.table$New_Species[1],'\n')
####  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
####  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
####  genus.table.genera<-unique(genus.table.genera)
####  genus.counts<-as.data.frame(table(as.character(genus.table$genus)))
####  colnames(genus.counts)<-c('genus','n.of.tips')
####  genus.table.genera<-merge(genus.table.genera,genus.counts)
####  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
####  genus.table$phylogroup.sampling.fraction<-sampling.fraction
####  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
####  tot_time<-max(node.age(phylogroup)$ages)
####  genus.table$phylogroup.size<-nrow(genus.table)
####  genus.table$phylogroup.name<-phylogroup$tip.label[1]
####  #skip for small or undersampled clades
####  if(genus.table$phylogroup.size[1]<mincladesize){
####    genus.table$lambda.rpanda1<-NA
####    genus.table$mu.rpanda1<-NA
####    genus.table$ndr.rpanda1<-NA
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####    genus.table$aicc<-NA
####    return(genus.table)
####  }
####  if(genus.table$phylogroup.sampling.fraction[1]<sampling){
####    genus.table$lambda.rpanda1<-NA
####    genus.table$mu.rpanda1<-NA
####    genus.table$ndr.rpanda1<-NA
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####    genus.table$aicc<-NA
####    return(genus.table)
####  }
####  f.lamb.cst<-function(t,y){y[1]}
####  f.lamb.exp<-function(t,y){y[1]*exp(y[2]*t)}
####  f.mu.cst.0<-function(t,y){0}
####  f.mu.cst<-function(t,y){y[1]}
####  f.mu.exp<-function(t,y){y[1]*exp(y[2]*t);y[1]*exp(y[2]*t)}
####  lamb_par_init.exp<-c(0.5,0.01)
####  mu_par_init.exp<-c(0.05,-0.01)
####  lamb_par_init.cst<-c(0.5)
####  mu_par_init.cst<-c(0.05)
####  mu_par_init.0<-c()
####  res.lambda.cst.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par_init.cst,mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE)
####  if(res.lambda.cst.mu.0$lamb_par[1]<0){
####    cat('negative','\n')
####    res.lambda.cst.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par = c(5),mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE)
####  }
####  res.lambda.exp.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par_init.exp,mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE)
####  if(res.lambda.exp.mu.0$lamb_par[1]<0){
####    cat('negative','\n')
####    res.lambda.exp.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par=c(5,0.05),mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE)
####  }
####  res.lambda.cst.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par_init.cst,mu_par_init.cst,f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE)
####  if(res.lambda.cst.mu.cst$mu_par[1]<0){
####    cat('negative','\n')
####    res.lambda.cst.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par=c(5),mu_par=c(0.01),f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE)
####  }  
####  #res.lambda.cst.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par_init.cst,mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE)
####  #if(res.lambda.cst.mu.exp$mu_par[1]<0){
####  #  cat('negative','\n')
####  #  res.lambda.cst.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par = c(5),mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE)
####  #}  
####  #res.lambda.exp.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.exp,lamb_par_init.exp,mu_par_init.exp,f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE)
####  #if(res.lambda.exp.mu.exp$mu_par[1]<0){
####  #  res.lambda.exp.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.exp,lamb_par=c(5,0.05),mu_par=c(0.05,-0.01),f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE)
####  ##  cat('negative','\n')
####  #}  
####  res.lambda.exp.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par_init.exp,mu_par_init.cst,f=sampling.fraction,expo.lamb=TRUE,cst.mu=TRUE)
####  if(res.lambda.exp.mu.cst$mu_par[1]<0){
####    cat('negative','\n')
####    res.lambda.exp.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par=c(5,0.05),mu_par=c(0.01),f=sampling.fraction,expo.lamb=TRUE,cst.mu=TRUE)
####  }  
####  #aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.exp.mu.0$aicc,res.lambda.exp.mu.exp$aicc,res.lambda.exp.mu.cst$aicc,res.lambda.cst.mu.exp$aicc,res.lambda.cst.mu.cst$aicc)
####  #names(aicc.vector)<-c('lambda.cst.mu0','lambda.exp.mu.0','lambda.exp.mu.exp','lambda.exp.mu.cst','lambda.cst.mu.exp','lambda.cst.mu.cst')
####  aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.exp.mu.0$aicc,res.lambda.exp.mu.cst$aicc,res.lambda.cst.mu.cst$aicc)
####  names(aicc.vector)<-c('lambda.cst.mu0','lambda.exp.mu.0','lambda.exp.mu.cst','lambda.cst.mu.cst')
####  aicc.weights<-sort(Weights(aicc.vector),decreasing=T)
####  if (names(aicc.weights)[1]=='lambda.cst.mu0'){
####    genus.table$lambda.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
####    genus.table$mu.rpanda1<-0
####    genus.table$ndr.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.exp.mu.0'){
####    genus.table$lambda.rpanda1<-res.lambda.exp.mu.0$lamb_par[1]
####    genus.table$mu.rpanda1<-0
####    genus.table$ndr.rpanda1<-res.lambda.exp.mu.0$lamb_par[1]
####  }
####  #else if (names(aicc.weights)[1]=='lambda.exp.mu.exp'){
####  #  genus.table$lambda.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]
####  #  genus.table$mu.rpanda1<-res.lambda.exp.mu.exp$mu_par[1]
####  #  genus.table$ndr.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]-res.lambda.exp.mu.exp$mu_par[1]
####  #}
####  else if (names(aicc.weights)[1]=='lambda.exp.mu.cst'){
####    genus.table$lambda.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.exp.mu.cst$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]-res.lambda.exp.mu.cst$mu_par[1]
####  }
####  #else if (names(aicc.weights)[1]=='lambda.cst.mu.exp'){
####  #  genus.table$lambda.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]
####  #  genus.table$mu.rpanda1<-res.lambda.cst.mu.exp$mu_par[1]
####  #  genus.table$ndr.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]-res.lambda.cst.mu.exp$mu_par[1]
####  #}
####  else if (names(aicc.weights)[1]=='lambda.cst.mu.cst'){
####    genus.table$lambda.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.cst.mu.cst$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]-res.lambda.cst.mu.cst$mu_par[1]
####  }
####  if (genus.table$mu.rpanda1[1]<0){
####    cat('*****NEGATIVE EXTINCTION RATE','\n')
####  }
####  genus.table$aicc<-names(aicc.weights)[1]
####  trait<-as.numeric(genus.table$Log10.Seed.Weight)
####  names(trait)<-genus.table$tip
####  trait<-trait[!is.na(trait)]
####  if (length(trait)<2){
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####  }else{
####    name.check<-name.check(phylogroup,trait)
####    if(name.check!='OK'){
####      tips.with.data<-drop.tip(phylogroup,name.check$tree_not_data)
####    }else{
####      tips.with.data<-phylogroup
####    }
####    if(length(tips.with.data$tip.label)>mincladesize){
####      res.trait1.BM<-fitContinuous(phylogroup,trait,model="BM")
####      res.trait1.EB<-fitContinuous(phylogroup,trait,model="EB")
####      genus.table$sigsq.BM<-res.trait1.BM$opt$sigsq
####      genus.table$aicc.BM<-res.trait1.BM$opt$aicc
####      genus.table$sigsq.EB<-res.trait1.EB$opt$sigsq
####      genus.table$aicc.EB<-res.trait1.EB$opt$aicc
####      genus.table$a.EB<-res.trait1.EB$opt$a
####    }else{
####      genus.table$sigsq.BM<-NA
####      genus.table$aicc.BM<-NA
####      genus.table$sigsq.EB<-NA
####      genus.table$aicc.EB<-NA
####      genus.table$a.EB<-NA
####    }
####    
####  }
####  return(genus.table)
####}
####
####run_RPANDA_EB_phylogroup_select_size_6models<-function(phylogroup,table,sampling,mincladesize){
####  genus.table<-table[match(phylogroup$tip.label,table$tip),]
####  cat(genus.table$New_Species[1],'\n')
####  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
####  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
####  genus.table.genera<-unique(genus.table.genera)
####  genus.counts<-as.data.frame(table(as.character(genus.table$genus)))
####  colnames(genus.counts)<-c('genus','n.of.tips')
####  genus.table.genera<-merge(genus.table.genera,genus.counts)
####  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
####  genus.table$phylogroup.sampling.fraction<-sampling.fraction
####  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
####  tot_time<-max(node.age(phylogroup)$ages)
####  genus.table$phylogroup.size<-nrow(genus.table)
####  genus.table$phylogroup.name<-phylogroup$tip.label[1]
####  #skip for small or undersampled clades
####  if(genus.table$phylogroup.size[1]<mincladesize){
####    genus.table$lambda.rpanda1<-NA
####    genus.table$mu.rpanda1<-NA
####    genus.table$ndr.rpanda1<-NA
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####    genus.table$aicc<-NA
####    return(genus.table)
####  }
####  if(genus.table$phylogroup.sampling.fraction[1]<sampling){
####    genus.table$lambda.rpanda1<-NA
####    genus.table$mu.rpanda1<-NA
####    genus.table$ndr.rpanda1<-NA
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####    genus.table$aicc<-NA
####    return(genus.table)
####  }
####  f.lamb.cst<-function(t,y){y[1]}
####  f.lamb.exp<-function(t,y){y[1]*exp(y[2]*t)}
####  f.mu.cst.0<-function(t,y){0}
####  f.mu.cst<-function(t,y){y[1]}
####  f.mu.exp<-function(t,y){y[1]*exp(y[2]*t)}
####  lamb_par_init.exp<-c(0.5,0.01)
####  mu_par_init.exp<-c(0.05,-0.01)
####  lamb_par_init.cst<-c(0.5)
####  mu_par_init.cst<-c(0.05)
####  mu_par_init.0<-c()
####  res.lambda.cst.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par_init.cst,mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE)
####  if(res.lambda.cst.mu.0$lamb_par[1]<0){
####    cat('negative','\n')
####    res.lambda.cst.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par = c(5),mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE)
####  }
####  if(res.lambda.cst.mu.0$lamb_par[1]>3){
####    lamb_par_init.exp<-c(5,0.01)
####    mu_par_init.exp<-c(5,0.01)
####    lamb_par_init.cst<-c(0.5)
####    mu_par_init.cst<-c(0.2)
####    mu_par_init.0<-c()
####    
####  }
####  cat('lambdaexpmu0','\n')
####  res.lambda.exp.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par_init.exp,mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE)
####  if(res.lambda.exp.mu.0$lamb_par[1]<0){
####    cat('negative','\n')
####    res.lambda.exp.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par=c(5,0.05),mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE)
####  }
####  cat('lambdacstmucst','\n')
####  res.lambda.cst.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par_init.cst,mu_par_init.cst,f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE)
####  if(res.lambda.cst.mu.cst$mu_par[1]<0){
####    cat('negative','\n')
####    res.lambda.cst.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par=c(5),mu_par=c(0.01),f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE)
####  }
####  cat('lambdacstmuexp','\n')
####  res.lambda.cst.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par_init.cst,mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE)
####  if(res.lambda.cst.mu.exp$mu_par[1]<0){
####    cat('negative','\n')
####    res.lambda.cst.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par = c(5),mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE)
####  }  
####  cat('lambdaexpmuexp','\n')
####  res.lambda.exp.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.exp,lamb_par_init.exp,mu_par_init.exp,f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE)
####  if(res.lambda.exp.mu.exp$mu_par[1]<0){
####    res.lambda.exp.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.exp,lamb_par=c(5,0.05),mu_par=c(0.05,-0.01),f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE)
####  #  cat('negative','\n')
####  }
####  cat('lambdaexpmucst','\n')
####  res.lambda.exp.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par_init.exp,mu_par_init.cst,f=sampling.fraction,expo.lamb=TRUE,cst.mu=TRUE)
####  if(res.lambda.exp.mu.cst$mu_par[1]<0){
####    cat('negative','\n')
####    res.lambda.exp.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par=c(5,0.05),mu_par=c(0.01),f=sampling.fraction,expo.lamb=TRUE,cst.mu=TRUE)
####  }  
####  aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.exp.mu.0$aicc,res.lambda.exp.mu.exp$aicc,res.lambda.exp.mu.cst$aicc,res.lambda.cst.mu.exp$aicc,res.lambda.cst.mu.cst$aicc)
####  names(aicc.vector)<-c('lambda.cst.mu0','lambda.exp.mu.0','lambda.exp.mu.exp','lambda.exp.mu.cst','lambda.cst.mu.exp','lambda.cst.mu.cst')
####  #aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.exp.mu.0$aicc,res.lambda.exp.mu.cst$aicc,res.lambda.cst.mu.cst$aicc)
####  #names(aicc.vector)<-c('lambda.cst.mu0','lambda.exp.mu.0','lambda.exp.mu.cst','lambda.cst.mu.cst')
####  aicc.weights<-sort(Weights(aicc.vector),decreasing=T)
####  if (names(aicc.weights)[1]=='lambda.cst.mu0'){
####    genus.table$lambda.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
####    genus.table$mu.rpanda1<-0
####    genus.table$ndr.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.exp.mu.0'){
####    genus.table$lambda.rpanda1<-res.lambda.exp.mu.0$lamb_par[1]
####    genus.table$mu.rpanda1<-0
####    genus.table$ndr.rpanda1<-res.lambda.exp.mu.0$lamb_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.exp.mu.exp'){
####    genus.table$lambda.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.exp.mu.exp$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]-res.lambda.exp.mu.exp$mu_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.exp.mu.cst'){
####    genus.table$lambda.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.exp.mu.cst$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]-res.lambda.exp.mu.cst$mu_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.cst.mu.exp'){
####    genus.table$lambda.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.cst.mu.exp$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]-res.lambda.cst.mu.exp$mu_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.cst.mu.cst'){
####    genus.table$lambda.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.cst.mu.cst$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]-res.lambda.cst.mu.cst$mu_par[1]
####  }
####  if (genus.table$mu.rpanda1[1]<0){
####    cat('*****NEGATIVE EXTINCTION RATE','\n')
####  }
####  genus.table$aicc<-names(aicc.weights)[1]
####  trait<-as.numeric(genus.table$Log10.Seed.Weight)
####  names(trait)<-genus.table$tip
####  trait<-trait[!is.na(trait)]
####  if (length(trait)<2){
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####  }else{
####    name.check<-name.check(phylogroup,trait)
####    if(name.check!='OK'){
####      tips.with.data<-drop.tip(phylogroup,name.check$tree_not_data)
####    }else{
####      tips.with.data<-phylogroup
####    }
####    if(length(tips.with.data$tip.label)>mincladesize){
####      res.trait1.BM<-fitContinuous(phylogroup,trait,model="BM")
####      res.trait1.EB<-fitContinuous(phylogroup,trait,model="EB")
####      genus.table$sigsq.BM<-res.trait1.BM$opt$sigsq
####      genus.table$aicc.BM<-res.trait1.BM$opt$aicc
####      genus.table$sigsq.EB<-res.trait1.EB$opt$sigsq
####      genus.table$aicc.EB<-res.trait1.EB$opt$aicc
####      genus.table$a.EB<-res.trait1.EB$opt$a
####    }else{
####      genus.table$sigsq.BM<-NA
####      genus.table$aicc.BM<-NA
####      genus.table$sigsq.EB<-NA
####      genus.table$aicc.EB<-NA
####      genus.table$a.EB<-NA
####    }
####    
####  }
####  return(genus.table)
####}
####
run_RPANDA_EB_z0_phylogroup_select_size_6models<-function(phylogroup,table,sampling,mincladesize){
  genus.table<-table[match(phylogroup$tip.label,table$tip),]
  cat(genus.table$New_Species[1],'\n')
  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
  genus.table.genera<-unique(genus.table.genera)
  genus.counts<-as.data.frame(table(as.character(genus.table$genus)))
  colnames(genus.counts)<-c('genus','n.of.tips')
  genus.table.genera<-merge(genus.table.genera,genus.counts)
  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
  genus.table$phylogroup.sampling.fraction<-sampling.fraction
  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
  tot_time<-max(node.age(phylogroup)$ages)
  genus.table$phylogroup.size<-nrow(genus.table)
  genus.table$phylogroup.name<-phylogroup$tip.label[1]
  #skip for small or undersampled clades
  if(genus.table$phylogroup.size[1]<mincladesize){
    genus.table$lambda.rpanda1<-NA
    genus.table$mu.rpanda1<-NA
    genus.table$ndr.rpanda1<-NA
    genus.table$sigsq.BM<-NA
    genus.table$aicc.BM<-NA
    genus.table$sigsq.EB<-NA
    genus.table$aicc.EB<-NA
    genus.table$a.EB<-NA
    genus.table$aicc<-NA
    genus.table$z0.EB<-NA
    genus.table$z0.BM<-NA
    
    return(genus.table)
  }
  if(genus.table$phylogroup.sampling.fraction[1]<sampling){
    genus.table$lambda.rpanda1<-NA
    genus.table$mu.rpanda1<-NA
    genus.table$ndr.rpanda1<-NA
    genus.table$sigsq.BM<-NA
    genus.table$aicc.BM<-NA
    genus.table$sigsq.EB<-NA
    genus.table$aicc.EB<-NA
    genus.table$a.EB<-NA
    genus.table$aicc<-NA
    genus.table$z0.EB<-NA
    genus.table$z0.BM<-NA
    return(genus.table)
  }
  f.lamb.cst<-function(t,y){y[1]}
  f.lamb.exp<-function(t,y){y[1]*exp(y[2]*t)}
  f.mu.cst.0<-function(t,y){0}
  f.mu.cst<-function(t,y){y[1]}
  f.mu.exp<-function(t,y){y[1]*exp(y[2]*t)}
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
  if(res.lambda.cst.mu.0$lamb_par[1]>3){
    lamb_par_init.exp<-c(5,0.01)
    mu_par_init.exp<-c(5,0.01)
    lamb_par_init.cst<-c(0.5)
    mu_par_init.cst<-c(0.2)
    mu_par_init.0<-c()
    
  }
  cat('lambdaexpmu0','\n')
  res.lambda.exp.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par_init.exp,mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE)
  if(res.lambda.exp.mu.0$lamb_par[1]<0){
    cat('negative','\n')
    res.lambda.exp.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par=c(5,0.05),mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE)
  }
  cat('lambdacstmucst','\n')
  res.lambda.cst.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par_init.cst,mu_par_init.cst,f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE)
  if(res.lambda.cst.mu.cst$mu_par[1]<0){
    cat('negative','\n')
    res.lambda.cst.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par=c(5),mu_par=c(0.01),f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE)
  }
  cat('lambdacstmuexp','\n')
  res.lambda.cst.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par_init.cst,mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE)
  if(res.lambda.cst.mu.exp$mu_par[1]<0){
    cat('negative','\n')
    res.lambda.cst.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par = c(5),mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE)
  }  
  cat('lambdaexpmuexp','\n')
  res.lambda.exp.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.exp,lamb_par_init.exp,mu_par_init.exp,f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE)
  if(res.lambda.exp.mu.exp$mu_par[1]<0){
    res.lambda.exp.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.exp,lamb_par=c(5,0.05),mu_par=c(0.05,-0.01),f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE)
    #  cat('negative','\n')
  }
  cat('lambdaexpmucst','\n')
  res.lambda.exp.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par_init.exp,mu_par_init.cst,f=sampling.fraction,expo.lamb=TRUE,cst.mu=TRUE)
  if(res.lambda.exp.mu.cst$mu_par[1]<0){
    cat('negative','\n')
    res.lambda.exp.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par=c(5,0.05),mu_par=c(0.01),f=sampling.fraction,expo.lamb=TRUE,cst.mu=TRUE)
  }  
  aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.exp.mu.0$aicc,res.lambda.exp.mu.exp$aicc,res.lambda.exp.mu.cst$aicc,res.lambda.cst.mu.exp$aicc,res.lambda.cst.mu.cst$aicc)
  names(aicc.vector)<-c('lambda.cst.mu0','lambda.exp.mu.0','lambda.exp.mu.exp','lambda.exp.mu.cst','lambda.cst.mu.exp','lambda.cst.mu.cst')
  #aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.exp.mu.0$aicc,res.lambda.exp.mu.cst$aicc,res.lambda.cst.mu.cst$aicc)
  #names(aicc.vector)<-c('lambda.cst.mu0','lambda.exp.mu.0','lambda.exp.mu.cst','lambda.cst.mu.cst')
  aicc.weights<-sort(Weights(aicc.vector),decreasing=T)
  if (names(aicc.weights)[1]=='lambda.cst.mu0'){
    genus.table$lambda.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
    genus.table$mu.rpanda1<-0
    genus.table$ndr.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
  }
  else if (names(aicc.weights)[1]=='lambda.exp.mu.0'){
    genus.table$lambda.rpanda1<-res.lambda.exp.mu.0$lamb_par[1]
    genus.table$mu.rpanda1<-0
    genus.table$ndr.rpanda1<-res.lambda.exp.mu.0$lamb_par[1]
  }
  else if (names(aicc.weights)[1]=='lambda.exp.mu.exp'){
    genus.table$lambda.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]
    genus.table$mu.rpanda1<-res.lambda.exp.mu.exp$mu_par[1]
    genus.table$ndr.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]-res.lambda.exp.mu.exp$mu_par[1]
  }
  else if (names(aicc.weights)[1]=='lambda.exp.mu.cst'){
    genus.table$lambda.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]
    genus.table$mu.rpanda1<-res.lambda.exp.mu.cst$mu_par[1]
    genus.table$ndr.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]-res.lambda.exp.mu.cst$mu_par[1]
  }
  else if (names(aicc.weights)[1]=='lambda.cst.mu.exp'){
    genus.table$lambda.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]
    genus.table$mu.rpanda1<-res.lambda.cst.mu.exp$mu_par[1]
    genus.table$ndr.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]-res.lambda.cst.mu.exp$mu_par[1]
  }
  else if (names(aicc.weights)[1]=='lambda.cst.mu.cst'){
    genus.table$lambda.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]
    genus.table$mu.rpanda1<-res.lambda.cst.mu.cst$mu_par[1]
    genus.table$ndr.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]-res.lambda.cst.mu.cst$mu_par[1]
  }
  if (genus.table$mu.rpanda1[1]<0){
    cat('*****NEGATIVE EXTINCTION RATE','\n')
  }
  genus.table$aicc<-names(aicc.weights)[1]
  trait<-as.numeric(genus.table$Log10.Seed.Weight)
  names(trait)<-genus.table$tip
  trait<-trait[!is.na(trait)]
  if (length(trait)<2){
    genus.table$sigsq.BM<-NA
    genus.table$aicc.BM<-NA
    genus.table$sigsq.EB<-NA
    genus.table$aicc.EB<-NA
    genus.table$a.EB<-NA
    genus.table$z0.EB<-NA
    genus.table$z0.BM<-NA
  }else{
    name.check<-name.check(phylogroup,trait)
    if(name.check!='OK'){
      tips.with.data<-drop.tip(phylogroup,name.check$tree_not_data)
    }else{
      tips.with.data<-phylogroup
    }
    if(length(tips.with.data$tip.label)>mincladesize){
      res.trait1.BM<-fitContinuous(phylogroup,trait,model="BM")
      res.trait1.EB<-fitContinuous(phylogroup,trait,model="EB")
      genus.table$sigsq.BM<-res.trait1.BM$opt$sigsq
      genus.table$aicc.BM<-res.trait1.BM$opt$aicc
      genus.table$sigsq.EB<-res.trait1.EB$opt$sigsq
      genus.table$aicc.EB<-res.trait1.EB$opt$aicc
      genus.table$a.EB<-res.trait1.EB$opt$a
      genus.table$z0.EB<-res.trait1.EB$opt$z0
      genus.table$z0.BM<-res.trait1.BM$opt$z0
    }else{
      genus.table$sigsq.BM<-NA
      genus.table$aicc.BM<-NA
      genus.table$sigsq.EB<-NA
      genus.table$aicc.EB<-NA
      genus.table$a.EB<-NA
      genus.table$z0.EB<-NA
      genus.table$z0.BM<-NA
    }
    
  }
  return(genus.table)
}

##########this runs MS across phylogroups
####run_MS_BM_phylogroup_select_size<-function(phylogroup,table,sampling,mincladesize){
####  genus.table<-table[match(phylogroup$tip.label,table$tip),]
####  cat(genus.table$New_Species[1],'\n')
####  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
####  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
####  genus.table.genera<-unique(genus.table.genera)
####  genus.counts<-as.data.frame(table(as.character(genus.table$genus)))
####  colnames(genus.counts)<-c('genus','n.of.tips')
####  genus.table.genera<-merge(genus.table.genera,genus.counts)
####  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
####  genus.table$phylogroup.sampling.fraction<-sampling.fraction
####  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
####  tot_time<-max(node.age(phylogroup)$ages)
####  genus.table$phylogroup.size<-nrow(genus.table)
####  genus.table$phylogroup.name<-phylogroup$tip.label[1]
####  #skip for small or undersampled clades
####  if(genus.table$phylogroup.size[1]<mincladesize){
####    genus.table$lambda.ms0<-NA
####    genus.table$lambda.ms05<-NA
####    genus.table$lambda.ms09<-NA
####    genus.table$ndr.ms05<-NA
####    genus.table$ndr.ms09<-NA
####    genus.table$sigsq<-NA
####    genus.table$z0<-NA
####    return(genus.table)
####  }
####  if(genus.table$phylogroup.sampling.fraction[1]<sampling){
####    genus.table$lambda.ms0<-NA
####    genus.table$lambda.ms05<-NA
####    genus.table$lambda.ms09<-NA
####    genus.table$ndr.ms05<-NA
####    genus.table$ndr.ms09<-NA
####    genus.table$sigsq<-NA
####    genus.table$z0<-NA
####    return(genus.table)
####  }
####  genus.table$lambda.ms0<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0)
####  genus.table$lambda.ms05<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.5)/(1-0.5)
####  genus.table$lambda.ms09<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.9)/(1-0.9)
####  genus.table$ndr.ms05<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.5)
####  genus.table$ndr.ms09<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.9)
####  trait<-as.numeric(genus.table$Log10.Seed.Weight)
####  names(trait)<-genus.table$tip
####  trait<-trait[!is.na(trait)]
####  if(length((trait))>1){
####    res.trait1<-fitContinuous(phylogroup,trait,model="BM")
####    genus.table$sigsq<-res.trait1$opt$sigsq
####    genus.table$z0<-res.trait1$opt$z0
####  }else{
####    genus.table$sigsq<-NA
####    genus.table$z0<-NA
####  }
####  return(genus.table)
####}
####
##########this runs MS + trait with BM & EB across phylogroups
####run_MS_EB_phylogroup_select_size<-function(phylogroup,table,sampling,mincladesize){
####  genus.table<-table[match(phylogroup$tip.label,table$tip),]
####  cat(genus.table$New_Species[1],'\n')
####  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
####  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
####  genus.table.genera<-unique(genus.table.genera)
####  genus.counts<-as.data.frame(table(as.character(genus.table$genus)))
####  colnames(genus.counts)<-c('genus','n.of.tips')
####  genus.table.genera<-merge(genus.table.genera,genus.counts)
####  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
####  genus.table$phylogroup.sampling.fraction<-sampling.fraction
####  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
####  tot_time<-max(node.age(phylogroup)$ages)
####  genus.table$phylogroup.size<-nrow(genus.table)
####  genus.table$phylogroup.name<-phylogroup$tip.label[1]
####  #skip for small or undersampled clades
####  if(genus.table$phylogroup.size[1]<mincladesize){
####    genus.table$lambda.ms0<-NA
####    genus.table$lambda.ms05<-NA
####    genus.table$lambda.ms09<-NA
####    genus.table$ndr.ms05<-NA
####    genus.table$ndr.ms09<-NA
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####    return(genus.table)
####  }
####  if(genus.table$phylogroup.sampling.fraction[1]<sampling){
####    genus.table$lambda.ms0<-NA
####    genus.table$lambda.ms05<-NA
####    genus.table$lambda.ms09<-NA
####    genus.table$ndr.ms05<-NA
####    genus.table$ndr.ms09<-NA
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####    return(genus.table)
####  }
####  genus.table$lambda.ms0<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0)
####  genus.table$lambda.ms05<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.5)/(1-0.5)
####  genus.table$lambda.ms09<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.9)/(1-0.9)
####  genus.table$ndr.ms05<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.5)
####  genus.table$ndr.ms09<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.9)
####  trait<-as.numeric(genus.table$Log10.Seed.Weight)
####  names(trait)<-genus.table$tip
####  trait<-trait[!is.na(trait)]
####  if (length(trait)<2){
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####  }else{
####    name.check<-name.check(phylogroup,trait)
####    if(name.check!='OK'){
####      tips.with.data<-drop.tip(phylogroup,name.check$tree_not_data)
####    }else{
####      tips.with.data<-phylogroup
####    }
####    if(length(tips.with.data$tip.label)>mincladesize){
####      res.trait1.BM<-fitContinuous(phylogroup,trait,model="BM")
####      res.trait1.EB<-fitContinuous(phylogroup,trait,model="EB")
####      genus.table$sigsq.BM<-res.trait1.BM$opt$sigsq
####      genus.table$aicc.BM<-res.trait1.BM$opt$aicc
####      genus.table$sigsq.EB<-res.trait1.EB$opt$sigsq
####      genus.table$aicc.EB<-res.trait1.EB$opt$aicc
####      genus.table$a.EB<-res.trait1.EB$opt$a
####    }else{
####      genus.table$sigsq.BM<-NA
####      genus.table$aicc.BM<-NA
####      genus.table$sigsq.EB<-NA
####      genus.table$aicc.EB<-NA
####      genus.table$a.EB<-NA
####    }
####      
####  }
####  return(genus.table)
####}

######this runs MS + trait with BM & EB across phylogroups
run_MS_EB_z0_phylogroup_select_size<-function(phylogroup,table,sampling,mincladesize){
  genus.table<-table[match(phylogroup$tip.label,table$tip),]
  cat(genus.table$New_Species[1],'\n')
  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
  genus.table.genera<-unique(genus.table.genera)
  genus.counts<-as.data.frame(table(as.character(genus.table$genus)))
  colnames(genus.counts)<-c('genus','n.of.tips')
  genus.table.genera<-merge(genus.table.genera,genus.counts)
  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
  genus.table$phylogroup.sampling.fraction<-sampling.fraction
  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
  tot_time<-max(node.age(phylogroup)$ages)
  genus.table$phylogroup.size<-nrow(genus.table)
  genus.table$phylogroup.name<-phylogroup$tip.label[1]
  #skip for small or undersampled clades
  if(genus.table$phylogroup.size[1]<mincladesize){
    genus.table$lambda.ms0<-NA
    genus.table$lambda.ms05<-NA
    genus.table$lambda.ms09<-NA
    genus.table$ndr.ms05<-NA
    genus.table$ndr.ms09<-NA
    genus.table$sigsq.BM<-NA
    genus.table$aicc.BM<-NA
    genus.table$sigsq.EB<-NA
    genus.table$aicc.EB<-NA
    genus.table$a.EB<-NA
    genus.table$z0.EB<-NA
    genus.table$z0.BM<-NA
    return(genus.table)
  }
  if(genus.table$phylogroup.sampling.fraction[1]<sampling){
    genus.table$lambda.ms0<-NA
    genus.table$lambda.ms05<-NA
    genus.table$lambda.ms09<-NA
    genus.table$ndr.ms05<-NA
    genus.table$ndr.ms09<-NA
    genus.table$sigsq.BM<-NA
    genus.table$aicc.BM<-NA
    genus.table$sigsq.EB<-NA
    genus.table$aicc.EB<-NA
    genus.table$a.EB<-NA
    genus.table$z0.EB<-NA
    genus.table$z0.BM<-NA
    return(genus.table)
  }
  genus.table$lambda.ms0<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0)
  genus.table$lambda.ms05<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.5)/(1-0.5)
  genus.table$lambda.ms09<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.9)/(1-0.9)
  genus.table$ndr.ms05<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.5)
  genus.table$ndr.ms09<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.9)
  trait<-as.numeric(genus.table$Log10.Seed.Weight)
  names(trait)<-genus.table$tip
  trait<-trait[!is.na(trait)]
  if (length(trait)<2){
    genus.table$sigsq.BM<-NA
    genus.table$aicc.BM<-NA
    genus.table$sigsq.EB<-NA
    genus.table$aicc.EB<-NA
    genus.table$a.EB<-NA
    genus.table$z0.EB<-NA
    genus.table$z0.BM<-NA
  }else{
    name.check<-name.check(phylogroup,trait)
    if(name.check!='OK'){
      tips.with.data<-drop.tip(phylogroup,name.check$tree_not_data)
    }else{
      tips.with.data<-phylogroup
    }
    if(length(tips.with.data$tip.label)>mincladesize){
      res.trait1.BM<-fitContinuous(phylogroup,trait,model="BM")
      res.trait1.EB<-fitContinuous(phylogroup,trait,model="EB")
      genus.table$sigsq.BM<-res.trait1.BM$opt$sigsq
      genus.table$aicc.BM<-res.trait1.BM$opt$aicc
      genus.table$sigsq.EB<-res.trait1.EB$opt$sigsq
      genus.table$aicc.EB<-res.trait1.EB$opt$aicc
      genus.table$a.EB<-res.trait1.EB$opt$a
      genus.table$z0.EB<-res.trait1.EB$opt$z0
      genus.table$z0.BM<-res.trait1.BM$opt$z0
      
    }else{
      genus.table$sigsq.BM<-NA
      genus.table$aicc.BM<-NA
      genus.table$sigsq.EB<-NA
      genus.table$aicc.EB<-NA
      genus.table$a.EB<-NA
      genus.table$z0.EB<-NA
      genus.table$z0.BM<-NA
      
    }
    
  }
  return(genus.table)
}

###wrapper to run clade analysis (with Magallon & Sanderson estimator) + pgls
#this creates the phylogroups (there's another script to run if phyogroups have previously been saved)
####run_clades_MS_pgls<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=mincladesize,ncores=ncores)
####  #save phylogroup object to Rsave file
####  save(phylogroups,file=paste('phylogroups_',minage,'_',maxage,'.Rsave',sep=''))
####  
####  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####  length(unique(phylogroups.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.MS<-lapply(phylogroups$phylogroups,function(x) run_MS_BM_phylogroup_select_size(x,table,mincladesize,sampling))
####  df.phylogroups <- do.call("rbind", phylogroups.MS)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.ms value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.ms0),]
####  
####  #count the number of seed size data points per phylogroup
####  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
####  }
####  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
####  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
####  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
####  
####  #get the mean seed size in each phylogroup
####  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
####  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  ###############this is for lambda
####  pdf(paste('./output/plots/clade_analyses_new/MSlambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$lambda.ms0,df.phylogroups$lambda.ms05,df.phylogroups$lambda.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.ms0','lambda.ms05','lambda.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedlambda<-try(pgls(log10(lambda.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  pgls.seedsigsq<-try(pgls(log10(lambda.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedlambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedlambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedlambda<-gls(log10(lambda.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(lambda.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSlambda_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  ###############this is for ndr
####  pdf(paste('./output/plots/clade_analyses_new/MSndr_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$ndr.ms0,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.ms0','ndr.ms05','ndr.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedndr<-try(pgls(log10(ndr.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedndr)
####  pgls.seedsigsq<-try(pgls(log10(ndr.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedndr)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedndr)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedndr<-gls(log10(ndr.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(ndr.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSndr_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####}
#this will load presaved phylogroups
####run_clades_MS_pgls_savedphylogroups<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  #phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=mincladesize,ncores=ncores)
####  load(file=paste('phylogroups_',minage,'_',maxage,'.Rsave',sep=''))
####  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####  length(unique(phylogroups.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.MS<-lapply(phylogroups$phylogroups,function(x) run_MS_BM_phylogroup_select_size(x,table,mincladesize=mincladesize,sampling=sampling))
####  df.phylogroups <- do.call("rbind", phylogroups.MS)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.ms value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.ms0),]
####  
####  #count the number of seed size data points per phylogroup
####  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
####  }
####  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
####  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
####  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
####  
####  #get the mean seed size in each phylogroup
####  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
####  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  ###############this is for lambda
####  pdf(paste('./output/plots/clade_analyses_new/MSlambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$lambda.ms0,df.phylogroups$lambda.ms05,df.phylogroups$lambda.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.ms0','lambda.ms05','lambda.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedlambda<-try(pgls(log10(lambda.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  pgls.seedsigsq<-try(pgls(log10(lambda.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedlambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedlambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedlambda<-gls(log10(lambda.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(lambda.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSlambda_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  ###############this is for ndr
####  pdf(paste('./output/plots/clade_analyses_new/MSndr_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.ms05','ndr.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedndr<-try(pgls(log10(ndr.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedndr)
####  pgls.seedsigsq<-try(pgls(log10(ndr.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedndr)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedndr)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedndr<-gls(log10(ndr.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(ndr.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSndr_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####}
####
#####this will load presaved phylogroups
####run_clades_MS_EB_pgls_savedphylogroups<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  #phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=mincladesize,ncores=ncores)
####  load(file=paste('phylogroups_',minage,'_',maxage,'.Rsave',sep=''))
####  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####  length(unique(phylogroups.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.MS<-lapply(phylogroups$phylogroups,function(x) run_MS_EB_phylogroup_select_size(x,table,mincladesize=mincladesize,sampling=sampling))
####  df.phylogroups <- do.call("rbind", phylogroups.MS)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.ms value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.ms0),]
####  
####  #count the number of seed size data points per phylogroup
####  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
####  }
####  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
####  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
####  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
####  
####  #get the mean seed size in each phylogroup
####  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
####  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  #AIC model selection,AICcBM-AICcEB
####  df.phylogroups$diff.AICc.BM.EB<-df.phylogroups$aicc.BM-df.phylogroups$aicc.EB
####  df.phylogroups$best.model.sigsq<-NA
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'best.model.sigsq']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'sigsq.BM']
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'best.model.sigsq']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'sigsq.EB']
####  df.phylogroups$output.best.model<-NA
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'output.best.model']<-'BM'
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'output.best.model']<-'EB'
####  output.best.model<-as.data.frame(table(df.phylogroups$output.best.model))
####  write.table(output.best.model,file=paste('./output/clade_analyses/MS_EB/MS_EB_lambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'best_fit_BMEB_model.txt',sep=''),quote=F,sep='\t',row.names=F)
####  
####  ###############this is for lambda
####  pdf(paste('./output/clade_analyses/MS_EB/MS_EB_lambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$lambda.ms0,df.phylogroups$lambda.ms05,df.phylogroups$lambda.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$best.model.sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.ms0','lambda.ms05','lambda.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedlambda<-try(pgls(log10(lambda.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  pgls.seedsigsq<-try(pgls(log10(lambda.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedlambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedlambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedlambda<-gls(log10(lambda.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(lambda.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/clade_analyses/MS_EB/MS_EB_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_lambda_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  ###############this is for ndr
####  pdf(paste('./output/clade_analyses/MS_EB/MS_EB_ndr_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$best.model.sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.ms05','ndr.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedndr<-try(pgls(log10(ndr.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedndr)
####  pgls.seedsigsq<-try(pgls(log10(ndr.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedndr)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedndr)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedndr<-gls(log10(ndr.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(ndr.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/clade_analyses/MS_EB/MS_EB_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_ndr_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####}



#this will load presaved phylogroups
run_clades_MS_EB_z0_pgls_savedphylogroups<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
  colnames(table)[4]<-'tip'
  colnames(table)[9]<-'genus.species.sampling.fraction'
  cat('getting clades','\n')
  #phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=mincladesize,ncores=ncores)
  load(file=paste('phylogroups_',minage,'_',maxage,'.Rsave',sep=''))
  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
  length(unique(phylogroups.species))
  cat('measuring diversification/trait evolution','\n')
  phylogroups.MS<-lapply(phylogroups$phylogroups,function(x) run_MS_EB_z0_phylogroup_select_size(x,table,mincladesize=mincladesize,sampling=sampling))
  df.phylogroups <- do.call("rbind", phylogroups.MS)
  df.phylogroups<-unique(df.phylogroups)
  #remove rows with no lambda.ms value
  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.ms0),]
  
  #count the number of seed size data points per phylogroup
  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
  }
  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
  
  #get the mean seed size in each phylogroup
  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
  #getting vector of representatives in each phylogroup
  #then keep just one tip per phylogroup
  phylogroup.clades.species<-vector('character')
  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
  }
  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
  colnames(names.phylogroups.tree)<-'Tip'
  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
  
  #another filter to check that size of clade is correct
  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
  #add another filter to select well sampled clades
  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
  #AIC model selection,AICcBM-AICcEB
  
  df.phylogroups$diff.AICc.BM.EB<-df.phylogroups$aicc.BM-df.phylogroups$aicc.EB
  df.phylogroups$best.model.sigsq<-NA
  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'best.model.sigsq']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'sigsq.BM']
  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'best.model.sigsq']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'sigsq.EB']
  df.phylogroups$output.best.model<-NA
  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'output.best.model']<-'BM'
  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'output.best.model']<-'EB'
  output.best.model<-as.data.frame(table(df.phylogroups$output.best.model))
  df.phylogroups$z0<-NA
  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'z0']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'z0.BM']
  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'z0']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'z0.EB']
  correlation.z0.SeedWeight<-cor.test(df.phylogroups$z0,df.phylogroups$phylogroup.Log10.Seed.Weight) 
  correlation.df<-cbind(correlation.z0.SeedWeight$estimate,correlation.z0.SeedWeight$p.value)
  colnames(correlation.df)<-c('estimate','p.value')
  write.table(correlation.df,file=paste('./output/clade_analyses/MS_EB_z0/MS_EB_z0_lambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'cor_z0_SeedWeight.txt',sep=''),quote=F,sep='\t',row.names=F)
  write.table(output.best.model,file=paste('./output/clade_analyses/MS_EB_z0/MS_EB_z0_lambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'best_fit_BMEB_model.txt',sep=''),quote=F,sep='\t',row.names=F)
  
  ###############this is for lambda
  pdf(paste('./output/clade_analyses/MS_EB_z0/MS_EB_z0_lambda_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
  par(mfrow=c(2,2))
  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$z0,df.phylogroups$lambda.ms0,df.phylogroups$lambda.ms05,df.phylogroups$lambda.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$best.model.sigsq)
  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.ms0','lambda.ms05','lambda.ms09','ndr.ms05','ndr.ms09','sigsq')
  phylogroups.tree$node.label<-NULL
  #this is for caper
  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
  #for mean family seed weight
  cat('running pgls','\n')
  pgls.seedlambda<-try(pgls(log10(lambda.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
  #summary(pgls.seedlambda)
  pgls.seedsigsq<-try(pgls(log10(lambda.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
  #summary(pgls.seedsigsq)
  ####try both optimisations (caper + nlme if caper gives an error)
  if(class(pgls.seedlambda)!='try-error'){
    plot(log10(phylogroups.tree.object$data$lambda.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
    abline(pgls.seedlambda)
    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
  }
  if(class(pgls.seedlambda)=='try-error'){
    row.names(phylogroups.data)<-phylogroups.data$Tip
    pgls.seedlambda<-gls(log10(lambda.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
    plot(log10(phylogroups.data$lambda.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
    abline(pgls.seedlambda)
    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
  }
  if(class(pgls.seedsigsq)!='try-error'){
    plot(log10(phylogroups.tree.object$data$lambda.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
    abline(pgls.seedsigsq)
    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
  }
  if(class(pgls.seedsigsq)=='try-error'){
    row.names(phylogroups.data)<-phylogroups.data$Tip
    pgls.seedsigsq<-gls(log10(lambda.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
    plot(log10(phylogroups.data$lambda.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
    abline(pgls.seedsigsq)
    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
  }
  
  dev.off()
  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
  colnames(results.df)<-c('analysis','slope','pvalue')
  rownames(results.df)<-NULL
  write.table(results.df,file=paste('./output/clade_analyses/MS_EB_z0/MS_EB_z0_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_lambda_results.txt',sep=''),quote=F,row.names=F,sep='\t')
  ###############this is for ndr
  pdf(paste('./output/clade_analyses/MS_EB_z0/MS_EB_z0_ndr_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
  par(mfrow=c(2,2))
  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$z0,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$best.model.sigsq)
  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.ms05','ndr.ms09','ndr.ms05','ndr.ms09','sigsq')
  phylogroups.tree$node.label<-NULL
  #this is for caper
  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
  #for mean family seed weight
  cat('running pgls','\n')
  pgls.seedndr<-try(pgls(log10(ndr.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
  #summary(pgls.seedndr)
  pgls.seedsigsq<-try(pgls(log10(ndr.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
  #summary(pgls.seedsigsq)
  ####try both optimisations (caper + nlme if caper gives an error)
  if(class(pgls.seedndr)!='try-error'){
    plot(log10(phylogroups.tree.object$data$ndr.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
    abline(pgls.seedndr)
    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
  }
  if(class(pgls.seedndr)=='try-error'){
    row.names(phylogroups.data)<-phylogroups.data$Tip
    pgls.seedndr<-gls(log10(ndr.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
    plot(log10(phylogroups.data$ndr.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
    abline(pgls.seedndr)
    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
  }
  if(class(pgls.seedsigsq)!='try-error'){
    plot(log10(phylogroups.tree.object$data$ndr.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
    abline(pgls.seedsigsq)
    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
  }
  if(class(pgls.seedsigsq)=='try-error'){
    row.names(phylogroups.data)<-phylogroups.data$Tip
    pgls.seedsigsq<-gls(log10(ndr.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
    plot(log10(phylogroups.data$ndr.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
    abline(pgls.seedsigsq)
    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
  }
  
  dev.off()
  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
  colnames(results.df)<-c('analysis','slope','pvalue')
  rownames(results.df)<-NULL
  write.table(results.df,file=paste('./output/clade_analyses/MS_EB_z0/MS_EB_z0_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_ndr_results.txt',sep=''),quote=F,row.names=F,sep='\t')
  
}

####run_clades_RPANDA_EB_pgls_savedphylogroups_6models<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  load(file=paste('phylogroups_',minage,'_',maxage,'.Rsave',sep=''))
####  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####  length(unique(phylogroups.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_EB_phylogroup_select_size_6models(x,table,sampling,mincladesize))
####  df.phylogroups <- do.call("rbind", phylogroups.RPANDA)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.rpanda value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.rpanda1),]
####  #count the number of seed size data points per phylogroup
####  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
####  }
####  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
####  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
####  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
####  #get the mean seed size in each phylogroup
####  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
####  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  #AIC model selection,AICcBM-AICcEB
####  df.phylogroups$diff.AICc.BM.EB<-df.phylogroups$aicc.BM-df.phylogroups$aicc.EB
####  df.phylogroups$best.model.sigsq<-NA
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'best.model.sigsq']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'sigsq.BM']
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'best.model.sigsq']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'sigsq.EB']
####  
####  df.phylogroups$output.best.model<-NA
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'output.best.model']<-'BM'
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'output.best.model']<-'EB'
####  output.best.model<-as.data.frame(table(df.phylogroups$output.best.model))
####  write.table(output.best.model,file=paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_lambda_6models_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'best_fit_BMEB_model.txt',sep=''),quote=F,sep='\t',row.names=F)
####  
####  
####  ###############correlation with lambda
####  pdf(paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_lambda_6models_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$best.model.sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.rpanda','mu.rpanda','ndr.rpanda','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #output a table with the best fitting models for RPANDA for each clade
####  best.fitting.RPANDA<-as.data.frame(table(df.phylogroups$aicc),stringsAsFactors = F)
####  write.table(best.fitting.RPANDA,file=paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_6models_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_bestfitmodel_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedlambda<-try(pgls(log10(lambda.rpanda)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  pgls.seedsigsq<-try(pgls(log10(lambda.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedlambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.rpanda)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedlambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedlambda<-gls(log10(lambda.rpanda) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.rpanda)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(lambda.rpanda) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.rpanda)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_6models_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_lambda_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####  ###########this is for ndr
####  pdf(paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_ndr_analysis_clades_6models_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$ndr.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$best.model.sigsq)
####  
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.rpanda','mu.rpanda','ndr.rpanda','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #####this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedndr<-try(pgls(log10(ndr.rpanda)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedndr)
####  pgls.seedsigsq<-try(pgls(log10(ndr.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  if(class(pgls.seedndr)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.rpanda)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedndr)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedndr<-gls(log10(ndr.rpanda) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.rpanda)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(ndr.rpanda) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.rpanda)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_6models_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_ndr_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####}
####run_clades_RPANDA_EB_pgls_savedphylogroups<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  load(file=paste('phylogroups_',minage,'_',maxage,'.Rsave',sep=''))
####  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####  length(unique(phylogroups.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_EB_phylogroup_select_size(x,table,sampling,mincladesize))
####  df.phylogroups <- do.call("rbind", phylogroups.RPANDA)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.rpanda value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.rpanda1),]
####  #count the number of seed size data points per phylogroup
####  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
####  }
####  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
####  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
####  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
####  #get the mean seed size in each phylogroup
####  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
####  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  #AIC model selection,AICcBM-AICcEB
####  df.phylogroups$diff.AICc.BM.EB<-df.phylogroups$aicc.BM-df.phylogroups$aicc.EB
####  df.phylogroups$best.model.sigsq<-NA
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'best.model.sigsq']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'sigsq.BM']
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'best.model.sigsq']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'sigsq.EB']
####  
####  df.phylogroups$output.best.model<-NA
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'output.best.model']<-'BM'
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'output.best.model']<-'EB'
####  output.best.model<-as.data.frame(table(df.phylogroups$output.best.model))
####  write.table(output.best.model,file=paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_lambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'best_fit_BMEB_model.txt',sep=''),quote=F,sep='\t',row.names=F)
####  
####  
####  ###############correlation with lambda
####  pdf(paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_lambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$best.model.sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.rpanda','mu.rpanda','ndr.rpanda','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #output a table with the best fitting models for RPANDA for each clade
####  best.fitting.RPANDA<-as.data.frame(table(df.phylogroups$aicc),stringsAsFactors = F)
####  write.table(best.fitting.RPANDA,file=paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_bestfitmodel_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedlambda<-try(pgls(log10(lambda.rpanda)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  pgls.seedsigsq<-try(pgls(log10(lambda.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedlambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.rpanda)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedlambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedlambda<-gls(log10(lambda.rpanda) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.rpanda)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(lambda.rpanda) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.rpanda)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_lambda_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####  ###########this is for ndr
####  pdf(paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_ndr_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$ndr.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$best.model.sigsq)
####  
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.rpanda','mu.rpanda','ndr.rpanda','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #####this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedndr<-try(pgls(log10(ndr.rpanda)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedndr)
####  pgls.seedsigsq<-try(pgls(log10(ndr.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  if(class(pgls.seedndr)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.rpanda)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedndr)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedndr<-gls(log10(ndr.rpanda) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.rpanda)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(ndr.rpanda) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.rpanda)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_ndr_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####}


run_clades_RPANDA_EB_z0_pgls_savedphylogroups_6models<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
  colnames(table)[4]<-'tip'
  colnames(table)[9]<-'genus.species.sampling.fraction'
  cat('getting clades','\n')
  load(file=paste('phylogroups_',minage,'_',maxage,'.Rsave',sep=''))
  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
  length(unique(phylogroups.species))
  cat('measuring diversification/trait evolution','\n')
  phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_EB_z0_phylogroup_select_size_6models(x,table,sampling,mincladesize))
  df.phylogroups <- do.call("rbind", phylogroups.RPANDA)
  df.phylogroups<-unique(df.phylogroups)
  #remove rows with no lambda.rpanda value
  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.rpanda1),]
  #count the number of seed size data points per phylogroup
  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
  }
  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
  #get the mean seed size in each phylogroup
  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
  #getting vector of representatives in each phylogroup
  #then keep just one tip per phylogroup
  phylogroup.clades.species<-vector('character')
  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
  }
  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
  colnames(names.phylogroups.tree)<-'Tip'
  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
  #another filter to check that size of clade is correct
  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
  #add another filter to select well sampled clades
  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
  #AIC model selection,AICcBM-AICcEB
  df.phylogroups$diff.AICc.BM.EB<-df.phylogroups$aicc.BM-df.phylogroups$aicc.EB
  df.phylogroups$best.model.sigsq<-NA
  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'best.model.sigsq']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'sigsq.BM']
  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'best.model.sigsq']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'sigsq.EB']
  
  df.phylogroups$output.best.model<-NA
  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'output.best.model']<-'BM'
  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'output.best.model']<-'EB'
  
  df.phylogroups$z0<-NA
  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'z0']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'z0.BM']
  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'z0']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'z0.EB']
  correlation.z0.SeedWeight<-cor.test(df.phylogroups$z0,df.phylogroups$phylogroup.Log10.Seed.Weight) 
  correlation.df<-cbind(correlation.z0.SeedWeight$estimate,correlation.z0.SeedWeight$p.value)
  colnames(correlation.df)<-c('estimate','p.value')
  write.table(correlation.df,file=paste('./output/clade_analyses/MS_EB_z0/MS_EB_z0_lambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'cor_z0_SeedWeight.txt',sep=''),quote=F,sep='\t',row.names=F)
  
  output.best.model<-as.data.frame(table(df.phylogroups$output.best.model))
  write.table(output.best.model,file=paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_lambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'best_fit_BMEB_model.txt',sep=''),quote=F,sep='\t',row.names=F)
  
  
  ###############correlation with lambda
  pdf(paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_lambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
  par(mfrow=c(2,2))
  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$z0,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$best.model.sigsq)
  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.rpanda','mu.rpanda','ndr.rpanda','sigsq')
  phylogroups.tree$node.label<-NULL
  #output a table with the best fitting models for RPANDA for each clade
  best.fitting.RPANDA<-as.data.frame(table(df.phylogroups$aicc),stringsAsFactors = F)
  write.table(best.fitting.RPANDA,file=paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_bestfitmodel_results.txt',sep=''),quote=F,row.names=F,sep='\t')
  #this is for caper
  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
  #for mean family seed weight
  cat('running pgls','\n')
  pgls.seedlambda<-try(pgls(log10(lambda.rpanda)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
  #summary(pgls.seedlambda)
  pgls.seedsigsq<-try(pgls(log10(lambda.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
  #summary(pgls.seedsigsq)
  ####try both optimisations (caper + nlme if caper gives an error)
  if(class(pgls.seedlambda)!='try-error'){
    plot(log10(phylogroups.tree.object$data$lambda.rpanda)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
    abline(pgls.seedlambda)
    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
  }
  if(class(pgls.seedlambda)=='try-error'){
    row.names(phylogroups.data)<-phylogroups.data$Tip
    pgls.seedlambda<-gls(log10(lambda.rpanda) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
    plot(log10(phylogroups.data$lambda.rpanda)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
    abline(pgls.seedlambda)
    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
  }
  if(class(pgls.seedsigsq)!='try-error'){
    plot(log10(phylogroups.tree.object$data$lambda.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
    abline(pgls.seedsigsq)
    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
  }
  if(class(pgls.seedsigsq)=='try-error'){
    row.names(phylogroups.data)<-phylogroups.data$Tip
    pgls.seedsigsq<-gls(log10(lambda.rpanda) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
    plot(log10(phylogroups.data$lambda.rpanda)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
    abline(pgls.seedsigsq)
    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
  }
  
  dev.off()
  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
  colnames(results.df)<-c('analysis','slope','pvalue')
  rownames(results.df)<-NULL
  write.table(results.df,file=paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_lambda_results.txt',sep=''),quote=F,row.names=F,sep='\t')
  
  ###########this is for ndr
  pdf(paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_ndr_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
  par(mfrow=c(2,2))
  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$z0,df.phylogroups$ndr.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$best.model.sigsq)
  
  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.rpanda','mu.rpanda','ndr.rpanda','sigsq')
  phylogroups.tree$node.label<-NULL
  #####this is for caper
  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
  #for mean family seed weight
  cat('running pgls','\n')
  pgls.seedndr<-try(pgls(log10(ndr.rpanda)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
  #summary(pgls.seedndr)
  pgls.seedsigsq<-try(pgls(log10(ndr.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
  #summary(pgls.seedsigsq)
  if(class(pgls.seedndr)!='try-error'){
    plot(log10(phylogroups.tree.object$data$ndr.rpanda)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
    abline(pgls.seedndr)
    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
  }
  if(class(pgls.seedndr)=='try-error'){
    row.names(phylogroups.data)<-phylogroups.data$Tip
    pgls.seedndr<-gls(log10(ndr.rpanda) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
    plot(log10(phylogroups.data$ndr.rpanda)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
    abline(pgls.seedndr)
    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
  }
  if(class(pgls.seedsigsq)!='try-error'){
    plot(log10(phylogroups.tree.object$data$ndr.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
    abline(pgls.seedsigsq)
    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
  }
  if(class(pgls.seedsigsq)=='try-error'){
    row.names(phylogroups.data)<-phylogroups.data$Tip
    pgls.seedsigsq<-gls(log10(ndr.rpanda) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
    plot(log10(phylogroups.data$ndr.rpanda)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
    abline(pgls.seedsigsq)
    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
  }
  
  dev.off()
  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
  colnames(results.df)<-c('analysis','slope','pvalue')
  rownames(results.df)<-NULL
  write.table(results.df,file=paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_ndr_results.txt',sep=''),quote=F,row.names=F,sep='\t')
}

###wrapper to run clade analysis with RPANDA + pgls
#this creates the phylogroups (there's another script to run if phyogroups have previously been saved)
####run_clades_RPANDA_pgls<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=mincladesize,ncores=ncores)
####  #save phylogroup object to Rsave file
####  save(phylogroups,file=paste('phylogroups_',minage,'_',maxage,'.Rsave',sep=''))
####  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####  length(unique(phylogroups.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_BM_phylogroup_select_size(x,table,sampling,mincladesize))
####  #phylogroups.RPANDA.nomuexp<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_BM_phylogroup_select_size(x,table,sampling,mincladesize))
####  df.phylogroups <- do.call("rbind", phylogroups.RPANDA)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.rpanda value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.rpanda1),]
####  #count the number of seed size data points per phylogroup
####  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
####  }
####  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
####  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
####  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
####  #get the mean seed size in each phylogroup
####  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
####  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  ###############correlation with lambda
####  pdf(paste('./output/plots/clade_analyses_new/RPANDAlambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.rpanda','mu.rpanda','ndr.rpanda','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedlambda<-try(pgls(log10(lambda.rpanda)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  pgls.seedsigsq<-try(pgls(log10(lambda.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedlambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.rpanda)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedlambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedlambda<-gls(log10(lambda.rpanda) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.rpanda)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(lambda.rpanda) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.rpanda)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_RPANDAlambda_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####  ###########this is for ndr
####  pdf(paste('./output/plots/clade_analyses_new/RPANDAndr_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$ndr.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$sigsq)
####  
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.rpanda','mu.rpanda','ndr.rpanda','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #####this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedndr<-try(pgls(log10(ndr.rpanda)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedndr)
####  pgls.seedsigsq<-try(pgls(log10(ndr.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  if(class(pgls.seedndr)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.rpanda)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedndr)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedndr<-gls(log10(ndr.rpanda) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.rpanda)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(ndr.rpanda) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.rpanda)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_RPANDAndr_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####}
#this will load presaved phylogroups
####run_clades_RPANDA_pgls_savedphylogroups<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  load(file=paste('phylogroups_',minage,'_',maxage,'.Rsave',sep=''))
####  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####  length(unique(phylogroups.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_BM_phylogroup_select_size(x,table,sampling,mincladesize))
####  #phylogroups.RPANDA.nomuexp<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_BM_phylogroup_select_size(x,table,sampling,mincladesize))
####  df.phylogroups <- do.call("rbind", phylogroups.RPANDA)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.rpanda value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.rpanda1),]
####  #count the number of seed size data points per phylogroup
####  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
####  }
####  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
####  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
####  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
####  #get the mean seed size in each phylogroup
####  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
####  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  ###############correlation with lambda
####  pdf(paste('./output/plots/clade_analyses_new/RPANDAlambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.rpanda','mu.rpanda','ndr.rpanda','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedlambda<-try(pgls(log10(lambda.rpanda)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  pgls.seedsigsq<-try(pgls(log10(lambda.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedlambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.rpanda)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedlambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedlambda<-gls(log10(lambda.rpanda) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.rpanda)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(lambda.rpanda) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.rpanda)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_RPANDAlambda_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####  ###########this is for ndr
####  pdf(paste('./output/plots/clade_analyses_new/RPANDAndr_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$ndr.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$sigsq)
####  
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.rpanda','mu.rpanda','ndr.rpanda','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #####this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedndr<-try(pgls(log10(ndr.rpanda)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedndr)
####  pgls.seedsigsq<-try(pgls(log10(ndr.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  if(class(pgls.seedndr)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.rpanda)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedndr)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedndr<-gls(log10(ndr.rpanda) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.rpanda)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(ndr.rpanda) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.rpanda)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_RPANDAndr_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####}

###wrapper to run clade analysis (with Magallon & Sanderson estimator) for congeneric species only
#this creates the phylogroups (there's another script to run if phyogroups have previously been saved)
####run_clades_MSlambda_congenerics_pgls<-function(minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=mincladesize,ncores=ncores)
####  #save phylogroup object to Rsave file
####  #check for phylogroups that only contain congeneric species
####  phylogroups.congenerics.positions<-unlist(lapply(phylogroups$phylogroups,function(x){genera<-length(unique(sapply(strsplit(as.character(x$tip.label),'_'),function(x) x[1]))==1)}))
####  phylogroups.congenerics<-phylogroups$phylogroups[which(phylogroups.congenerics.positions==1)]
####  phylogroups.congenerics.species<-unlist(lapply(phylogroups.congenerics,function(x)x$tip.label))
####  length(unique(phylogroups.congenerics.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.MS<-lapply(phylogroups$phylogroups,function(x) run_MS_BM_phylogroup_select_size(x,table))
####  df.phylogroups <- do.call("rbind", phylogroups.MS)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.ms value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.ms0),]
####  
####  #count the number of seed size data points per phylogroup
####  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
####  }
####  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
####  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
####  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
####  
####  #get the mean seed size in each phylogroup
####  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
####  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  ###############this is for lambda
####  pdf(paste('./output/plots/clade_analyses_new/MSlambda_congenerics_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$lambda.ms0,df.phylogroups$lambda.ms05,df.phylogroups$lambda.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.ms0','lambda.ms05','lambda.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedlambda<-try(pgls(log10(lambda.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  pgls.seedsigsq<-try(pgls(log10(lambda.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedlambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedlambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedlambda<-gls(log10(lambda.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(lambda.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSlambda_congenerics_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  ###############this is for ndr
####  pdf(paste('./output/plots/clade_analyses_new/MSndr_congenerics_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$ndr.ms0,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.ms0','ndr.ms05','ndr.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedndr<-try(pgls(log10(ndr.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedndr)
####  pgls.seedsigsq<-try(pgls(log10(ndr.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedndr)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedndr)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedndr<-gls(log10(ndr.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(ndr.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSndr_congenerics_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####}
#####this will load presaved phylogroups
####run_clades_MSlambda_congenerics_pgls_savedphylogroups<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  #phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=mincladesize,ncores=ncores)
####  #save phylogroup object to Rsave file
####  load(file=paste('phylogroups_',minage,'_',maxage,'.Rsave',sep=''))
####  #check for phylogroups that only contain congeneric species
####  phylogroups.congenerics.positions<-unlist(lapply(phylogroups$phylogroups,function(x){genera<-length(unique(sapply(strsplit(as.character(x$tip.label),'_'),function(x) x[1]))==1)}))
####  phylogroups.congenerics<-phylogroups$phylogroups[which(phylogroups.congenerics.positions==1)]
####  phylogroups.congenerics.species<-unlist(lapply(phylogroups.congenerics,function(x)x$tip.label))
####  length(unique(phylogroups.congenerics.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.MS<-lapply(phylogroups$phylogroups,function(x) run_MS_BM_phylogroup_select_size(x,table,mincladesize = mincladesize,sampling=sampling))
####  df.phylogroups <- do.call("rbind", phylogroups.MS)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.ms value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.ms0),]
####  
####  #count the number of seed size data points per phylogroup
####  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
####  }
####  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
####  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
####  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
####  
####  #get the mean seed size in each phylogroup
####  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
####  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  ###############this is for lambda
####  pdf(paste('./output/plots/clade_analyses_new/MSlambda_congenerics_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$lambda.ms0,df.phylogroups$lambda.ms05,df.phylogroups$lambda.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.ms0','lambda.ms05','lambda.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedlambda<-try(pgls(log10(lambda.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  pgls.seedsigsq<-try(pgls(log10(lambda.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedlambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedlambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedlambda<-gls(log10(lambda.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(lambda.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSlambda_congenerics_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  ###############this is for ndr
####  pdf(paste('./output/plots/clade_analyses_new/MSndr_congenerics_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.ms05','ndr.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedndr<-try(pgls(log10(ndr.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedndr)
####  pgls.seedsigsq<-try(pgls(log10(ndr.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedndr)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedndr)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedndr<-gls(log10(ndr.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(ndr.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSndr_congenerics_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
}
run_clades_MSlambda_z0_congenerics_pgls_savedphylogroups<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
  colnames(table)[4]<-'tip'
  colnames(table)[9]<-'genus.species.sampling.fraction'
  cat('getting clades','\n')
  #phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=mincladesize,ncores=ncores)
  #save phylogroup object to Rsave file
  load(file=paste('phylogroups_',minage,'_',maxage,'.Rsave',sep=''))
  #check for phylogroups that only contain congeneric species
  phylogroups.congenerics.positions<-unlist(lapply(phylogroups$phylogroups,function(x){genera<-length(unique(sapply(strsplit(as.character(x$tip.label),'_'),function(x) x[1]))==1)}))
  phylogroups.congenerics<-phylogroups$phylogroups[which(phylogroups.congenerics.positions==1)]
  phylogroups.congenerics.species<-unlist(lapply(phylogroups.congenerics,function(x)x$tip.label))
  length(unique(phylogroups.congenerics.species))
  cat('measuring diversification/trait evolution','\n')
  phylogroups.MS<-lapply(phylogroups.congenerics,function(x) run_MS_EB_z0_phylogroup_select_size(x,table,mincladesize = mincladesize,sampling=sampling))
  df.phylogroups <- do.call("rbind", phylogroups.MS)
  df.phylogroups<-unique(df.phylogroups)
  #remove rows with no lambda.ms value
  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.ms0),]
  
  #count the number of seed size data points per phylogroup
  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
  }
  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
  
  #get the mean seed size in each phylogroup
  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
  #getting vector of representatives in each phylogroup
  #then keep just one tip per phylogroup
  phylogroup.clades.species<-vector('character')
  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
  }
  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
  colnames(names.phylogroups.tree)<-'Tip'
  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
  #another filter to check that size of clade is correct
  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
  #add another filter to select well sampled clades
  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
  #AIC model selection,AICcBM-AICcEB
  
  df.phylogroups$diff.AICc.BM.EB<-df.phylogroups$aicc.BM-df.phylogroups$aicc.EB
  df.phylogroups$best.model.sigsq<-NA
  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'best.model.sigsq']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'sigsq.BM']
  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'best.model.sigsq']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'sigsq.EB']
  df.phylogroups$output.best.model<-NA
  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'output.best.model']<-'BM'
  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'output.best.model']<-'EB'
  output.best.model<-as.data.frame(table(df.phylogroups$output.best.model))
  df.phylogroups$z0<-NA
  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'z0']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'z0.BM']
  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'z0']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'z0.EB']
  correlation.z0.SeedWeight<-cor.test(df.phylogroups$z0,df.phylogroups$phylogroup.Log10.Seed.Weight) 
  correlation.df<-cbind(correlation.z0.SeedWeight$estimate,correlation.z0.SeedWeight$p.value)
  colnames(correlation.df)<-c('estimate','p.value')
  write.table(correlation.df,file=paste('./output/clade_analyses/MS_EB_z0/MS_EB_z0_lambda_congenerics_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'cor_z0_SeedWeight.txt',sep=''),quote=F,sep='\t',row.names=F)
  write.table(output.best.model,file=paste('./output/clade_analyses/MS_EB_z0/MS_EB_z0_lambda_congenerics_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'best_fit_BMEB_model.txt',sep=''),quote=F,sep='\t',row.names=F)
  
  
  ###############this is for lambda
  pdf(paste('./output/plots/MSlambda_congenerics_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
  par(mfrow=c(2,2))
  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$z0,df.phylogroups$lambda.ms0,df.phylogroups$lambda.ms05,df.phylogroups$lambda.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$best.model.sigsq)
  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.ms0','lambda.ms05','lambda.ms09','ndr.ms05','ndr.ms09','sigsq')
  phylogroups.tree$node.label<-NULL
  #this is for caper
  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
  #for mean family seed weight
  cat('running pgls','\n')
  pgls.seedlambda<-try(pgls(log10(lambda.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
  #summary(pgls.seedlambda)
  pgls.seedsigsq<-try(pgls(log10(lambda.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
  #summary(pgls.seedsigsq)
  ####try both optimisations (caper + nlme if caper gives an error)
  if(class(pgls.seedlambda)!='try-error'){
    plot(log10(phylogroups.tree.object$data$lambda.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
    abline(pgls.seedlambda)
    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
  }
  if(class(pgls.seedlambda)=='try-error'){
    row.names(phylogroups.data)<-phylogroups.data$Tip
    pgls.seedlambda<-gls(log10(lambda.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
    plot(log10(phylogroups.data$lambda.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
    abline(pgls.seedlambda)
    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
  }
  if(class(pgls.seedsigsq)!='try-error'){
    plot(log10(phylogroups.tree.object$data$lambda.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
    abline(pgls.seedsigsq)
    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
  }
  if(class(pgls.seedsigsq)=='try-error'){
    row.names(phylogroups.data)<-phylogroups.data$Tip
    pgls.seedsigsq<-gls(log10(lambda.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
    plot(log10(phylogroups.data$lambda.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
    abline(pgls.seedsigsq)
    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
  }
  
  dev.off()
  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
  colnames(results.df)<-c('analysis','slope','pvalue')
  rownames(results.df)<-NULL
  write.table(results.df,file=paste('./output/plots/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSlambda_congenerics_results.txt',sep=''),quote=F,row.names=F,sep='\t')
  ###############this is for ndr
  pdf(paste('./output/clade_analyses/MS_EB_z0/MSndr_congenerics_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
  par(mfrow=c(2,2))
  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$z0,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$best.model.sigsq)
  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.ms05','ndr.ms09','ndr.ms05','ndr.ms09','sigsq')
  phylogroups.tree$node.label<-NULL
  #this is for caper
  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
  #for mean family seed weight
  cat('running pgls','\n')
  pgls.seedndr<-try(pgls(log10(ndr.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
  #summary(pgls.seedndr)
  pgls.seedsigsq<-try(pgls(log10(ndr.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
  #summary(pgls.seedsigsq)
  ####try both optimisations (caper + nlme if caper gives an error)
  if(class(pgls.seedndr)!='try-error'){
    plot(log10(phylogroups.tree.object$data$ndr.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
    abline(pgls.seedndr)
    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
  }
  if(class(pgls.seedndr)=='try-error'){
    row.names(phylogroups.data)<-phylogroups.data$Tip
    pgls.seedndr<-gls(log10(ndr.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
    plot(log10(phylogroups.data$ndr.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
    abline(pgls.seedndr)
    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
  }
  if(class(pgls.seedsigsq)!='try-error'){
    plot(log10(phylogroups.tree.object$data$ndr.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
    abline(pgls.seedsigsq)
    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
  }
  if(class(pgls.seedsigsq)=='try-error'){
    row.names(phylogroups.data)<-phylogroups.data$Tip
    pgls.seedsigsq<-gls(log10(ndr.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
    plot(log10(phylogroups.data$ndr.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
    abline(pgls.seedsigsq)
    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
  }
  
  dev.off()
  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
  colnames(results.df)<-c('analysis','slope','pvalue')
  rownames(results.df)<-NULL
  write.table(results.df,file=paste('./output/clade_analyses/MS_EB_z0/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSndr_congenerics_results.txt',sep=''),quote=F,row.names=F,sep='\t')
  
}

cat('preparing dataset to run clade based analyses','\n')
#Qian and Jin created an updated version of Zanne tree - 'Appendix S3-v2' or 'QianTree.txt'
#reference doi: 10.1093/jpe/rtv047
Qian<-read.tree('./raw_data/QianTree.txt')
Qian_TPL<-read.table('./output/tables/Qian_dataset_TPL.txt',header=T,sep='\t',stringsAsFactors = F)
Qian_TPL$Merged<-paste(Qian_TPL$New.Genus,Qian_TPL$New.Species,sep='_')
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
Qian_Lookup_Angios_TPL<-Qian_Lookup_Angios_TPL[,c(c(1:6),23)]
#drop tips in Qian tree (not angios, duplicated, hybrids)
remove_tips<-setdiff(Qian$tip.label,Qian_Lookup_Angios_TPL$Old_Species)
Qian_dropped<-drop.tip(Qian,remove_tips)
#remove "_sp" tips
remove_sp_tips<-Qian_Lookup_Angios_TPL[grep('_sp$',Qian_Lookup_Angios_TPL$Old_Species),]$Old_Species
Qian_dropped<-drop.tip(Qian_dropped,remove_sp_tips)
Qian_Lookup_Angios_TPL<-Qian_Lookup_Angios_TPL[-grep('_sp$',Qian_Lookup_Angios_TPL$New_Species),]
Qiandf<-data.frame(Qian_dropped$tip.label)
colnames(Qiandf)<-'Tipname'
#replace tip names in Qian Tree with TPL+taxonlookup alternative
Qian_Lookup_Angios_TPL<-merge(Qiandf,Qian_Lookup_Angios_TPL,by.x='Tipname',by.y='Old_Species')
Qian_Lookup_Angios_TPL$Tipname<-as.character(Qian_Lookup_Angios_TPL$Tipname)

species.counts<-as.data.frame(table(Qian_Lookup_Angios_TPL$genus),stringsAsFactors = F)
colnames(species.counts)<-c('genus','number.of.species.tree')
Qian_Lookup_Angios_TPL<-merge(Qian_Lookup_Angios_TPL,species.counts)
Qian_Lookup_Angios_TPL$genus.sampling.fraction<-Qian_Lookup_Angios_TPL$number.of.species.tree/Qian_Lookup_Angios_TPL$number.of.species
Qian_Lookup_Angios_TPL<-Qian_Lookup_Angios_TPL[match(Qian_dropped$tip.label,Qian_Lookup_Angios_TPL$Tipname),]
Qian_dropped$tip.label<-Qian_Lookup_Angios_TPL$New_Species

tree<-Qian_dropped
#get total info from TaxonLookUp
plant_lookup_version_current()
PlantLookup<-plant_lookup(include_counts = TRUE)
#get number of genera in each family
genus.counts<-as.data.frame(table(PlantLookup$family),stringsAsFactors = F)
colnames(genus.counts)<-c('family','number.of.genera')

#run TaxonLookup on taxa from the tree
tree.table<-lookup_table(tree$tip.label)
genus.counts.tree<-as.data.frame(table(tree.table$family),stringsAsFactors = F)
colnames(genus.counts.tree)<-c('family','number.of.genera.tree')

#merge all info
genus.counts<-merge(genus.counts,genus.counts.tree)
genus.counts$family.sampling.fraction<-genus.counts$number.of.genera.tree/genus.counts$number.of.genera
#merge with all tree info
Qian_Lookup_Angios_TPL<-merge(Qian_Lookup_Angios_TPL,genus.counts,all.x=TRUE)
#add info on seed size (not discarding species WITHOUT seed size data)
kewsid.TPL<-read.table('./output/tables/kewsid_TPL_seed.txt',header=T,sep='\t')
kewsid.TPL<-kewsid.TPL[,c('New.Merged',"Log10.Seed.Weight")]
colnames(kewsid.TPL)[1]<-c('New_Species')
kewsid.TPL<-unique(kewsid.TPL)
#this solves an issue with a duplicate Chenopodium_murale and Lessingia_germanorum
kewsid.TPL<-subset(kewsid.TPL,!duplicated(kewsid.TPL$New_Species))
#this adds seed size data to Qian dataset BUT doesn't discard species with no seed size data
Qian_Lookup_Angios_TPL_seed<-merge(Qian_Lookup_Angios_TPL,kewsid.TPL,all.x=TRUE)
#adding info on taxonomic discordance (when there are more tips in the tree than accepted species according to taxonlookup)
#it's only 555 species(1% of the total) // 79 genera (1% of the total)
Qian_Lookup_Angios_TPL_seed$taxonomic.discordance<-NA
Qian_Lookup_Angios_TPL_seed[Qian_Lookup_Angios_TPL_seed$number.of.species<Qian_Lookup_Angios_TPL_seed$number.of.species.tree,]$taxonomic.discordance<-1
Qian_Lookup_Angios_TPL_seed[Qian_Lookup_Angios_TPL_seed$number.of.species>=Qian_Lookup_Angios_TPL_seed$number.of.species.tree,]$taxonomic.discordance<-0
#to simplify I make the genus.species.sampling.fraction for these species = 1. But they can be removed as well for further analyses
#Qian_Lookup_Angios_TPL_seed[Qian_Lookup_Angios_TPL_seed$taxonomic.discordance==1,]$genus.sampling.fraction<-1
#remove species with taxonomic discordance
Qian_Lookup_Angios_TPL_seed<-Qian_Lookup_Angios_TPL_seed[!Qian_Lookup_Angios_TPL_seed$taxonomic.discordance==1,]
Qian_dropped<-drop.tip(Qian_dropped,setdiff(Qian_dropped$tip.label,Qian_Lookup_Angios_TPL_seed$New_Species))
write.tree(Qian_dropped,file='./Qian_dropped.tree')






