library(caper)

#traitdatacov: "BAMMPhyndr_GenusSeed_cov_trait_data.txt"

mean_cov_seed_mass<-function(treefile,traitdatacov){
  tree<-read.tree(treefile)
  trait_input.cov<-read.table(traitdatacov,sep='\t',header=T)
  trait_input.cov<-na.omit(trait_input.cov)
  remove.tips<-setdiff(tree$tip.label,trait_input.cov$tip.label)
  tree.cov<-drop.tip(tree,remove.tips)
  tree.cov$node.label<-NULL
  trait_input.cov$mean<-log10(trait_input.cov$mean)
  #plot(trait_input.cov$mean,log10(trait_input.cov$cov),las=1)
  cov.data<-comparative.data(tree.cov,trait_input.cov,names.col = 'tip.label')
  #plot(trait_input.cov$mean,log10(trait_input.cov$cov),las=1)
  pgls.cov.seedweight<-pgls(log10(mean)~cov, data=cov.data)
  write(as.character(summary(pgls.cov.seedweight)),file='./output/tables/pgls_mean_cov.txt',sep='\t')
}

