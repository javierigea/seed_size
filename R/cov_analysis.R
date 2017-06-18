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

