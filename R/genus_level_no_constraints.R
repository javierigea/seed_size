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

#####this script runs the genera analysis with the species level tree
#####it selects monophyletic genera with seed size data in the Zanne tree
#####and measures diversification with Magallon & Sanderson
#####and runs pgls
cat ('running genus level analysis - with no constraints for species with seed size data','\n')
cat ('getting phylogenetic data','\n')

#Qian and Jin created an updated version of Zanne tree - 'Appendix S3-v2' or 'QianTree.txt'
#reference doi: 10.1093/jpe/rtv047
Qian<-read.tree('./raw_data/QianTree.txt')
Qian_TPL<-read.table('./output/tables/Qian_dataset_TPL.txt',header=T,sep='\t',stringsAsFactors = F)
Qian_TPL$Merged<-paste(Qian_TPL$Genus,Qian_TPL$Species,sep='_')

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
Qiandf<-data.frame(Qian_dropped$tip.label)
colnames(Qiandf)<-'Tipname'
#replace tip names in Qian Tree with TPL+taxonlookup alternative
Qian_Lookup_Angios_TPL<-merge(Qiandf,Qian_Lookup_Angios_TPL,by.x='Tipname',by.y='Old_Species')
Qian_Lookup_Angios_TPL$Tipname<-as.character(Qian_Lookup_Angios_TPL$Tipname)
Qian_Lookup_Angios_TPL<-Qian_Lookup_Angios_TPL[match(Qian_dropped$tip.label,Qian_Lookup_Angios_TPL$Tipname),]
Qian_dropped$tip.label<-Qian_Lookup_Angios_TPL$New_Species
#count the number of species per genera in the tree
table.Genus<-table(Qian_Lookup_Angios_TPL$genus)
table.Genus<-data.frame(names(table.Genus),unname(table.Genus))
table.Genus<-table.Genus[,c(1,3)]
colnames(table.Genus)<-c('genus','n.species.tree')
Qian_Lookup_Angios_TPL<-merge(Qian_Lookup_Angios_TPL,table.Genus)
#check monophyly of genera
cat ('checking monophyly of genera','\n')
monophyly.table<-matrix(NA,nrow=0,ncol=2)
for (i in 1:length(unique(Qian_Lookup_Angios_TPL$genus))){
  #cat(unique(Qian_Lookup_Angios_TPL$genus)[i],'\n')
  monophyly.genus<-is.monophyletic(Qian_dropped,Qian_Lookup_Angios_TPL[Qian_Lookup_Angios_TPL$genus %in% unique(Qian_Lookup_Angios_TPL$genus)[i],]$New_Species)
  monophyly.table<-rbind(monophyly.table,c(unique(Qian_Lookup_Angios_TPL$genus)[i],monophyly.genus))
}
#add the monophyly info to the table
monophyly.table<-as.data.frame(monophyly.table)
colnames(monophyly.table)<-c('genus','monophyly')
Qian_Lookup_Angios_TPL<-merge(Qian_Lookup_Angios_TPL,monophyly.table)
#add info on the number of described species in each genera
plant_lookup_version_current()
PlantLookup<-plant_lookup(include_counts = TRUE)
PlantLookup<-PlantLookup[,c(2,1)]
colnames(PlantLookup)[2]<-c('n.described.species')
Qian_Lookup_Angios_TPL<-merge(Qian_Lookup_Angios_TPL,PlantLookup)
Qian_Lookup_Angios_TPL$genus.sampling.fraction<-Qian_Lookup_Angios_TPL$n.species.tree/Qian_Lookup_Angios_TPL$n.described.species
#there are some genera with more species in the tree than described species
#create an extra column for this: mismatch_plantlookup
Qian_Lookup_Angios_TPL$mismatch_plantlookup<-0
Qian_Lookup_Angios_TPL[Qian_Lookup_Angios_TPL$n.species.tree>Qian_Lookup_Angios_TPL$n.described.species,]$mismatch_plantlookup<-1
#select genus with > 3 tips in tree
Qian_Lookup_Angios_TPL.selected<-Qian_Lookup_Angios_TPL[Qian_Lookup_Angios_TPL$n.species.tree>3,]
#select monophyletic genera
Qian_Lookup_Angios_TPL.selected<-Qian_Lookup_Angios_TPL.selected[(Qian_Lookup_Angios_TPL.selected$monophyly==TRUE),]
#select genera with > 0.3 sampling
Qian_Lookup_Angios_TPL.selected<-Qian_Lookup_Angios_TPL.selected[Qian_Lookup_Angios_TPL.selected$genus.sampling.fraction>0.3,]
#calculate mrca age
cat ('calculating genus ages','\n')
Qian_dropped$node.label<-NULL
bt <- branching.times(Qian_dropped)
genus.ages.matrix<-matrix(NA,nrow=0,ncol=2)
for (i in 1:length(unique(Qian_Lookup_Angios_TPL.selected$genus))){
  #cat(unique(Qian_Lookup_Angios_TPL.selected$genus)[i],'\n')
  mrca.node<-getMRCA(Qian_dropped,Qian_Lookup_Angios_TPL.selected[Qian_Lookup_Angios_TPL.selected$genus %in% unique(Qian_Lookup_Angios_TPL.selected$genus)[i],]$New_Species)
  genus.age<-bt[as.character(mrca.node)]
  genus.ages.matrix<-rbind(genus.ages.matrix,c(unique(Qian_Lookup_Angios_TPL.selected$genus)[i],genus.age))
}
#add the genus age info to the table
genus.ages.table<-as.data.frame(genus.ages.matrix)
colnames(genus.ages.table)<-c('genus','mrca.age')
Qian_Lookup_Angios_TPL.selected<-merge(Qian_Lookup_Angios_TPL.selected,genus.ages.table)
#filter genera with taxonlookup mismatch
Qian_Lookup_Angios_TPL.selected<-Qian_Lookup_Angios_TPL.selected[Qian_Lookup_Angios_TPL.selected$mismatch_plantlookup<1,]
#run MS on selected genera
cat('measuring diversification','\n')
ms.matrix<-matrix(NA,nrow=0,ncol=7)
for (i in 1:length(unique(Qian_Lookup_Angios_TPL.selected$genus))){
  #cat(unique(Qian_Lookup_Angios_TPL.selected$genus)[i],'\n')
  age<-as.numeric(Qian_Lookup_Angios_TPL.selected[match(unique(Qian_Lookup_Angios_TPL.selected$genus)[i],Qian_Lookup_Angios_TPL.selected$genus),]$mrca.age)
  species<-as.numeric(Qian_Lookup_Angios_TPL.selected[match(unique(Qian_Lookup_Angios_TPL.selected$genus)[i],Qian_Lookup_Angios_TPL.selected$genus),]$n.described.species)
  ms0<-bd.ms(time=age,n=species,epsilon=0)
  lambda.ms0<-ms0
  ms05<-bd.ms(time=age,n=species,epsilon=0.5)
  lambda.ms05<-ms05/(1-0.5)
  ms09<-bd.ms(time=age,n=species,epsilon=0.9)
  lambda.ms09<-ms09/(1-0.9)
  ms.matrix<-rbind(ms.matrix,c(unique(Qian_Lookup_Angios_TPL.selected$genus)[i],ms0,lambda.ms0,ms05,lambda.ms05,ms09,lambda.ms09))
}
ms.table<-as.data.frame(ms.matrix,stringsAsFactors = F)
colnames(ms.table)<-c('genus','ms0','lambda.ms0','ms05','lambda.ms05','ms09','lambda.ms09')
Qian_Lookup_Angios_TPL.selected<-merge(Qian_Lookup_Angios_TPL.selected,ms.table)

##measure RPANDA in genera
##rpanda.matrix<-matrix(NA,nrow=0,ncol=5)
#for (i in 1:length(unique(Qian_Lookup_Angios_TPL.selected$genus))){
#  mrca.node<-getMRCA(Qian_dropped,Qian_Lookup_Angios_TPL.selected[Qian_Lookup_Angios_TPL.selected$genus %in% unique(Qian_Lookup_Angios_TPL.selected$genus)[i],]$New_Species)
#  genus.clade<-extract.clade(Qian_dropped,mrca.node)
#  cat(unique(Qian_Lookup_Angios_TPL.selected$genus)[i],'\n')
#  tot_time<-max(node.age(genus.clade)$ages)
#  sampling.fraction<-Qian_Lookup_Angios_TPL.selected[Qian_Lookup_Angios_TPL.selected$genus %in% unique(Qian_Lookup_Angios_TPL.selected$genus)[i],]$genus.sampling.fraction[1]
#  f.lamb.cst<-function(t,y){y[1]}
#  f.lamb.exp<-function(t,y){y[1]*exp(y[2]*t)}
#  f.mu.cst.0<-function(t,y){0}
#  f.mu.cst<-function(t,y){y[1]}
#  f.mu.exp<-function(t,y){y[1]*exp(y[2]*t)}
#  lamb_par_init.exp<-c(0.05,0.01)
#  mu_par_init.exp<-c(0.01,0.01)
#  lamb_par_init.cst<-c(0.05)
#  mu_par_init.cst<-c(0.01)
#  mu_par_init.0<-c()
#  res.lambda.cst.mu.0<-fit_bd(genus.clade,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par_init.cst,mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE)
#  res.lambda.cst.mu.cst<-fit_bd(genus.clade,tot_time,f.lamb.cst,f.mu.cst,lamb_par_init.cst,mu_par_init.cst,f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE)
#  res.lambda.exp.mu.exp<-fit_bd(genus.clade,tot_time,f.lamb.exp,f.mu.exp,lamb_par_init.exp,mu_par_init.exp,f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE)
#  res.lambda.exp.mu.cst<-fit_bd(genus.clade,tot_time,f.lamb.exp,f.mu.cst,lamb_par_init.exp,mu_par_init.cst,f=sampling.fraction,expo.lamb=TRUE,cst.mu=TRUE)
#  res.lambda.cst.mu.exp<-fit_bd(genus.clade,tot_time,f.lamb.cst,f.mu.exp,lamb_par_init.cst,mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE)
#  aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.exp.mu.exp$aicc,res.lambda.exp.mu.cst$aicc,res.lambda.cst.mu.exp$aicc,res.lambda.cst.mu.cst$aicc)
#  names(aicc.vector)<-c('lambda.cst.mu.0','lambda.exp.mu.exp','lambda.exp.mu.cst','lambda.cst.mu.exp','lambda.cst.mu.cst')
#  aicc.weights<-sort(Weights(aicc.vector),decreasing=T)
#  if (names(aicc.weights)[1]=='lambda.cst.mu.0'){
#    lambda.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
#    mu.rpanda1<-0
#    ndr.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
#  }
#  if (names(aicc.weights)[1]=='lambda.exp.mu.exp'){
#    lambda.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]
#    mu.rpanda1<-res.lambda.exp.mu.exp$mu_par[1]
#    ndr.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]-res.lambda.exp.mu.exp$mu_par[1]
#  }
#  if (names(aicc.weights)[1]=='lambda.exp.mu.cst'){
#    lambda.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]
#    mu.rpanda1<-res.lambda.exp.mu.cst$mu_par[1]
#    ndr.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]-res.lambda.exp.mu.cst$mu_par[1]
#  }
#  if (names(aicc.weights)[1]=='lambda.cst.mu.exp'){
#    lambda.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]
#    mu.rpanda1<-res.lambda.cst.mu.exp$mu_par[1]
#    ndr.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]-res.lambda.cst.mu.exp$mu_par[1]
#  }
#  if (names(aicc.weights)[1]=='lambda.cst.mu.cst'){
#    lambda.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]
#    mu.rpanda1<-res.lambda.cst.mu.cst$mu_par[1]
#    ndr.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]-res.lambda.cst.mu.cst$mu_par[1]
#  }
#  aicc<-names(aicc.weights)[1]
#  cat(aicc,'\n')
#  rpanda.matrix<-rbind(rpanda.matrix,c(unique(Qian_Lookup_Angios_TPL.selected$genus)[i],lambda.rpanda1,mu.rpanda1,ndr.rpanda1,aicc))
#}
##add the ms info to the table
#rpanda.table<-as.data.frame(rpanda.matrix)

#colnames(rpanda.table)<-c('genus','lambda.rpanda','mu.rpanda','ndr.rpanda','aicc')
#Qian_Lookup_Angios_TPL.selected<-merge(Qian_Lookup_Angios_TPL.selected,rpanda.table)
#add genus seed weights
genusseed<-read.table('./output/tables/BAMM_trait_rates_seed.txt',header=T)
genusseed<-genusseed[,c(1,3)]
colnames(genusseed)<-c('genus','log.Seed.Weight.genus')
Qian_Lookup_Angios_TPL.selected<-merge(Qian_Lookup_Angios_TPL.selected,genusseed,all.x=TRUE)
#select data with seed mass
Qian_Lookup_Angios_TPL.selected<-na.omit(Qian_Lookup_Angios_TPL.selected)
#droptips in tree & table to one representative per genus
cat('preparing dataset for pgls','\n')
sample.matrix<-matrix(NA,nrow=0,ncol=0)
for (i in 1:length(unique(Qian_Lookup_Angios_TPL.selected$genus))){
  cat(unique(Qian_Lookup_Angios_TPL.selected$genus)[i],'\n')
  sample<-Qian_Lookup_Angios_TPL.selected[match(unique(Qian_Lookup_Angios_TPL.selected$genus)[i],Qian_Lookup_Angios_TPL.selected$genus),]
  sample.matrix<-rbind(sample.matrix,sample)
}
sample<-as.data.frame(sample.matrix)
Qian_dropped.selected<-drop.tip(Qian_dropped,setdiff(Qian_dropped$tip.label,sample$New_Species))
pgls.df<-data.frame(sample$New_Species,sample$log.Seed.Weight.genus,sample$ms0,sample$ms05,sample$ms09,stringsAsFactors = F)
colnames(pgls.df)<-c('New_Species','log.Seed.Weight.genus','ms0','ms05','ms09')
cat ('running pgls','\n')
pgls.df$ms09<-as.numeric(pgls.df$ms09)
pgls.df$ms05<-as.numeric(pgls.df$ms05)
pgls.df$ms0<-as.numeric(pgls.df$ms0)
pgls.object<-comparative.data(phy=Qian_dropped.selected,data=pgls.df,names='New_Species',vcv=TRUE)
pgls.seedMS05<-pgls(log10(ms05)~log.Seed.Weight.genus, data = pgls.object, lambda='ML')
cat('pgls of Magallon&Sanderson estimator with epsilon 0.5 ~ seed size','\n')
cat('slope = ',summary(pgls.seedMS05)$coefficients[2,1],', p-value = ',summary(pgls.seedMS05)$coefficients[2,4])


