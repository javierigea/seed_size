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
#####it selects monophyletic genera with 4 or more species with seed size data in the Zanne tree
#####and measures diversification with Magallon & Sanderson and fits BM with fitContinuous
cat('preparing dataset to run clade based analyses','\n')
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
cat ('filtering genera','\n')
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
ms.matrix<-matrix(NA,nrow=0,ncol=4)
for (i in 1:length(unique(Qian_Lookup_Angios_TPL.selected$genus))){
  #cat(unique(Qian_Lookup_Angios_TPL.selected$genus)[i],'\n')
  age<-as.numeric(Qian_Lookup_Angios_TPL.selected[match(unique(Qian_Lookup_Angios_TPL.selected$genus)[i],Qian_Lookup_Angios_TPL.selected$genus),]$mrca.age)
  species<-as.numeric(Qian_Lookup_Angios_TPL.selected[match(unique(Qian_Lookup_Angios_TPL.selected$genus)[i],Qian_Lookup_Angios_TPL.selected$genus),]$n.described.species)
  ms0<-bd.ms(time=age,n=species,epsilon=0)
  ms05<-bd.ms(time=age,n=species,epsilon=0.5)
  ms09<-bd.ms(time=age,n=species,epsilon=0.9)
  ms.matrix<-rbind(ms.matrix,c(unique(Qian_Lookup_Angios_TPL.selected$genus)[i],ms0,ms05,ms09))
}
#add the ms info to the table
ms.table<-as.data.frame(ms.matrix)
colnames(ms.table)<-c('genus','ms0','ms05','ms09')
Qian_Lookup_Angios_TPL.selected<-merge(Qian_Lookup_Angios_TPL.selected,ms.table)
#overlap seed and tree datasets
Kew.TaxLookUp<-read.table('./output/tables/KewSIDdata_taxonlookup.txt',header=T,quote='',sep='\t')
Seed_Data<-Kew.TaxLookUp$Log10.Seed.Weight
names(Seed_Data)<-Kew.TaxLookUp$Species
name.check.Seed<-name.check(Qian_dropped,Seed_Data)
Qian_dropped.Seed<-drop.tip(Qian_dropped,name.check.Seed$tree_not_data)
Kew.TaxLookUp.Qiandropped<-merge(data.frame(Qian_dropped.Seed$tip.label),Kew.TaxLookUp,by.x='Qian_dropped.Seed.tip.label',by.y='Species')
Kew.TaxLookUp.Qiandropped<-Kew.TaxLookUp.Qiandropped[,c(1,8)]
colnames(Kew.TaxLookUp.Qiandropped)<-c('Species','Log.Seed.Weight')
Qian_Lookup_Angios_TPL.selected<-merge(Qian_Lookup_Angios_TPL.selected,Kew.TaxLookUp.Qiandropped,by.x='New_Species',by.y='Species',all.x=TRUE)
Qian_Lookup_Angios_TPL.selected<-na.omit(Qian_Lookup_Angios_TPL.selected)
#count the number of species per genera with seed data
genus.species.seedmass.table<-table(Qian_Lookup_Angios_TPL.selected$genus)
genus.species.seedmass.table<-data.frame(genus.species.seedmass.table,stringsAsFactors = F)
colnames(genus.species.seedmass.table)<-c('genus','n.species.seedmass')
Qian_Lookup_Angios_TPL.selected<-merge(Qian_Lookup_Angios_TPL.selected,genus.species.seedmass.table)
#select genera with > 3 species
Qian_Lookup_Angios_TPL.selected<-Qian_Lookup_Angios_TPL.selected[Qian_Lookup_Angios_TPL.selected$n.species.seedmass>3,]
#run fitcontinuous on selected
cat ('measuring trait evolution','\n')
fitcontinuous.matrix<-matrix(NA,nrow=0,ncol=2)
for (i in 1:length(unique(Qian_Lookup_Angios_TPL.selected$genus))){
  #cat(unique(Qian_Lookup_Angios_TPL.selected$genus)[i],'\n')
  mrca.node<-getMRCA(Qian_dropped,Qian_Lookup_Angios_TPL.selected[Qian_Lookup_Angios_TPL.selected$genus %in% unique(Qian_Lookup_Angios_TPL.selected$genus)[i],]$New_Species)
  genus.clade<-extract.clade(Qian_dropped,mrca.node)
  trait.vector<-Qian_Lookup_Angios_TPL.selected[Qian_Lookup_Angios_TPL.selected$genus %in% unique(Qian_Lookup_Angios_TPL.selected$genus)[i],]$Log.Seed.Weight
  names(trait.vector)<-Qian_Lookup_Angios_TPL.selected[Qian_Lookup_Angios_TPL.selected$genus %in% unique(Qian_Lookup_Angios_TPL.selected$genus)[i],]$New_Species
  res.trait<-fitContinuous(genus.clade,trait.vector,model="BM")
  fitcontinuous.matrix<-rbind(fitcontinuous.matrix,c(unique(Qian_Lookup_Angios_TPL.selected$genus)[i],res.trait$opt$sigsq))
}
#add the traitevo info to the table
fitcontinuous.matrix<-as.data.frame(fitcontinuous.matrix,stringsAsFactors = F)
colnames(fitcontinuous.matrix)<-c('genus','sigsq')
Qian_Lookup_Angios_TPL.selected<-merge(Qian_Lookup_Angios_TPL.selected,fitcontinuous.matrix)
#droptips in tree & table to one representative per genus
cat('preparing dataset for pgls','\n')
sample.matrix<-matrix(NA,nrow=0,ncol=0)
for (i in 1:length(unique(Qian_Lookup_Angios_TPL.selected$genus))){
  #cat(unique(Qian_Lookup_Angios_TPL.selected$genus)[i],'\n')
  sample<-Qian_Lookup_Angios_TPL.selected[match(unique(Qian_Lookup_Angios_TPL.selected$genus)[i],Qian_Lookup_Angios_TPL.selected$genus),]
  sample.matrix<-rbind(sample.matrix,sample)
}
sample<-as.data.frame(sample.matrix)
Qian_dropped.selected<-drop.tip(Qian_dropped,setdiff(Qian_dropped$tip.label,sample$New_Species))
pgls.df<-data.frame(sample$New_Species,sample$Log.Seed.Weight,sample$ms0,sample$ms05,sample$ms09,sample$sigsq,stringsAsFactors = F)
cat ('running pgls','\n')
colnames(pgls.df)<-c('New_Species','log.Seed.Weight','ms0','ms05','ms09','sigsq')
pgls.df$log.Seed.Weight<-as.numeric(pgls.df$log.Seed.Weight)
pgls.df$ms09<-as.numeric(pgls.df$ms09)
pgls.df$ms05<-as.numeric(pgls.df$ms05)
pgls.df$ms0<-as.numeric(pgls.df$ms0)
pgls.df$sigsq<-as.numeric(pgls.df$sigsq)
pgls.object<-comparative.data(phy=Qian_dropped.selected,data=pgls.df,names='New_Species',vcv=TRUE)
pgls.seedMS05<-pgls(log10(ms05)~log.Seed.Weight, data = pgls.object, lambda='ML')
cat('pgls of Magallon&Sanderson estimator with epsilon 0.5 ~ seed size','\n')
cat('slope = ',summary(pgls.seedMS05)$coefficients[2,1],', p-value = ',summary(pgls.seedMS05)$coefficients[2,4])
pgls.seedrateMS05<-pgls(log10(ms05)~log10(sigsq), data = pgls.object, lambda='ML')
cat('pgls of Magallon&Sanderson estimator with epsilon 0.5 ~ seed size rate','\n')
cat('slope = ',summary(pgls.seedrateMS05)$coefficients[2,1],', p-value = ',summary(pgls.seedrateMS05)$coefficients[2,4])
