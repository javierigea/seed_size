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

######################################################
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
#remove "_sp" tips
remove_sp_tips<-Qian_Lookup_Angios_TPL[grep('_sp$',Qian_Lookup_Angios_TPL$Old_Species),]$Old_Species
Qian_dropped<-drop.tip(Qian_dropped,remove_sp_tips)
Qian_Lookup_Angios_TPL<-Qian_Lookup_Angios_TPL[-grep('_sp$',Qian_Lookup_Angios_TPL$New_Species),]
Qiandf<-data.frame(Qian_dropped$tip.label)
colnames(Qiandf)<-'Tipname'
#replace tip names in Qian Tree with TPL+taxonlookup alternative
Qian_Lookup_Angios_TPL<-merge(Qiandf,Qian_Lookup_Angios_TPL,by.x='Tipname',by.y='Old_Species')
Qian_Lookup_Angios_TPL$Tipname<-as.character(Qian_Lookup_Angios_TPL$Tipname)
Qian_Lookup_Angios_TPL<-Qian_Lookup_Angios_TPL[match(Qian_dropped$tip.label,Qian_Lookup_Angios_TPL$Tipname),]
Qian_dropped$tip.label<-Qian_Lookup_Angios_TPL$New_Species
#overlap seed and tree datasets
#kew.sid<-read.table('./output/tables/KewSIDdata_taxonlookup.txt',header=T,quote='',sep='\t')
#Seed_Data<-kew.sid$Log10.Seed.Weight
#names(Seed_Data)<-as.character(kew.sid$Species)
#name.check.Seed<-name.check(Qian_dropped,Seed_Data)
#Qian_dropped.Seed<-drop.tip(Qian_dropped,name.check.Seed$tree_not_data)
#kew.sid.Qiandropped<-merge(data.frame(Qian_dropped.Seed$tip.label),kew.sid,by.x='Qian_dropped.Seed.tip.label',by.y='Species')
tree<-Qian_dropped

#####calculating sampling fraction
#using family level incomplete sampling (another option is to assume genus level incomplete sampling, but assumes monophyly of genera...)
plant_lookup_version_current()
PlantLookup<-plant_lookup(include_counts = TRUE)
PlantLookup.2<-aggregate(number.of.species~family, data=PlantLookup, sum)
colnames(PlantLookup.2)[2]<-'number.of.species.family'
PlantLookup<-merge(PlantLookup, PlantLookup.2, by.x = 'family', by.y = 'family', all.x = TRUE)
PlantLookup<-unique(PlantLookup[,c(1,6)])
table.sampling<-merge(Qian_Lookup_Angios_TPL,PlantLookup,by='family',all.x=TRUE)
family.count.table<-as.data.frame(table(table.sampling$family),stringsAsFactors = F)
colnames(family.count.table)<-c('family','number.of.species.tree')
family.count.table<-family.count.table[family.count.table$number.of.species.tree>0,]
table.sampling<-merge(table.sampling,family.count.table,by='family',all.x=TRUE)
table.sampling$family.sampling.fraction<-table.sampling$number.of.species.tree/table.sampling$number.of.species.family
#adding info on taxonomic discordance (when there are more tips in the tree than accepted species in a family according to taxonlookup)
#it's only 1 family (Nepenthaceae) with 6 described species
table.sampling$taxonomic.discordance<-NA
table.sampling[table.sampling$number.of.species.family<table.sampling$number.of.species.tree,]$taxonomic.discordance<-1
table.sampling[table.sampling$number.of.species.family>=table.sampling$number.of.species.tree,]$taxonomic.discordance<-0
#remove taxonomic discordant species
#remove species with taxonomic discordance
table.sampling<-table.sampling[!table.sampling$taxonomic.discordance==1,]
tree<-drop.tip(tree,setdiff(tree$tip.label,table.sampling$New_Species))

Family_BAMM_Sampling<-table.sampling[,c("Tipname","family","family.sampling.fraction")]
colnames(Family_BAMM_Sampling)[1]<-c('species')

#write output files
write.table(Family_BAMM_Sampling, file = './output/tables/BAMM_Species_FamilySamplingFractions_noseed.txt', row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
table.sampling_backbone<-nrow(table.sampling)/sum(PlantLookup$number.of.species.family)
write.table(table.sampling_backbone,file = './output/tables/BAMM_Species_backbone_sampling_noseed.txt',row.names=FALSE,col.names = FALSE,quote=FALSE)
#write trait file
write.table(table.sampling[,c('Qian_dropped.Seed.tip.label','Log10.Seed.Weight')],file='./output/tables/BAMM_Species_Seeddata.txt',row.names=F,quote=F,col.names=F,sep='\t')


#checks for BAMM
is.binary.tree(tree)
is.ultrametric(tree)
min(tree$edge.length) 
summary(duplicated(tree$tip.label))

#write tree for BAMM to file
write.tree(tree,file='./output/trees/BAMM_Species_tree_noseed.tree')

###BAMM priors
setBAMMpriors(tree,outfile='./output/tables/BAMM_Species_diversification_Priors_noseed.txt')
setBAMMpriors(tree,traits='./output/tables/BAMM_Species_Seeddata.txt',outfile='./output/tables/BAMM_Species_trait_Priors.txt')

########
library(phytools)
trees<-getCladesofSize(tree,5000)
