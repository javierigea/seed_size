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
#script to collate Zanne tree with phenotypic traits individually
library(Taxonstand)
library(taxonlookup)
library(ape)
library(BAMMtools)


#####prepare_height_dataset includes TPL (as opposed to the rest of the traits where this is in a separate function)
prepare_height_dataset<-function(){
  table<-read.table('./raw_data/height_pruned_header.txt',header=T,sep='\t',fill=TRUE,stringsAsFactors = F)
  table<-table[,c(5,7,10,15,16,25)]
  table<-table[table$TraitID==18,]
  table$SIunit<-NA
  #drop records with errors in measure field (no values)
  table$SIunit.test<-sub(table$OrigValueStr,pattern=',',replacement='')
  table$SIunit.test<-as.numeric(table$SIunit.test)
  nas<-which(is.na(table$SIunit.test))
  table<-table[-nas,]
  table<-table[,-ncol(table)]
  #change all units to m (SI Unit)
  table$SIunit<-as.numeric(table$SIunit)
  table[table$OrigUnitStr=='cm',]$SIunit<-as.numeric(table[table$OrigUnitStr=='cm',]$OrigValueStr)/100
  table[table$OrigUnitStr=='m',]$SIunit<-as.numeric(table[table$OrigUnitStr=='m',]$OrigValueStr)
  table[table$OrigUnitStr=='feet',]$SIunit<-as.numeric(table[table$OrigUnitStr=='feet',]$OrigValueStr)*0.3048
  #drop records where ErrorRisk>4: see Try_Data_Release_Notes
  table<-table[table$ErrorRisk<4,]
  species.median<-aggregate(SIunit~AccSpeciesName,data=table,median)
  colnames(species.median)<-c('species','median.Plant.Height')
  species.sd<-aggregate(SIunit~AccSpeciesName,data=table,sd)
  colnames(species.sd)<-c('species','sd.Plant.Height')
  species.COV<-aggregate(SIunit~AccSpeciesName,data=table,function(x)sqrt(var(x))/mean(x))
  colnames(species.median)<-c('species','median.Plant.Height')
  
  TPL.plantheight<-lapply(species.median$species,function(x)try(TPL(x)))
  classes<-unlist(lapply(TPL.plantheight,function(x)class(x)))
  errors<-which(classes=='try-error')
  if(length(errors)>0){
    #repeat TPL with errors, there are problems with parentheses
    TPL.plantheight.errors<-lapply(sub(species.median[errors,]$species,pattern=' \\(.+',replacement=''),function(x)try(TPL(x)))
    #get them into TPL.plantheight list
    TPL.plantheight[errors]<-TPL.plantheight.errors
  }
  #TPL results to dataframe
  TPL.df<- do.call("rbind", TPL.plantheight)
  TPL.df<-TPL.df[,c('Genus','Species','Plant.Name.Index','New.Genus','New.Hybrid.marker','New.Species')]
  #remove hybrids
  TPL.df<-TPL.df[-which(TPL.df$New.Hybrid.marker=='×'),]
  TPL.species.df<-as.data.frame(cbind(paste(TPL.df$Genus,TPL.df$Species,sep='_'),paste(TPL.df$New.Genus,TPL.df$New.Species,sep='_'),TPL.df$Plant.Name.Index))
  
  colnames(TPL.species.df)<-c('OldSpecies','NewSpecies','Plant.Name.Index')
  TPL.species.df$OldSpecies<-as.character(TPL.species.df$OldSpecies)
  TPL.species.df$NewSpecies<-as.character(TPL.species.df$NewSpecies)
  #group measures by taxon names post TPL
  #add TPL info to table
  TPL.species.df$OldSpecies<-sub(TPL.species.df$OldSpecies,pattern='_',replacement=' ')
  table.TPL<-merge(table,TPL.species.df,by.x='AccSpeciesName',by.y='OldSpecies',all.x=TRUE)
  
  TPL.species.median<-aggregate(SIunit~NewSpecies,data=table.TPL,median)
  colnames(TPL.species.median)<-c('species','median.Plant.Height')
  TPL.species.sd<-aggregate(SIunit~NewSpecies,data=table.TPL,sd)
  colnames(TPL.species.sd)<-c('species','sd.Plant.Height')
  
  #run through taxonlookup
  TPL.species.Lookup<-lookup_table(unique(TPL.species.df$NewSpecies), by_species = TRUE)
  TPL.species.Lookup$species<-row.names(TPL.species.Lookup)
  row.names(TPL.species.Lookup)<-NULL
  #get angiosperms only
  TPL.species.Lookup.angios<-TPL.species.Lookup[TPL.species.Lookup$group=='Angiosperms',]
  TPL.species.angios.height<-merge(TPL.species.Lookup.angios,TPL.species.median,all.x=TRUE)
  TPL.species.angios.height<-merge(TPL.species.angios.height,TPL.species.sd,all.x=TRUE)
  plot(TPL.species.angios.height$median.Plant.Height,TPL.species.angios.height$sd,xlab='median.plant.height',ylab='sd.plant.height',main='TPL.dataset')
  
  tree<-read.tree('./output/trees/BAMM_Species_tree_noseed.tree')
  height.tree<-drop.tip(tree,setdiff(tree$tip.label,TPL.species.angios.height$species))
  TPL.species.angios.height.tree<-TPL.species.angios.height[TPL.species.angios.height$species%in%height.tree$tip.label,]
  
  plot(TPL.species.angios.height.tree$median.Plant.Height,TPL.species.angios.height.tree$sd,xlab='median.plant.height',ylab='sd.plant.height',main='TPL.tree.dataset')
  write.tree(height.tree,file='./output/trees/BAMM_Species_tree_height.tree')
  height.data<-TPL.species.angios.height.tree[,c('species','median.Plant.Height')]
  height.data$median.Plant.Height<-log10(height.data$median.Plant.Height)
  write.table(height.data,file='./output/tables/BAMM_Species_tree_height_table.txt',sep='\t',quote=F,row.names=F,col.names = F)
  setBAMMpriors(height.tree,traits='./output/tables/BAMM_Species_tree_height_table.txt',outfile='./height_priors.txt')
}
  


#######prepare_woodiness_dataset<-function(){
#######  #############WOODINESS
#######  woodiness<-read.csv('./raw_data/GlobalWoodinessDatabase.csv',stringsAsFactors = F)
#######  #just using the Zanne et al coding here
#######  woodiness$woody.state<-NA
#######  woodiness[woodiness$woodiness=='H',]$woody.state<-0
#######  woodiness[woodiness$woodiness=='W',]$woody.state<-1
#######  #dropping the "variable species"
#######  woodiness<-na.omit(woodiness)
#######  woodiness$Merged<-sub(woodiness$gs,pattern=' ',replacement='_')
#######  woodiness<-woodiness[,c("Merged","woody.state")]
#######  woodiness.TPL<-read.table('./output/tables/woodiness_dataset_TPL.txt',header=T,sep='\t')
#######  woodiness.TPL$Old.Species<-paste(woodiness.TPL$Genus,woodiness.TPL$Species,sep='_')
#######  woodiness.TPL<-merge(woodiness,woodiness.TPL,by.x='Merged',by.y='Old.Species')
#######  colnames(woodiness.TPL)[c(1,2)]<-c('Old.Species','woodiness.state')
#######  woodiness.TPL$Merged<-paste(woodiness.TPL$New.Genus,woodiness.TPL$New.Species,sep='_')
#######  woodiness.TPL<-woodiness.TPL[!duplicated(woodiness.TPL$Merged),]
#######  woodiness.TPL<-woodiness.TPL[,c('Merged','woodiness.state')]
#######  tree<-read.tree('./output/trees/BAMM_Species_tree_noseed.tree')
#######  woodiness.tree<-drop.tip(tree,setdiff(tree$tip.label,woodiness.TPL$Merged))
#######  write.tree(woodiness.tree,file='./output/trees/BAMM_Species_tree_woodiness.tree')
#######  woodiness.data<-woodiness.TPL[woodiness.TPL$Merged%in%woodiness.tree$tip.label,]
#######  write.table(woodiness.data,file='./output/tables/BAMM_Species_tree_woodiness_table.txt',sep='\t',quote=F,row.names=F,colna)
#######}  

prepare_lifecycle_Cvalue_dataset<-function(){
  #############LIFE CYCLE
  lifecycle.TPL<-read.table('./output/tables/lifecycle_dataset_TPL.txt',header=T,sep='\t')
  annuals<-read.csv('./raw_data/annual_CvalueDB.txt', header = T, sep = '\t')
  perennials<-read.csv('./raw_data/perennial_CvalueDB.txt',header = T, sep = '\t')
  #annuals
  annuals$Merged<-paste(annuals$Genus, annuals$Species, sep = '_')
  perennials$Merged<-paste(perennials$Genus, perennials$Species, sep = '_')
  lifecycle.TPL$Old.Species<-paste(lifecycle.TPL$Genus,lifecycle.TPL$Species,sep='_')
  annuals<-annuals[,c(6,3)]
  annuals.TPL<-merge(annuals,lifecycle.TPL,by.x='Merged',by.y='Old.Species')
  colnames(annuals.TPL)[c(1,2)]<-c('Old.Species','C.value')
  annuals.TPL$Merged<-paste(annuals.TPL$New.Genus,annuals.TPL$New.Species,sep='_')
  annuals.TPL<-annuals.TPL[!duplicated(annuals.TPL$Merged),]
  annuals.TPL<-annuals.TPL[,c('Merged','C.value')]
  annuals.TPL$cycle<-'annual'
  #perennials
  perennials$Merged<-paste(perennials$Genus, perennials$Species, sep = '_')
  perennials$Merged<-paste(perennials$Genus, perennials$Species, sep = '_')
  lifecycle.TPL$Old.Species<-paste(lifecycle.TPL$Genus,lifecycle.TPL$Species,sep='_')
  perennials<-perennials[,c(6,3)]
  perennials.TPL<-merge(perennials,lifecycle.TPL,by.x='Merged',by.y='Old.Species')
  colnames(perennials.TPL)[c(1,2)]<-c('Old.Species','C.value')
  perennials.TPL$Merged<-paste(perennials.TPL$New.Genus,perennials.TPL$New.Species,sep='_')
  perennials.TPL<-perennials.TPL[!duplicated(perennials.TPL$Merged),]
  perennials.TPL<-perennials.TPL[,c('Merged','C.value')]
  perennials.TPL$cycle<-'perennial'
  lifecycle.TPL<-rbind(annuals.TPL,perennials.TPL)
  #sort out duplicates
  lifecycle.TPL<-lifecycle.TPL[!duplicated(lifecycle.TPL$Merged),]
  lifecycle.TPL<-unique(lifecycle.TPL)
  tree<-read.tree('./output/trees/BAMM_Species_tree_noseed.tree')
  lifecycle.tree<-drop.tip(tree,setdiff(tree$tip.label,lifecycle.TPL$Merged))
  write.tree(lifecycle.tree,file='./output/trees/BAMM_Species_tree_lifecycle.tree')
  write.tree(lifecycle.tree,file='./output/trees/BAMM_Species_tree_Cvalue.tree')
  lifecycle.data<-lifecycle.TPL[lifecycle.TPL$Merged%in%lifecycle.tree$tip.label,c('Merged','cycle')]
  write.table(lifecycle.data,file='./output/tables/BAMM_Species_tree_lifecycle_table.txt',sep='\t',quote=F,row.names=F,col.names=F)
  
  cvalue.data<-lifecycle.TPL[lifecycle.TPL$Merged%in%lifecycle.tree$tip.label,c('Merged','C.value')]
  cvalue.data$C.value<-log10(cvalue.data$C.value)
  write.table(cvalue.data,file='./output/tables/BAMM_Species_tree_Cvalue_table.txt',sep='\t',quote=F,row.names=F,col.names = F)
  setBAMMpriors(lifecycle.tree,traits='./output/tables/BAMM_Species_tree_Cvalue_table.txt',outfile='./Cvalue_priors.txt')
}  



######################################################
#function to collate tree + seed + lifecycle + Cvalue + woodiness datasets
#load the KewSID + Qian tree dataset
################prepare_lifecycle_Cvalue_woodiness_dataset<-function(){
################  kewsid.TPL<-read.table('./output/tables/BAMM_Species_Seeddata.txt',sep='\t',header=F)
################  colnames(kewsid.TPL)<-c('Merged','log10.SeedWeight')
################  tree<-read.tree('./output/trees/BAMM_Species_tree.tree')
################  #############LIFE CYCLE
################  lifecycle.TPL<-read.table('./output/tables/lifecycle_dataset_TPL.txt',header=T,sep='\t')
################  annuals<-read.csv('./raw_data/annual_CvalueDB.txt', header = T, sep = '\t')
################  perennials<-read.csv('./raw_data/perennial_CvalueDB.txt',header = T, sep = '\t')
################  #annuals
################  annuals$Merged<-paste(annuals$Genus, annuals$Species, sep = '_')
################  perennials$Merged<-paste(perennials$Genus, perennials$Species, sep = '_')
################  lifecycle.TPL$Old.Species<-paste(lifecycle.TPL$Genus,lifecycle.TPL$Species,sep='_')
################  annuals<-annuals[,c(6,3)]
################  annuals.TPL<-merge(annuals,lifecycle.TPL,by.x='Merged',by.y='Old.Species')
################  colnames(annuals.TPL)[c(1,2)]<-c('Old.Species','C.value')
################  annuals.TPL$Merged<-paste(annuals.TPL$New.Genus,annuals.TPL$New.Species,sep='_')
################  annuals.TPL<-annuals.TPL[!duplicated(annuals.TPL$Merged),]
################  annuals.TPL<-annuals.TPL[,c('Merged','C.value')]
################  annuals.TPL$cycle<-'annual'
################  #perennials
################  perennials$Merged<-paste(perennials$Genus, perennials$Species, sep = '_')
################  perennials$Merged<-paste(perennials$Genus, perennials$Species, sep = '_')
################  lifecycle.TPL$Old.Species<-paste(lifecycle.TPL$Genus,lifecycle.TPL$Species,sep='_')
################  perennials<-perennials[,c(6,3)]
################  perennials.TPL<-merge(perennials,lifecycle.TPL,by.x='Merged',by.y='Old.Species')
################  colnames(perennials.TPL)[c(1,2)]<-c('Old.Species','C.value')
################  perennials.TPL$Merged<-paste(perennials.TPL$New.Genus,perennials.TPL$New.Species,sep='_')
################  perennials.TPL<-perennials.TPL[!duplicated(perennials.TPL$Merged),]
################  perennials.TPL<-perennials.TPL[,c('Merged','C.value')]
################  perennials.TPL$cycle<-'perennial'
################  lifecycle.TPL<-rbind(annuals.TPL,perennials.TPL)
################  #sort out duplicates
################  lifecycle.TPL<-lifecycle.TPL[!duplicated(lifecycle.TPL$Merged),]
################  lifecycle.TPL<-unique(lifecycle.TPL)
################  
################  #############WOODINESS
################  woodiness<-read.csv('./raw_data/GlobalWoodinessDatabase.csv',stringsAsFactors = F)
################  #just using the Zanne et al coding here
################  woodiness$woody.state<-NA
################  woodiness[woodiness$woodiness=='H',]$woody.state<-0
################  woodiness[woodiness$woodiness=='W',]$woody.state<-1
################  #dropping the "variable species"
################  woodiness<-na.omit(woodiness)
################  woodiness$Merged<-sub(woodiness$gs,pattern=' ',replacement='_')
################  woodiness<-woodiness[,c("Merged","woody.state")]
################  woodiness.TPL<-read.table('./output/tables/woodiness_dataset_TPL.txt',header=T,sep='\t')
################  woodiness.TPL$Old.Species<-paste(woodiness.TPL$Genus,woodiness.TPL$Species,sep='_')
################  woodiness.TPL<-merge(woodiness,woodiness.TPL,by.x='Merged',by.y='Old.Species')
################  colnames(woodiness.TPL)[c(1,2)]<-c('Old.Species','woodiness.state')
################  woodiness.TPL$Merged<-paste(woodiness.TPL$New.Genus,woodiness.TPL$New.Species,sep='_')
################  woodiness.TPL<-woodiness.TPL[!duplicated(woodiness.TPL$Merged),]
################  woodiness.TPL<-woodiness.TPL[,c('Merged','woodiness.state')]
################  
################  #merge all datasets together
################  kewsid.TPL.lifecycle<-merge(kewsid.TPL,lifecycle.TPL)
################  kewsid.TPL.lifecycle.woodiness<-merge(kewsid.TPL.lifecycle,woodiness.TPL)
################  colnames(kewsid.TPL.lifecycle.woodiness)[1]<-c('Species')
################  write.table(kewsid.TPL.lifecycle.woodiness,'./output/tables/kewsid_lifecycle_woodiness_TPL.txt',sep='\t',quote=F,row.names=F)
################  
################}