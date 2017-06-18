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

library(ape)
library(geiger)
library(Taxonstand)
library(BAMMtools)
library(coda)
library(taxonlookup)
library(phytools)
library(diversitree)
library(phyndr)
library(gdata)

################################################################################
###this function runs Qian tree tip labels through TPL (iterates through TPL three times to circumvent weird errors with Taxonstand)
##treefile: './raw_data/Qian_tree.txt'
##output to ./output/tables/Qian_dataset_TPL.txt
################
run_tree_through_TPL<-function(treefile){
  Qian<-read.tree(treefile)
  Qian.names<-sub(Qian$tip.label,pattern='_',replacement=' ')
  results.list.Qian<-list()
  errors.Qian<-vector('numeric')
  #run through TPL and add to original table
  for (i in 1:length(Qian.names)){
    tpl.output<-try(TPL(Qian.names[i]))
    if(class(tpl.output)=='try-error'){
      errors.Qian<-c(errors.Qian,i)
      cat(i, 'error','\n')
      next
    }else{
      results.list.Qian[[i]]<-tpl.output
      cat(i,'\n')
    }
  }
  results.list.Qian.2<-list()
  errors.Qian.2<-vector('numeric')
  for (i in 1:length(errors.Qian)){
    tpl.output<-try(TPL(Qian.names[i]))
    if(class(tpl.output)=='try-error'){
      errors.Qian.2<-c(errors.Qian.2,i)
      cat(i, 'error','\n')
      next
    }else{
      results.list.Qian.2[[i]]<-tpl.output
      cat(i,'\n')
    }
  }
  for (i in 1:length(errors.Qian)){
    results.list.Qian[[errors.Qian[i]]]<-results.list.Qian.2[[i]]
  }
  results.list.Qian.df <- do.call("rbind", results.list.Qian)
  
  species.false<-paste(results.list.Qian.df[results.list.Qian.df$Plant.Name.Index=='FALSE',]$Genus,results.list.Qian.df[results.list.Qian.df$Plant.Name.Index=='FALSE',]$Species,sep=' ')
  results.list.Qian.df<-results.list.Qian.df[!results.list.Qian.df$Plant.Name.Index==FALSE,]
  results.species.false.Qian<-list()
  errors.species.false.Qian<-vector('numeric')
  for (i in 1:length(species.false)){
    tpl.output<-try(TPL(species.false[i]))
    if(class(tpl.output)=='try-error'){
      errors.species.false.Qian<-c(errors.species.false.Qian,i)
      cat(i, 'error','\n')
      next
    }else{
      results.species.false.Qian[[i]]<-tpl.output
      cat(i,'\n')
    }
  }
  results.list.Qian.df.false<-do.call("rbind", results.species.false.Qian)
  results.list.Qian.df<-rbind(results.list.Qian.df,results.list.Qian.df.false)
  nrow(results.list.Qian.df)
  table(results.list.Qian.df$Plant.Name.Index)
  for (i in 1:length(species.false)){
    tpl.output<-try(TPL(species.false[i]))
    if(class(tpl.output)=='try-error'){
      errors.species.false.Qian<-c(errors.species.false.Qian,i)
      cat(i, 'error','\n')
      next
    }else{
      results.species.false.Qian[[i]]<-tpl.output
      cat(i,'\n')
    }
  }
  results.list.Qian.df.false<-do.call("rbind", results.species.false.Qian)
  results.list.Qian.df<-rbind(results.list.Qian.df,results.list.Qian.df.false)
  nrow(results.list.Qian.df)
  table(results.list.Qian.df$Plant.Name.Index)
  drop.columns<-c('Authority')
  results.list.Qian.df<-results.list.Qian.df[,!colnames(results.list.Qian.df)%in%drop.columns]
  write.table(results.list.Qian.df,file='./output/tables/Qian_dataset_TPL.txt',sep='\t',quote=F,row.names=F)
}

################################################################################
###this function runs KEW SID through TPL (iterates through TPL three times to circumvent weird errors with Taxonstand)
##traitfile: './raw_data/AllKewdata.txt
##output to ./output/tables/kewsid_dataset_TPL.txt
################
run_kewSID_through_TPL<-function(traitfile){
  Kew<-read.table(traitfile, sep="\t", header = T)
  #remove duplicates
  Kew<-unique(Kew)
  kewsid.names<-paste(Kew$Genus,Kew$Species,sep=' ')
  results.list.kewsid<-list()
  errors.kewsid<-vector('numeric')
  #run through TPL and add to original table
  for (i in 1:length(kewsid.names)){
    tpl.output<-try(TPL(kewsid.names[i]))
    if(class(tpl.output)=='try-error'){
      errors.kewsid<-c(errors.kewsid,i)
      cat(i, 'error','\n')
      next
    }else{
      results.list.kewsid[[i]]<-tpl.output
      cat(i,'\n')
    }
  }
  results.list.kewsid.2<-list()
  errors.kewsid.2<-vector('numeric')
  for (i in 1:length(errors.kewsid)){
    tpl.output<-try(TPL(kewsid.names[i]))
    if(class(tpl.output)=='try-error'){
      errors.kewsid.2<-c(errors.kewsid.2,i)
      cat(i, 'error','\n')
      next
    }else{
      results.list.kewsid.2[[i]]<-tpl.output
      cat(i,'\n')
    }
  }
  for (i in 1:length(errors.kewsid)){
    results.list.kewsid[[errors.kewsid[i]]]<-results.list.kewsid.2[[i]]
  }
  results.list.kewsid.df <- do.call("rbind", results.list.kewsid)
  
  
  species.false<-paste(results.list.kewsid.df[results.list.kewsid.df$Plant.Name.Index=='FALSE',]$Genus,results.list.kewsid.df[results.list.kewsid.df$Plant.Name.Index=='FALSE',]$Species,sep=' ')
  results.list.kewsid.df<-results.list.kewsid.df[!results.list.kewsid.df$Plant.Name.Index==FALSE,]
  results.species.false.kewsid<-list()
  errors.species.false.kewsid<-vector('numeric')
  for (i in 1:length(species.false)){
    tpl.output<-try(TPL(species.false[i]))
    if(class(tpl.output)=='try-error'){
      errors.species.false.kewsid<-c(errors.species.false.kewsid,i)
      cat(i, 'error','\n')
      next
    }else{
      results.species.false.kewsid[[i]]<-tpl.output
      cat(i,'\n')
    }
  }
  results.list.kewsid.df.false<-do.call("rbind", results.species.false.kewsid)
  results.list.kewsid.df<-rbind(results.list.kewsid.df,results.list.kewsid.df.false)
  species.false<-paste(results.list.kewsid.df[results.list.kewsid.df$Plant.Name.Index=='FALSE',]$Genus,results.list.kewsid.df[results.list.kewsid.df$Plant.Name.Index=='FALSE',]$Species,sep=' ')
  results.list.kewsid.df<-results.list.kewsid.df[!results.list.kewsid.df$Plant.Name.Index==FALSE,]
  results.species.false.kewsid<-list()
  errors.species.false.kewsid<-vector('numeric')
  for (i in 1:length(species.false)){
    tpl.output<-try(TPL(species.false[i]))
    if(class(tpl.output)=='try-error'){
      errors.species.false.kewsid<-c(errors.species.false.kewsid,i)
      cat(i, 'error','\n')
      next
    }else{
      results.species.false.kewsid[[i]]<-tpl.output
      cat(i,'\n')
    }
  }
  results.list.kewsid.df.false<-do.call("rbind", results.species.false.kewsid)
  results.list.kewsid.df<-rbind(results.list.kewsid.df,results.list.kewsid.df.false)
  nrow(results.list.kewsid.df)
  table(results.list.kewsid.df$Plant.Name.Index)
  
  species.false<-paste(results.list.kewsid.df[results.list.kewsid.df$Plant.Name.Index=='FALSE',]$Genus,results.list.kewsid.df[results.list.kewsid.df$Plant.Name.Index=='FALSE',]$Species,sep=' ')
  results.list.kewsid.df<-results.list.kewsid.df[!results.list.kewsid.df$Plant.Name.Index==FALSE,]
  results.species.false.kewsid<-list()
  errors.species.false.kewsid<-vector('numeric')
  for (i in 1:length(species.false)){
    tpl.output<-try(TPL(species.false[i]))
    if(class(tpl.output)=='try-error'){
      errors.species.false.kewsid<-c(errors.species.false.kewsid,i)
      cat(i, 'error','\n')
      next
    }else{
      results.species.false.kewsid[[i]]<-tpl.output
      cat(i,'\n')
    }
  }
  results.list.kewsid.df.false<-do.call("rbind", results.species.false.kewsid)
  results.list.kewsid.df<-rbind(results.list.kewsid.df,results.list.kewsid.df.false)
  nrow(results.list.kewsid.df)
  table(results.list.kewsid.df$Plant.Name.Index)
  drop.columns<-c('Authority')
  #add seed weight info
  results.list.kewsid.df$Old.Species<-paste(results.list.kewsid.df$Genus,results.list.kewsid.df$Species,sep='_')
  Kew$Old.Species<-paste(Kew$Genus,Kew$Species,sep='_')
  Kew<-Kew[,c('Old.Species','SeedWeight_g')]
  results.list.kewsid.df2<-merge(results.list.kewsid.df,Kew,by.x='Old.Species',by.y='Old.Species')
  results.list.kewsid.df<-results.list.kewsid.df[,!colnames(results.list.kewsid.df)%in%drop.columns]
  results.list.kewsid.df<-unique(results.list.kewsid.df)
  write.table(results.list.kewsid.df,file='./output/tables/kewsid_dataset_TPL.txt',sep='\t',quote=F,row.names=F)
}

##########################
###this function runs lifecycle dataset through TPL (iterates through TPL three times to circumvent weird errors with Taxonstand)
##annualfile: './raw_data/annual_CvalueDB.txt
##perennialfile: './raw_data/perennial_CvalueDB.txt
##output to ./output/tables/lifecycle_dataset_TPL.txt
#######
run_lifecycle_through_TPL<-function(annualfile,perennialfile){
  annual<-read.tabler(annualfile,sep='\t',header=T)
  perennial<-read.tabler(perennialfile,sep='\t',header=T)
  lifecycle.table<-rbind(annual,perennial)
  #remove duplicates
  lifecycle.table<-unique(lifecycle.table)
  lifecycle.names<-paste(lifecycle.table$Genus,lifecycle.table$Species,sep=' ')
  results.list.lifecycle<-list()
  errors.lifecycle<-vector('numeric')
  for (i in 1:length(lifecycle.names)){
    tpl.output<-try(TPL(lifecycle.names[i]))
    if(class(tpl.output)=='try-error'){
      errors.lifecycle<-c(errors.lifecycle,i)
      cat(i, 'error','\n')
      next
    }else{
      results.list.lifecycle[[i]]<-tpl.output
      cat(i,'\n')
    }
  }
  results.list.lifecycle.2<-list()
  errors.lifecycle.2<-vector('numeric')
  for (i in 1:length(errors.lifecycle)){
    tpl.output<-try(TPL(lifecycle.names[i]))
    if(class(tpl.output)=='try-error'){
      errors.lifecycle.2<-c(errors.lifecycle.2,i)
      cat(i, 'error','\n')
      next
    }else{
      results.list.lifecycle.2[[i]]<-tpl.output
      cat(i,'\n')
    }
  }
  
  for (i in 1:length(errors.lifecycle)){
    results.list.lifecycle[[errors.lifecycle[i]]]<-results.list.lifecycle.2[[i]]
  }
  results.list.lifecycle.df <- do.call("rbind", results.list.lifecycle)
  #there are still "FALSE" entries in Plant.Index (these were warnings) to be sorted out
  
  species.false<-paste(results.list.lifecycle.df[results.list.lifecycle.df$Plant.Name.Index=='FALSE',]$Genus,results.list.lifecycle.df[results.list.lifecycle.df$Plant.Name.Index=='FALSE',]$Species,sep=' ')
  results.list.lifecycle.df<-results.list.lifecycle.df[!results.list.lifecycle.df$Plant.Name.Index==FALSE,]
  results.species.false.lifecycle<-list()
  errors.species.false.lifecycle<-vector('numeric')
  for (i in 1:length(species.false)){
    tpl.output<-try(TPL(species.false[i]))
    if(class(tpl.output)=='try-error'){
      errors.species.false.lifecycle<-c(errors.species.false.lifecycle,i)
      cat(i, 'error','\n')
      next
    }else{
      results.species.false.lifecycle[[i]]<-tpl.output
      cat(i,'\n')
    }
  }
  results.list.lifecycle.df.false<-do.call("rbind", results.species.false.lifecycle)
  results.list.lifecycle.df<-rbind(results.list.lifecycle.df,results.list.lifecycle.df.false)
  species.false<-paste(results.list.lifecycle.df[results.list.lifecycle.df$Plant.Name.Index=='FALSE',]$Genus,results.list.lifecycle.df[results.list.lifecycle.df$Plant.Name.Index=='FALSE',]$Species,sep=' ')
  results.list.lifecycle.df<-results.list.lifecycle.df[!results.list.lifecycle.df$Plant.Name.Index==FALSE,]
  results.species.false.lifecycle<-list()
  errors.species.false.lifecycle<-vector('numeric')
  for (i in 1:length(species.false)){
    tpl.output<-try(TPL(species.false[i]))
    if(class(tpl.output)=='try-error'){
      errors.species.false.lifecycle<-c(errors.species.false.lifecycle,i)
      cat(i, 'error','\n')
      next
    }else{
      results.species.false.lifecycle[[i]]<-tpl.output
      cat(i,'\n')
    }
  }
  results.list.lifecycle.df.false<-do.call("rbind", results.species.false.lifecycle)
  results.list.lifecycle.df<-rbind(results.list.lifecycle.df,results.list.lifecycle.df.false)
  nrow(results.list.lifecycle.df)
  table(results.list.lifecycle.df$Plant.Name.Index)
  
  species.false<-paste(results.list.lifecycle.df[results.list.lifecycle.df$Plant.Name.Index=='FALSE',]$Genus,results.list.lifecycle.df[results.list.lifecycle.df$Plant.Name.Index=='FALSE',]$Species,sep=' ')
  results.list.lifecycle.df<-results.list.lifecycle.df[!results.list.lifecycle.df$Plant.Name.Index==FALSE,]
  results.species.false.lifecycle<-list()
  errors.species.false.lifecycle<-vector('numeric')
  for (i in 1:length(species.false)){
    tpl.output<-try(TPL(species.false[i]))
    if(class(tpl.output)=='try-error'){
      errors.species.false.lifecycle<-c(errors.species.false.lifecycle,i)
      cat(i, 'error','\n')
      next
    }else{
      results.species.false.lifecycle[[i]]<-tpl.output
      cat(i,'\n')
    }
  }
  results.list.lifecycle.df.false<-do.call("rbind", results.species.false.lifecycle)
  results.list.lifecycle.df<-rbind(results.list.lifecycle.df,results.list.lifecycle.df.false)
  nrow(results.list.lifecycle.df)
  table(results.list.lifecycle.df$Plant.Name.Index)
  drop.columns<-c('Authority')
  results.list.lifecycle.df<-results.list.lifecycle.df[,!colnames(results.list.lifecycle.df)%in%drop.columns]
  write.table(results.list.lifecycle.df,file='./output/tables/lifecycle_dataset_TPL.txt',sep='\t',quote=F,row.names=F)
}

##########################
###this function runs lifecycle dataset through TPL (iterates through TPL three times to circumvent weird errors with Taxonstand)
##traitfile: './raw_data/GlobalWoodinessDatabase.csv
##output to ./output/tables/woodiness_dataset_TPL.txt
#######
run_woodiness_through_TPL<-function(traitfile){
  woodiness.table<-read.csv(traitfile)
  woodiness.names<-woodiness.table$gs
  woodiness.TPL<-TPL(woodiness.names)
  results.list.woodiness<-list()
  errors.woodiness<-vector('numeric')
  for (i in 1:length(woodiness.names)){
    tpl.output<-try(TPL(woodiness.names[i]))
    if(class(tpl.output)=='try-error'){
      errors.woodiness<-c(errors.woodiness,i)
      cat(i, 'error','\n')
      next
    }else{
      results.list.woodiness[[i]]<-tpl.output
      cat(i,'\n')
    }
  }
  results.list.woodiness.2<-list()
  errors.woodiness.2<-vector('numeric')
  for (i in 1:length(errors.woodiness)){
    tpl.output<-try(TPL(woodiness.names[i]))
    if(class(tpl.output)=='try-error'){
      errors.woodiness.2<-c(errors.woodiness.2,i)
      cat(i, 'error','\n')
      next
    }else{
      results.list.woodiness.2[[i]]<-tpl.output
      cat(i,'\n')
    }
  }
  for (i in 1:length(errors.woodiness)){
    results.list.woodiness[[errors.woodiness[i]]]<-results.list.woodiness.2[[i]]
  }
  results.list.woodiness.df <- do.call("rbind", results.list.woodiness)
  
  species.false<-paste(results.list.woodiness.df[results.list.woodiness.df$Plant.Name.Index=='FALSE',]$Genus,results.list.woodiness.df[results.list.woodiness.df$Plant.Name.Index=='FALSE',]$Species,sep=' ')
  results.list.woodiness.df<-results.list.woodiness.df[!results.list.woodiness.df$Plant.Name.Index==FALSE,]
  results.species.false.woodiness<-list()
  errors.species.false.woodiness<-vector('numeric')
  for (i in 1:length(species.false)){
    tpl.output<-try(TPL(species.false[i]))
    if(class(tpl.output)=='try-error'){
      errors.species.false.woodiness<-c(errors.species.false.woodiness,i)
      cat(i, 'error','\n')
      next
    }else{
      results.species.false.woodiness[[i]]<-tpl.output
      cat(i,'\n')
    }
  }
  results.list.woodiness.df.false<-do.call("rbind", results.species.false.woodiness)
  results.list.woodiness.df<-rbind(results.list.woodiness.df,results.list.woodiness.df.false)
  species.false<-paste(results.list.woodiness.df[results.list.woodiness.df$Plant.Name.Index=='FALSE',]$Genus,results.list.woodiness.df[results.list.woodiness.df$Plant.Name.Index=='FALSE',]$Species,sep=' ')
  results.list.woodiness.df<-results.list.woodiness.df[!results.list.woodiness.df$Plant.Name.Index==FALSE,]
  results.species.false.woodiness<-list()
  errors.species.false.woodiness<-vector('numeric')
  for (i in 1:length(species.false)){
    tpl.output<-try(TPL(species.false[i]))
    if(class(tpl.output)=='try-error'){
      errors.species.false.woodiness<-c(errors.species.false.woodiness,i)
      cat(i, 'error','\n')
      next
    }else{
      results.species.false.woodiness[[i]]<-tpl.output
      cat(i,'\n')
    }
  }
  results.list.woodiness.df.false<-do.call("rbind", results.species.false.woodiness)
  results.list.woodiness.df<-rbind(results.list.woodiness.df,results.list.woodiness.df.false)
  nrow(results.list.woodiness.df)
  table(results.list.woodiness.df$Plant.Name.Index)
  
  species.false<-paste(results.list.woodiness.df[results.list.woodiness.df$Plant.Name.Index=='FALSE',]$Genus,results.list.woodiness.df[results.list.woodiness.df$Plant.Name.Index=='FALSE',]$Species,sep=' ')
  results.list.woodiness.df<-results.list.woodiness.df[!results.list.woodiness.df$Plant.Name.Index==FALSE,]
  results.species.false.woodiness<-list()
  errors.species.false.woodiness<-vector('numeric')
  for (i in 1:length(species.false)){
    tpl.output<-try(TPL(species.false[i]))
    if(class(tpl.output)=='try-error'){
      errors.species.false.woodiness<-c(errors.species.false.woodiness,i)
      cat(i, 'error','\n')
      next
    }else{
      results.species.false.woodiness[[i]]<-tpl.output
      cat(i,'\n')
    }
  }
  results.list.woodiness.df.false<-do.call("rbind", results.species.false.woodiness)
  results.list.woodiness.df<-rbind(results.list.woodiness.df,results.list.woodiness.df.false)
  nrow(results.list.woodiness.df)
  table(results.list.woodiness.df$Plant.Name.Index)
  drop.columns<-c('Authority')
  results.list.woodiness.df<-results.list.woodiness.df[,!colnames(results.list.woodiness.df)%in%drop.columns]
  write.table(results.list.woodiness.df,file='./output/tables/woodiness_dataset_TPL.txt',sep='\t',quote=F,row.names=F)
  
}

