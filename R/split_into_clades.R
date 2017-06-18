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
library(phangorn)
library(MuMIn)
library(ape)
library(picante)
library(geiger)
library(adephylo)

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

alltips<-tree$tip.label
all.tips.length<-lapply(alltips,function(x) length(x))
all.tips.length.df<-as.data.frame(unlist(all.tips.length))
all.tips.length.df$node.number<-c(1:nrow(all.tips.length.df))
all.tips.length.df$node.number<-all.tips.length.df$node.number+length(tree$tip.label)
colnames(all.tips.length.df)[1]<-'n.descendant.tips'

node.numbers<-c(1:length(tree$node.label))
node.numbers<-node.numbers+length(tree$tip.label)
length.node.descendants<-sapply(node.numbers,function(x) length(unlist(Descendants(tree,x,type='tips'))))
length.node.descendants<-as.data.frame(length.node.descendants)
length.node.descendants$parental.node<-node.numbers
colnames(length.node.descendants)[1]<-'n.descendant.tips'

node.age(tree)->phy.age
cbind(phy.age$edge,phy.age$age, tree$edge.length)->BL.position
max(phy.age$age)-BL.position[,3]->dist.tip
cbind(BL.position,dist.tip)->BL.positions
BL.positions[,5]+BL.positions[,4]->ages
cbind(BL.positions,ages)->BL.positions
as.data.frame(BL.positions)->node.ages
names(node.ages)<-c("parental.node","daughter.node","dist.root","BL","dist.tip","mrca.age")
node.ages$node.label<-NA

node.labels.df<-as.data.frame(c(1:length(tree$node.label)),tree$node.label)
node.labels.df$node.label<-row.names(node.labels.df)
row.names(node.labels.df)<-NULL
colnames(node.labels.df)[1]<-c('node.number')
node.labels.df$node.number<-node.labels.df$node.number+length(tree$tip.label)
node.ages<-merge(node.labels.df,node.ages,by.x='node.number','parental.node')
colnames(node.ages)[1]<-'parental.node'
node.ages<-merge(node.ages,length.node.descendants,by.x='parental.node',by.y='parental.node',all.x=TRUE)

nodes.to.check<-unique(node.ages$parental.node)
#start with Lamiidae, it's 4126 species in a monophyletic clade

lamiidae.node<-node.ages[node.ages$node.label=='Lamiidae',]$parental.node[1]
chosen.nodes<-vector()
length.chosen.clades<-vector()

chosen.nodes[1]<-lamiidae.node
length.chosen.clades<-node.ages[node.ages$node.label=='Lamiidae',]$n.descendant.tips[1]

if(node.ages[node.ages$parental.node==Siblings(tree,chosen.nodes[1])]$n.descendant.tips)

parent.node<-node.ages[(node.ages$daughter.node==chosen.nodes[length(chosen.nodes)]),]$parental.node[1]
sister.node<-node.ages[(node.ages$parental.node==parent.node)&(node.ages$daughter.node!=chosen.nodes[length(chosen.nodes)]),]$daughter.node

while(sum(length.chosen.clades)<length(tree$tip.label)){
  parent.node<-node.ages[(node.ages$daughter.node==chosen.nodes[length(chosen.nodes)]),]$parent.node[1]
  parent.node.df<-node.ages[node.ages$daughter.node==chosen.nodes[length(chosen.nodes)],]$parent.node[1]
}

#this function creates a set of non inclusive non monophyletic clades of <size
#starting with lamiidae.node
nodenumber<-lamiidae.node<-node.ages[node.ages$node.label=='Lamiidae',]$parental.node[1]
size<-5000

recursive_collapse_clades<-function(nodenumber,tree,size,table){
  chosen.nodes<-vector()
  length.chosen.clades<-vector()
  collapsed.species<-vector()
  while(sum(length.chosen.clades)<length(tree$tip.label)){
    #get number of descendants of node
    length.descendants<-length(unlist(Descendants(tree,nodenumber,type='tips')))
    #if number of descendants of node > size, go to left daughter
    sister<-Siblings(tree,nodenumber)
    
    length.sister.descendants<-length(unlist(Descendants(tree,sister,type='tips')))
    if (length.descendants>size){
      nodenumber<-Children(tree,nodenumber)[1]
    }
    #check if sister is already in chosen.nodes vector
    #if it is, add "descendants of"cousins"
    if(length(which(chosen.nodes==sister))>0){
      
      
    
    #get number of descendants of sister node
    
    #if sum of descendants of node + sister > size
    #keep node as one clade
    if(length.descendants+length.sister.descendants>size){
      chosen.nodes<-c(chosen.nodes,nodenumber)
      length.chosen.clades<-c(length.chosen.clades,length.descendants)
      
      node.number<-sister
    }else{
      #if sum of descendants of node + sister < size, collapse (go to ancestor)
      nodenumber<-Ancestors(tree,nodenumber,type='parent')
    }
    
     
    }else{
      
    }
    
    
  }
  
  
  #check n of descendants of sister in table
  #if sum of descendants of node and sister < size, collapse and nodenumber is the ancestor node
  #if not, store nodenumber and size in vector and use sister node as node number
  
  
}





Commelinidae<-getDescendants(tree, node=which(tree$node.label=='Commelinidae')+length(tree$tip.label))
length(Commelinidae[Commelinidae<=length(tree$tip.label)])
Commelinidae.species<-tree$tip.label[Commelinidae[Commelinidae<=length(tree$tip.label)]]
# 3367 species
Monocotyledoneae<-getDesister.nodescendants(tree, node=which(tree$node.label=='Monocotyledoneae')+length(tree$tip.label))
length(Monocotyledoneae[Monocotyledoneae<length(tree$tip.label)])
Monocotyledoneae.species<-tree$tip.label[Monocotyledoneae[Monocotyledoneae<=length(tree$tip.label)]]

Monocotyledoneae.species.noCom<-Monocotyledoneae.species[!(Monocotyledoneae.species %in% Commelinidae.species)]
#3591 species

which(tree$node.label=='Commelinidae')+length(tree$tip.label)



###select clades < 5000 species
node.ages.5000<-node.ages[node.ages$n.descendant.tips<5000,]
node.ages.5000<-node.ages.6000[,c('parental.node','n.descendant.tips')]
node.ages.5000<-unique(node.ages.5000)
node.ages.5000<-node.ages.5000[order(-node.ages.5000$n.descendant.tips),]
selected.species<-vector()
selected.nodes<-vector()
selected.lengths<-vector()
for (i in 1:nrow(node.ages.5000)){
  descendant.species<-unlist(Descendants(tree,node.ages.5000[i,'parental.node'],type='tips'))
  if(length(intersect(descendant.species,selected.species))==0){
    selected.species<-c(selected.species,descendant.species)
    selected.nodes<-c(selected.nodes,node.ages.5000[i,'parental.node'])
    selected.lengths<-c(selected.lengths,length(descendant.species))
    
  }
  cat(sum(selected.lengths),'\n')
  if(sum(selected.lengths)==length(tree$tip.label)){
    cat('done','\n')
    break
  }
}
#get one species per clade
species.selected<-tree$tip.label[sapply(selected.nodes,function(x) unlist(Descendants(tree,x,type='tips'))[1])]
tree.clades<-drop.tip(tree,setdiff(tree$tip.label,species.selected))
#add the number of species in each clade
species.selected.length<-paste(species.selected,selected.lengths,sep='_')
order.vector<-sapply(species.selected[match(tree.clades$tip.label,species.selected)],function(x) grep(x,species.selected.length))
species.selected.length<-species.selected.length[order(order)]
tree.clades$tip.label<- species.selected.length[order]
plot(tree.clades)


###select clades < 6000 species
node.ages.6000<-node.ages[node.ages$n.descendant.tips<6000,]
node.ages.6000<-node.ages.6000[,c('parental.node','n.descendant.tips')]
node.ages.6000<-unique(node.ages.6000)
node.ages.6000<-node.ages.6000[order(-node.ages.6000$n.descendant.tips),]
selected.species<-vector()
selected.nodes<-vector()
selected.lengths<-vector()
#this looks through the whole tree
for (i in 1:nrow(node.ages.6000)){
  descendant.species<-unlist(Descendants(tree,node.ages.6000[i,'parental.node'],type='tips'))
  if(length(intersect(descendant.species,selected.species))==0){
    selected.species<-c(selected.species,descendant.species)
    selected.nodes<-c(selected.nodes,node.ages.6000[i,'parental.node'])
    selected.lengths<-c(selected.lengths,length(descendant.species))
    
  }
  cat(sum(selected.lengths),'\n')
  if(sum(selected.lengths)==length(tree$tip.label)){
    cat('done','\n')
    break
  }
}
#add Amborella 
amborella.node<-which(tree$tip.label=='Amborella_trichopoda')
selected.nodes.Ambo<-c(selected.nodes,amborella.node)
selected.lengths.Ambo<-c(selected.lengths,1)
#get one species per clade
species.selected<-tree$tip.label[sapply(selected.nodes.Ambo,function(x) unlist(Descendants(tree,x,type='tips'))[1])]
tree.clades<-drop.tip(tree,setdiff(tree$tip.label,species.selected))
#add the number of species in each clade
species.selected.length<-paste(species.selected,selected.lengths.Ambo,sep='_')
order.vector<-sapply(species.selected[match(tree.clades$tip.label,species.selected)],function(x) grep(x,species.selected.length))

tree.clades$tip.label<- species.selected.length[order.vector]
plot(tree.clades)

#check monophyly of clades
species.selected.all<-sapply(selected.nodes.Ambo,function(x) unlist(Descendants(tree,x,type='tips')))
monophyly.check<-sapply(species.selected.all,function(x) is.monophyletic(tree,tree$tip.label[x]))
#they're all monophyletic
#then extract the clades > 3000 spp
big.monophyletic.clades<-lapply(selected.nodes.Ambo[selected.lengths.Ambo>3000], function(x) extract.clade(tree,node=x))
#and extract the rest of the species with the backbone analysis
big.monophyletic.clade.species<-lapply(selected.nodes.Ambo[selected.lengths.Ambo>3000],function(x) tree$tip.label[unlist(Descendants(tree,x,'tips'))])
#keep the first species on each big monophyletic clade (not dropping those to get the full backbone)
monophyletic.clade.species.to.drop<-unlist(lapply(big.monophyletic.clade.species,function(x) x[-1]))
backbone.clade<-drop.tip(tree,monophyletic.clade.species.to.drop)
#double check that there's no duplicated species in big monophyletic clades
length(unique(unlist(sapply(big.monophyletic.clades,function(x) x$tip.label))))
length(unlist(sapply(big.monophyletic.clades,function(x) x$tip.label)))
#get sampling fraction for each clade
BAMM.clades.to.run<-big.monophyletic.clades
BAMM.clades.to.run[[7]]<-backbone.clade
length(unique(unlist(sapply(BAMM.clades.to.run,function(x) x$tip.label))))
#there's 6 duplicates (one for each clade)
length(unlist(sapply(BAMM.clades.to.run,function(x) x$tip.label)))
#get the dataframes with sampling fractions
#it's easy for the big monopohyletic clades
big.monophyletic.clades.sampling<-lapply(big.monophyletic.clades,function(x) table.sampling[(table.sampling$New_Species %in% x$tip.label),])
big.monophyletic.clades.sampling.BAMM<-lapply(big.monophyletic.clades.sampling,function(x) x[,c("Tipname","family","family.sampling.fraction")])
big.monophyletic.clades.sampling.BAMM<-lapply(big.monophyletic.clades.sampling.BAMM,function(x) {colnames(x)[1]<-'species';return(x)})

#length of tips in each monophyletic clade
unlist(sapply(big.monophyletic.clades,function(x)length(x$tip.label)))
unlist(sapply(big.monophyletic.clades.sampling,function(x)nrow(x)))
#calculate the backbone sampling for these as well (this will be useful for the backbone clade as well)
big.monophyletic.clades.backbone.sampling<-lapply(big.monophyletic.clades.sampling,function(x) unique(x[,c('family','number.of.species.family','number.of.species.tree')]))
big.monophyletic.clades.backbone.samplingfraction<-unlist(lapply(big.monophyletic.clades.backbone.sampling,function(x) sum(x$number.of.species.tree)/sum(x$number.of.species.family)))
#now calculate the sampling fractions for the backbone clade
backbone.clades.sampling<-table.sampling[(table.sampling$New_Species %in% backbone.clade$tip.label),]
backbone.clades.sampling.BAMM<-backbone.clades.sampling[,c("Tipname","family","family.sampling.fraction")]
colnames(backbone.clades.sampling.BAMM)[1]<-c('species')
#substitute the sampling in the species representative of each monophyletic clade
big.monophyletic.clade.species.representative<-unlist(sapply(big.monophyletic.clade.species,function(x)x[1]))
new.sampling.fractions.clades<-1/(unlist(lapply(big.monophyletic.clades.backbone.sampling,function(x) sum(x$number.of.species.family))))
#replace in dataframe
for (i in 1:length(big.monophyletic.clade.species.representative)){
  backbone.clades.sampling.BAMM[backbone.clades.sampling.BAMM$species==big.monophyletic.clade.species.representative[i],]$family.sampling.fraction<-new.sampling.fractions.clades[i]
}
BAMM.clades.sampling.df<-big.monophyletic.clades.sampling.BAMM
BAMM.clades.sampling.df[[7]]<-backbone.clades.sampling.BAMM
#check that lengths of tree & sampling df are identical
unlist(lapply(BAMM.clades.to.run,function(x) length(x$tip.label)))
unlist(lapply(BAMM.clades.sampling.df,function(x) nrow(x)))
#calculate the sampling fraction of the backbone clade analysis
BAMM.clades.sampling.backbone.fractions<-c(big.monophyletic.clades.backbone.samplingfraction,nrow(BAMM.clades.sampling.df[[7]])/sum(PlantLookup$number.of.species.family))
#write output for BAMM
#first do checks and runsetBAMMpriors
lapply(BAMM.clades.to.run,function(x){c(is.binary.tree(x),is.ultrametric(x),min(x$edge.length),summary(duplicated(x$tip.label)))})
#they're all ok
#now run setBAMMpriors and store in files
for (i in 1:length(BAMM.clades.to.run)){
  setBAMMpriors(BAMM.clades.to.run[[i]],outfile = paste('./output/trees/Zanne_clades/clade_',i,'_priors.txt',sep=''))
}
#write treefiles
for (i in 1:length(BAMM.clades.to.run)){
  write.tree(BAMM.clades.to.run[[i]],file = paste('./output/trees/Zanne_clades/clade_',i,'.tree',sep=''))
}
#write sampling fractions,first backbones
for (i in 1:length(BAMM.clades.sampling.df)){
  write.table(BAMM.clades.sampling.backbone.fractions[i],file = paste('./output/trees/Zanne_clades/clade_',i,'_sampling.txt',sep=''),quote=F,row.names=F,sep='\t',col.names = F)
}

#write sampling fractions,then table
for (i in 1:length(BAMM.clades.sampling.df)){
  write.table(BAMM.clades.sampling.df[[i]],file = paste('./output/trees/Zanne_clades/clade_',i,'_sampling.txt',sep=''),quote=F,row.names=F,sep='\t',col.names = F,append = T)
}

#had to add extinctionProbMax = 0.99999 to the control file of clade 7 (the backbone)
#explained here: https://groups.google.com/forum/#!topic/bamm-project/X4m5BQkl7ls




