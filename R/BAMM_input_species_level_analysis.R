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
Qian_Lookup_Angios_TPL<-Qian_Lookup_Angios_TPL[,c(c(1:6),24)]
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
#####calculating sampling fraction
#using family level incomplete sampling (another option is to assume genus level incomplete sampling, but assumes monophyly of genera...)
plant_lookup_version_current()
PlantLookup<-plant_lookup(include_counts = TRUE)
PlantLookup.2<-aggregate(number.of.species~family, data=PlantLookup, sum)
colnames(PlantLookup.2)[2]<-'number.of.species.family'
PlantLookup<-merge(PlantLookup, PlantLookup.2, by.x = 'family', by.y = 'family', all.x = TRUE)
PlantLookup<-unique(PlantLookup[,c(1,6)])
table.sampling<-merge(Kew.TaxLookUp.Qiandropped,PlantLookup,by='family',all.x=TRUE)
family.count.table<-as.data.frame(table(table.sampling$family),stringsAsFactors = F)
colnames(family.count.table)<-c('family','number.of.species.tree')
family.count.table<-family.count.table[family.count.table$number.of.species.tree>0,]
table.sampling<-merge(table.sampling,family.count.table,by='family',all.x=TRUE)
table.sampling$family.sampling.fraction<-table.sampling$number.of.species.tree/table.sampling$number.of.species.family
Family_BAMM_Sampling<-table.sampling[,c("Qian_dropped.Seed.tip.label","family","family.sampling.fraction")]
colnames(Family_BAMM_Sampling)[1]<-c('species')
#write output files
write.table(Family_BAMM_Sampling, file = './output/tables/BAMM_Species_FamilySamplingFractions.txt', row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
table.sampling_backbone<-nrow(table.sampling)/sum(PlantLookup$number.of.species.family)
write.table(table.sampling_backbone,file = './output/tables/BAMM_Species_backbone_sampling.txt',row.names=FALSE,col.names = FALSE,quote=FALSE)
#write trait file
write.table(table.sampling[,c(2,7)],file='./output/tables/BAMM_Species_Seeddata.txt',row.names=F,quote=F,col.names=F,sep='\t')

#checks for BAMM
is.ultrametric(tree)
min(tree$edge.length) 
summary(duplicated(tree$tip.label))

#write tree for BAMM to file
write.tree(tree,file='./output/trees/BAMM_Species_tree.tree')

###BAMM priors
setBAMMpriors(tree,outfile='./output/tables/BAMM_Species_diversification_Priors.txt')
setBAMMpriors(tree,traits='./output/tables/BAMM_Species_Seeddata.txt',outfile='./output/tables/BAMM_Species_trait_Priors.txt')

