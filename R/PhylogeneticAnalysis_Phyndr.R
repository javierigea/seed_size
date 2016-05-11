library(ape)
library(Taxonstand)
library(taxonlookup)
library(phytools)
library(diversitree)
library(geiger)
library(BAMMtools)
library(nlme)
library(phyndr)
library(MonoPhy)
#Qian and Jin created an updated version of Zanne tree - 'Appendix S3-v2' or 'QianTree.txt'
#reference doi: 10.1093/jpe/rtv047
Qian<-read.tree('./raw_data/QianTree.txt')
Qian_names<-Qian$tip.label
Qian_names<-sub(Qian_names,pattern='_',replacement=' ')
#run TPL for species names
Qian_TPL<-TPL(Qian_names)
Qian_TPL[,'Merged']<-paste(Qian_TPL$New.Genus,Qian_TPL$New.Species, sep = '_')
write.table(Qian_TPL,file = './output/tables/Qian_TPL.txt', quote = FALSE, row.names = FALSE, sep = '\t')
#run through taxonlookup
Qian_TPLunique<-Qian_TPL[!duplicated(Qian_TPL$Merged),]
Qian_Lookup<-lookup_table(Qian_TPLunique$Merged,by_species = TRUE)
#get table with angiosperms only
Qian_Lookup_Angios<-subset(Qian_Lookup,Qian_Lookup$group=='Angiosperms')
Qian_Lookup_Angios$Fullspecies<-row.names(Qian_Lookup_Angios)
row.names(Qian_Lookup_Angios)<-NULL
Qian_Lookup_Angios_TPL<-merge(Qian_Lookup_Angios,Qian_TPLunique,by.x='Fullspecies',by.y='Merged')
colnames(Qian_Lookup_Angios_TPL)[1]<-'New_Species'
Qian_Lookup_Angios_TPL$Old_Species<-paste(Qian_Lookup_Angios_TPL$Genus,Qian_Lookup_Angios_TPL$Species,sep='_')
Qian_Lookup_Angios_TPL<-Qian_Lookup_Angios_TPL[,c(c(1:5),23)]
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
#run phyndr
phyndrgenus<-phyndr_genus(Qian_dropped,names(Seed_Data))
phyndr.sample<-phyndr_sample(phyndrgenus)
#checking for monophyly
monoreport<-AssessMonophyly(phyndr.sample)
outliertips<-GetOutlierTips(monoreport)
outliertips<-unlist(outliertips)
outliertips<-unname(outliertips)
#drop the outliers from tree
QianSeedTreePhyndr_nooutliers<-drop.tip(phyndr.sample,outliertips)
#(still getting around 10% non monophyletic genera)
monoreport_nooutliers<-AssessMonophyly(QianSeedTreePhyndr_nooutliers)
GetSummaryMonophyly(monoreport_nooutliers)
#genus level tree,get genera only
generaPhyndr<-unique(sapply(strsplit(as.character(QianSeedTreePhyndr_nooutliers$tip.label),'_'),function(x) x[1]))
generaPhyndr.2<-sapply(generaPhyndr, function(x) x<- paste(x,'_',sep=''))
ii<-sapply(generaPhyndr.2,function(x,y) grep(x,y,fixed=TRUE)[1],y=QianSeedTreePhyndr_nooutliers$tip.label)
#drop all but one of every genera
QianSeedTreePhyndr_nooutliers_GenusTree<-drop.tip(QianSeedTreePhyndr_nooutliers,setdiff(QianSeedTreePhyndr_nooutliers$tip.label,QianSeedTreePhyndr_nooutliers$tip.label[ii]))
#remove full species name and leave only genera as tip label
QianSeedTreePhyndr_nooutliers_GenusTree$tip.label<-sapply(strsplit(QianSeedTreePhyndr_nooutliers_GenusTree$tip.label,"_"),function(x) x[1])
#write output tree
write.tree(QianSeedTreePhyndr_nooutliers_GenusTree,file='./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree')

