#Seed Genus Level data-prepare trait data for BAMM****MEAN
#read Seed Table
Kew.Genera.Seed<-read.table('./output/tables/KewSID_GenusLevelSeeddata_taxonlookup.txt',sep='\t',header=T)
#read phylogenetic tree
QianSeedTreePhyndr_nooutliers_GenusTree<-read.tree('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree')
treeQianPhyndrdf<-data.frame(QianSeedTreePhyndr_nooutliers_GenusTree$tip.label)
colnames(treeQianPhyndrdf)<-'tip.label'
BAMM.Phyndr.Trait<-merge(treeQianPhyndrdf,Kew.Genera.Seed,by.x='tip.label',by.y='genus')
BAMM.Phyndr.Trait$log.Seed.Weight<-log10(BAMM.Phyndr.Trait$Seed.Weight)
BAMM.Phyndr.Trait$Seed.Weight<-NULL
write.table(BAMM.Phyndr.Trait, file = './output/tables/BAMMPhyndr_GenusSeed_trait_data.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names=FALSE)

#Seed Genus Level data-prepare trait data for BAMM ***COEFFICIENT OF VARIATION
#read Seed Table
Kew.Genera.Seed.cov<-read.table('./output/tables/KewSID_GenusLevelSeeddata_Cov_taxonlookup.txt',sep='\t',header=T)
#read phylogenetic tree
QianSeedTreePhyndr_nooutliers_GenusTree<-read.tree('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree')
treeQianPhyndrdf<-data.frame(QianSeedTreePhyndr_nooutliers_GenusTree$tip.label)
colnames(treeQianPhyndrdf)<-'tip.label'
BAMM.Phyndr.Trait.cov<-merge(treeQianPhyndrdf,Kew.Genera.Seed.cov,by.x='tip.label',by.y='genus')
write.table(BAMM.Phyndr.Trait.cov, file = './output/tables/BAMMPhyndr_GenusSeed_cov_trait_data.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names=T)


#prepare incomplete sampling file for BAMM
#get total info from TaxonLookUp
plant_lookup_version_current()
PlantLookup<-plant_lookup(include_counts = TRUE)
#get total number of genera in each family
PlantLookup.2<-aggregate(genus~family, data=PlantLookup, length)
PlantLookup<-merge(PlantLookup, PlantLookup.2, by.x = 'family', by.y = 'family', all.x = TRUE)
colnames(PlantLookup)[3]<-'Genus'
colnames(PlantLookup)[6]<-'number.of.genera'
#run TaxonLookup on taxa from the tree
QianSeedTreePhyndr_nooutliers_GenusTree.table<-lookup_table(QianSeedTreePhyndr_nooutliers_GenusTree$tip.label)
QianPhyndrGeneracount<-aggregate(genus~family, data=QianSeedTreePhyndr_nooutliers_GenusTree.table, length)
QianPhyndrSeednooutilers_GenusTree.table<-merge(QianSeedTreePhyndr_nooutliers_GenusTree.table,QianPhyndrGeneracount,by.x='family',by.y='family')
colnames(QianPhyndrSeednooutilers_GenusTree.table)[2]<-'Genus'
colnames(QianPhyndrSeednooutilers_GenusTree.table)[5]<-'number.of.genera'
#merge tree info with full taxonlookup info
Phyndr_Genus_Incomplete<-merge(QianPhyndrSeednooutilers_GenusTree.table, PlantLookup, by.x = 'Genus', by.y = 'Genus')
Phyndr_Genus_Incomplete$family.y<-NULL
Phyndr_Genus_Incomplete$order.y<-NULL
Phyndr_Genus_Incomplete$group.y<-NULL
colnames(Phyndr_Genus_Incomplete)[c(2,3,4,5,6,7)]<-c('family','order','group','n.genera.tree','n.species','n.total.genera')
#calculate sampling fractions
Phyndr_Genus_Incomplete$Sampling_Fraction<-Phyndr_Genus_Incomplete$n.genera.tree/Phyndr_Genus_Incomplete$n.total.genera
Phyndr_BAMM_Sampling<-Phyndr_Genus_Incomplete[,c("Genus","family","Sampling_Fraction")]
#write output files
write.table(Phyndr_BAMM_Sampling, file = './output/tables/Phyndr_BAMM_Genus_SamplingFractions.txt', row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
Phyndr_backbone<-nrow(Phyndr_Genus_Incomplete)/nrow(PlantLookup)
write.table(Phyndr_backbone,file = './output/tables/Phyndr_BAMM_backbone_sampling.txt',row.names=FALSE,col.names = FALSE,quote=FALSE)
#get Family Sampling info for Supplementary Figure
Phyndr_Family_Sampling<-unique(data.frame(Phyndr_BAMM_Sampling$family,Phyndr_BAMM_Sampling$Sampling_Fraction))
write.table(Phyndr_Family_Sampling, file = './output/tables/Phyndr_BAMM_Family_Sampling_Fractions.txt', quote = FALSE, sep = '\t', row.names = FALSE)

#checks for BAMM
is.ultrametric(QianSeedTreePhyndr_nooutliers_GenusTree)
min(QianSeedTreePhyndr_nooutliers_GenusTree$edge.length) 
summary(duplicated(QianSeedTreePhyndr_nooutliers_GenusTree$tip.label))

###BAMM priors
setBAMMpriors(QianSeedTreePhyndr_nooutliers_GenusTree,outfile='./output/tables/BAMM_diversification_Priors.txt')
setBAMMpriors(QianSeedTreePhyndr_nooutliers_GenusTree,traits=,outfile='./output/tables/BAMM_trait_Priors.txt')

