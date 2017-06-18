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

#upload Kew data from 'All Kew Data - unformatted2'####J: changed name here to AllKewdata.txt
Kew<-read.table('./raw_data/AllKewdata.txt', sep="\t", header = T)
#remove duplicates
Kew<-unique(Kew)
#run through TPL and add to original table
Kew.corrected<-TPL(genus = Kew$Genus, species = Kew$Species)

write.table(Kew.corrected, file = './output/tables/AllKewdata_duplicated_removed_TPL.txt', quote = FALSE, row.names = FALSE, sep = '\t')
#add data to Kew table
Kew[,'Confirmed.Family']<-paste(Kew.corrected$Family)
Kew[,'Confirmed.Genus']<-paste(Kew.corrected$New.Genus)
Kew[,'Confirmed.Species']<-paste(Kew.corrected$New.Species)
Kew[,'Merged']<-paste(Kew.corrected$New.Genus, Kew.corrected$New.Species, sep = '_')
#remove hybrids 'x' and species = 'sp' etc
Kew<-subset(Kew,Kew$Confirmed.Species !='sp')
Kew<-subset(Kew,Kew$Confirmed.Species !='sp.')
Kew<-subset(Kew,Kew$Confirmed.Species !='spp.')
Kew<-subset(Kew,Kew$Confirmed.Species !='NA')
Kew<-subset(Kew,Kew$Confirmed.Species !='af')
Kew<-subset(Kew,Kew$Confirmed.Species !='c')
Kew<-subset(Kew,Kew$Confirmed.Species !='sect.')
Kew<-subset(Kew,Kew$Confirmed.Species !='x')
#remove blanks and 00.0000e0g from seed weight col - these are recorded as 0g on Kew Database
Kew<-subset(Kew, Kew$Merged !='Aneilema_setiferum')
Kew<-subset(Kew, Kew$Merged !='Heimia_salicifolia')
Kew<-subset(Kew, Kew$Merged !='Pluchea_odorata')
Kew<-subset(Kew, Kew$Merged !='Psychotria_vellosiana')
Kew<-subset(Kew, Kew$SeedWeight_g !='')
Kew<-subset(Kew, Kew$SeedWeight_g !='g.')
Kew<-subset(Kew, Kew$SeedWeight_g !='Goetgh.,')
Kew<-subset(Kew, Kew$SeedWeight_g !='G.o.')
#get one value of seed weight for each species (the mean)
Kew[, 'Seed.Weight']<-sub(Kew$SeedWeight_g, pattern = 'g', replacement = '')
Kew$Seed.Weight<-as.numeric(Kew$Seed.Weight)
uniqueKew<-unique(Kew)
uniqueKew$Mergedfactor<-as.factor(uniqueKew$Merged)
Kewuniqueseed<-aggregate(uniqueKew$Seed.Weight~uniqueKew$Mergedfactor, FUN=mean)

#length(unique(Kewuniqueseed$`uniqueKew$Mergedfactor`))
colnames(Kewuniqueseed)[1] <- 'Mergedfactor'
colnames(Kewuniqueseed)[2] <- 'Seed.Weight'
uniqueKewdf<-data.frame(uniqueKew$Mergedfactor,uniqueKew$Confirmed.Family,uniqueKew$Confirmed.Genus,uniqueKew$Confirmed.Species)
uniqueKewdf<-unique(uniqueKewdf)
uniqueKewdf<-uniqueKewdf[!duplicated(uniqueKewdf$uniqueKew.Mergedfactor),]
Kew2<-merge(Kewuniqueseed,uniqueKewdf, by.x = 'Mergedfactor', by.y = 'uniqueKew.Mergedfactor')
colnames(Kew2)[1]<-'Name'
colnames(Kew2)[3]<-'Confirmed.Family'
colnames(Kew2)[4]<-'Confirmed.Genus'
colnames(Kew2)[5]<-'Confirmed.Species'
#Kew SID database lists SeedWeight as SeedWeight of 1000 seeds
Kew2$Seed.Weight<-Kew2$Seed.Weight/1000
#run through taxonlookup
Kew2$Name<-as.character(Kew2$Name)
Kew_Lookup<-lookup_table(Kew2$Name, by_species = TRUE)
####################
Kew_Lookup$Species<-row.names(Kew_Lookup)
row.names(Kew_Lookup)<-NULL
Kew.TaxLookUp<-merge(Kew_Lookup,Kew2, by.x = 'Species', by.y = 'Name')
Kew.TaxLookUp$Confirmed.Family<-NULL
Kew.TaxLookUp$Confirmed.Genus<-NULL
Kew.TaxLookUp$Confirmed.Species<-NULL
Kew.TaxLookUp$Seed.Weight<-as.numeric(Kew.TaxLookUp$Seed.Weight)
Kew.TaxLookUp$Log10.Seed.Weight<-log10(Kew.TaxLookUp$Seed.Weight)

write.table(Kew.TaxLookUp, file = './output/tables/KewSIDdata_taxonlookup.txt', quote = FALSE, sep = '\t', row.names = FALSE)

#get Genus-level Seed Data****MEAN
Kew.Genera.Seed<-aggregate(Seed.Weight~genus, data=Kew.TaxLookUp, mean)
write.table(Kew.Genera.Seed, file = './output/tables/KewSID_GenusLevelSeeddata_taxonlookup.txt', quote = FALSE, sep = '\t', row.names = FALSE)
#get Genus-level Seed Data****STDEV
Kew.Genera.Seed.sd<-aggregate(Seed.Weight~genus, data=Kew.TaxLookUp, sd)
Kew.Genera.Seed.sd$log10<-log10(Kew.Genera.Seed.sd$Seed.Weight)
colnames(Kew.Genera.Seed.sd)<-c('genus','sd.Seed.Weight','log10.Seed.Weight')
write.table(Kew.Genera.Seed.sd, file = './output/tables/KewSID_GenusLevelSeeddata_SD_taxonlookup.txt', quote = FALSE, sep = '\t', row.names = FALSE)
#get Genus-level Seed Data****Coefficient of Variation
Kew.Genera.Seed.mean<-aggregate(Seed.Weight~genus, data=Kew.TaxLookUp, mean)
Kew.Genera.Seed.stdev<-aggregate(Seed.Weight~genus, data=Kew.TaxLookUp, sd)
Kew.Genera.Seed.CoV<-merge(Kew.Genera.Seed.mean,Kew.Genera.Seed.stdev,by.x='genus',by.y='genus')
colnames(Kew.Genera.Seed.CoV)[2]<-'mean'
colnames(Kew.Genera.Seed.CoV)[3]<-'sd'
Kew.Genera.Seed.CoV$cov<-(Kew.Genera.Seed.CoV$sd/Kew.Genera.Seed.CoV$mean)*100
write.table(Kew.Genera.Seed.CoV, file = './output/tables/KewSID_GenusLevelSeeddata_Cov_taxonlookup.txt', quote = FALSE, sep = '\t', row.names = FALSE)


