library(ape)
library(BAMMtools)
#script to paste the subclades into the backbone analysis

#this function takes a subclade.tree file and its corresponding eventfile and paste the event data into the backbone event file
#age.root.Zanne is the age of the root in Zanne tree
add_subclade_to_backbone_BAMM<-function(subclade.treefile,subclade.eventfile,backbone.treefile,backbone.eventfile,newname,age.root.Zanne){
  subclade.tree<-read.tree(subclade.treefile)
  backbone.tree<-read.tree(backbone.treefile)
  clade.species.backbone<-intersect(backbone.tree$tip.label,subclade.tree$tip.label)
  #get age of clade root
  age.clade.root<-max(branching.times(subclade.tree))
  age.clade.root<-age.root.Zanne-age.clade.root
  #load subclade analysis (eventfile)
  subclade.event<-read.csv(subclade.eventfile,stringsAsFactors = F)
  #add age of clade root to absolute time
  subclade.event$abstime<-subclade.event$abstime+age.clade.root
  #load backbone analysis (eventfile)
  backbone.event<-read.csv(backbone.eventfile,stringsAsFactors = F)
  #check for shifts in branch leading to clade 
  grep.clade.species.backbone<-grep(clade.species.backbone,c(backbone.event$leftchild,backbone.event$rightchild))
  if(length(grep.clade.species.backbone)>0){
    #drop all shifts in the branch leading to the clade
    backbone.event<-backbone.event[-grep.clade.species.backbone,]
  }
  #this method of splicing things together will place a rate shift at the exact crown node of your clipped out subtrees. There is some possibility that this will cause numerical issues if the shift occurs at the exact time of an internal node
  # It would be safest to change the time of the rate shift at the  base of the subclade by an arbitrary but tiny amount, so that the shift moves a small bit down the parent branch.
  #add a small negative number (say, -0.001) to the root event for the analysis for subclade X, which will just push this shift a tiny bit closer to the root and make sure it doesn’t sit exactly on the node.
  correction<-sample(seq(from=0.0005,to=0.0015,by=0.0001),size=1)
  cat(correction,'\n')
  subclade.event[subclade.event$abstime==min(subclade.event$abstime),]$abstime<-subclade.event[subclade.event$abstime==min(subclade.event$abstime),]$abstime-correction
  merged.event<-rbind(backbone.event,subclade.event)
  merged.event<-merged.event[order(merged.event$generation),]
  #return(merged.event)
  write.csv(merged.event,file=paste('./output/BAMM_results/Zanne_clades/',newname,'.txt',sep=''),quote=F,row.names=F)
  
}
tree<-read.tree('./output/trees/BAMM_Species_tree_noseed.tree')
age.root.Zanne<-max(branching.times(tree))
#this combines the BAMM event files (the eventfiles are already burned in)
aa<-add_subclade_to_backbone_BAMM(subclade.treefile="./output/trees/Zanne_clades/clade_1.tree",subclade.eventfile="./output/BAMM_results/Zanne_clades/run2_event_data_clade_1_50.txt",backbone.treefile="./output/trees/clade_7.tree",backbone.eventfile="./output/BAMM_results/Zanne_clades/run2_event_data_clade_7_50.txt",newname = 'backbone1',age.root.Zanne)
aa<-add_subclade_to_backbone_BAMM(subclade.treefile="./output/trees/Zanne_clades/clade_2.tree",subclade.eventfile="./output/BAMM_results/Zanne_clades/run2_event_data_clade_2_50.txt",backbone.treefile="./output/trees/clade_7.tree",backbone.eventfile="./output/BAMM_results/Zanne_clades/backbone1.txt",newname = 'backbone12',age.root.Zanne)
aa<-add_subclade_to_backbone_BAMM(subclade.treefile="./output/trees/Zanne_clades/clade_3.tree",subclade.eventfile="./output/BAMM_results/Zanne_clades/run2_event_data_clade_3_50.txt",backbone.treefile="./output/trees/clade_7.tree",backbone.eventfile="./output/BAMM_results/Zanne_clades/backbone12.txt",newname = 'backbone123',age.root.Zanne)
aa<-add_subclade_to_backbone_BAMM(subclade.treefile="./output/trees/Zanne_clades/clade_4.tree",subclade.eventfile="./output/BAMM_results/Zanne_clades/run2_event_data_clade_4_50.txt",backbone.treefile="./output/trees/clade_7.tree",backbone.eventfile="./output/BAMM_results/Zanne_clades/backbone123.txt",newname = 'backbone1234',age.root.Zanne)
aa<-add_subclade_to_backbone_BAMM(subclade.treefile="./output/trees/Zanne_clades/clade_5.tree",subclade.eventfile="./output/BAMM_results/Zanne_clades/run2_event_data_clade_5_50.txt",backbone.treefile="./output/trees/clade_7.tree",backbone.eventfile="./output/BAMM_results/Zanne_clades/backbone1234.txt",newname = 'backbone12345',age.root.Zanne)
aa<-add_subclade_to_backbone_BAMM(subclade.treefile="./output/trees/Zanne_clades/clade_6.tree",subclade.eventfile="./output/BAMM_results/Zanne_clades/run2_event_data_clade_6_50.txt",backbone.treefile="./output/trees/clade_7.tree",backbone.eventfile="./output/BAMM_results/Zanne_clades/backbone12345.txt",newname = 'backbone123456',age.root.Zanne)


