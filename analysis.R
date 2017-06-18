#analysis script for seed size analysis
#recreates the analyses in the manuscript and crude versions of the plot
#warning: takes >1 day to run the full analysis

#raw_data folder has to be included with:
#AllKewdata.txt
#QianTree.txt
#genome_size_allKew.txt
#annual_CvalueDB.txt
#perennial_CvalueDB.txt
#GlobalWoodinessDatabase.csv
#height_pruned_header.txt

dir.create('./raw_data/BAMM_results/')
#create folder structure
dir.create('./output/')
dir.create('./output/tables/')
dir.create('./output/trees/')
dir.create('./output/plots/')
dir.create('./output/plots/clade_analyses/')
dir.create('./output/BAMM_results/')
dir.create('./output/BAMM_results/Zanne_clades/')
dir.create('./output/tables/Zanne_clades/')

#./output/BAMM_results/ with bamm mcmcout file and event file

###standardise datasets through TPL, and save output to tables
source('./R/clean_datasets_through_TPL.R')
#warning, takes a few hours
run_tree_through_TPL(treefile='./raw_data/Qian_tree.txt')
run_kewSID_through_TPL(traitfile='./raw_data/AllKewdata.txt')
run_lifecycle_through_TPL(annualfile='./raw_data/annual_CvalueDB.txt',perennialfile='./raw_data/perennial_CvalueDB.txt')
run_woodiness_through_TPL(traitfile='./raw_data/GlobalWoodinessDatabase.csv')

###collate Zanne tree with other phenotypic trait datasets (to run BAMM phenotypic evolution, on continuous traits) 
source('./R/prepare_phenotypic_traits_individual_analyses.R')
prepare_height_dataset()
prepare_lifecycle_Cvalue_dataset()

########Zanne phylogeny + seed size data (13,777 species)
#process Zanne phylogeny, overlap with seed mass dataset and get BAMM input files
#this is for the 13ktip (species in Zanne tree + seed data)
source('./R/BAMM_input_species_level_analysis.R')
#run BAMM with input,priors and tree file generated above
#BAMM results files have to go in /output/BAMM_results
#for speciation, we discard run 1 & run 2 and use run 3 and 4 as below
#combining runs of BAMM mcmcout
combine_multiple_BAMM_runs(c("./output/BAMM_results/run_3_mcmc_out_div_BAMM_Species_tree_50_laptop.txt","./output/BAMM_results/run_4_mcmc_out_div_BAMM_Species_tree_50_laptop.txt"),burnin=c(0,0))
combine_multiple_BAMM_runs_output(c("/Users/javier/Desktop/hydrogen/seed_size/div/run_3_mcmc_out_div_BAMM_Species_tree_100_laptop.txt","/Users/javier/Desktop/hydrogen/seed_size/div/run_4_mcmc_out_div_BAMM_Species_tree_100_laptop.txt"),burnin=c(0.2,0),name='50_div_mcmc_out.txt')
#combining runs of BAMM event data
concatenate_event_data(c("./output/BAMM_results/run_3_event_data_div_BAMM_Species_tree_50_laptop.txt","./output/BAMM_results/run_4_event_data_div_BAMM_Species_tree_50_laptop.txt"),c(0.2,0),name = '50_speciation')


#######Full Zanne phylogeny (29,703 species)
#process Zanne phylogeny, split into clades of ca 5000 species (see Fig.S6 in the paper)
source('./R/split_into_clades.R')
#BAMM results files have to go in /output/BAMM_results/Zanne_clades/
#checking convergence of Zanne clades runs
combine_multiple_BAMM_runs(c("./output/BAMM_results/Zanne_clades/mcmc_out_clade_1_50.txt","./output/BAMM_results/Zanne_clades/run1_mcmc_out_clade_1_50.txt","./output/BAMM_results/Zanne_clades/run2_mcmc_out_clade_1_50.txt"),burnin=c(1,0.99,0))
combine_multiple_BAMM_runs(c("./output/BAMM_results/Zanne_clades/mcmc_out_clade_2_50.txt","./output/BAMM_results/Zanne_clades/run1_mcmc_out_clade_2_50.txt","./output/BAMM_results/Zanne_clades/run2_mcmc_out_clade_2_50.txt"),burnin=c(1,0.99,0))
combine_multiple_BAMM_runs(c("./output/BAMM_results/Zanne_clades/mcmc_out_clade_3_50.txt","./output/BAMM_results/Zanne_clades/run1_mcmc_out_clade_3_50.txt","./output/BAMM_results/Zanne_clades/run2_mcmc_out_clade_3_50.txt"),burnin=c(1,0.99,0))
combine_multiple_BAMM_runs(c("./output/BAMM_results/Zanne_clades/mcmc_out_clade_4_50.txt","./output/BAMM_results/Zanne_clades/run1_mcmc_out_clade_4_50.txt","./output/BAMM_results/Zanne_clades/run2_mcmc_out_clade_4_50.txt"),burnin=c(1,0.99,0))
combine_multiple_BAMM_runs(c("./output/BAMM_results/Zanne_clades/mcmc_out_clade_5_50.txt","./output/BAMM_results/Zanne_clades/run1_mcmc_out_clade_5_50.txt","./output/BAMM_results/Zanne_clades/run2_mcmc_out_clade_5_50.txt"),burnin=c(1,0.99,0))
combine_multiple_BAMM_runs(c("./output/BAMM_results/Zanne_clades/mcmc_out_clade_6_50.txt","./output/BAMM_results/Zanne_clades/run1_mcmc_out_clade_6_50.txt","./output/BAMM_results/Zanne_clades/run2_mcmc_out_clade_6_50.txt"),burnin=c(1,0.99,0))
combine_multiple_BAMM_runs(c("./output/BAMM_results/Zanne_clades/mcmc_out_clade_7_50.txt","./output/BAMM_results/Zanne_clades/run1_mcmc_out_clade_7_50.txt","./output/BAMM_results/Zanne_clades/run2_mcmc_out_clade_7_50.txt"),burnin=c(1,0.99,0))
#taking the run2_files (last 50 million generations of the runs)
#this combines the BAMM event files (the eventfiles are already burned in)
source('./R/merge_BAMMevents.R')


#analyse BAMM diversification and trait evolution
source('./R/analyse_BAMM_convergence.R')
source('./R/BAMM_analysis.R')
source('./R/plots.R')

#getEventdata of speciation (no burnin because we have applied it earlier)
#save the object for later

#this is for the 13,777 species tree
tree<-read.tree('./output/trees/BAMM_Species_tree.tree')
Div_edata <- getEventData(tree, eventdata="./50_00_speciation_concatenate_event_data.txt", burnin = 0, nsamples = 1000,type='diversification')
saveRDS(Div_edata,file='./output/Div_edata_50shifts_1000samples.RDS')
#get diversification rates
analyse_BAMM_output('./output/trees/BAMM_Species_tree.tree',eventfile=Div_edata, traitdata='./output/tables/BAMM_Species_Seeddata.txt',burnin=0.0,mode='diversification',name='all_seed')

#getEventdata of trait evolution
tree<-read.tree('./output/trees/BAMM_Species_tree.tree')
Div_edata.trait <- getEventData(tree, eventdata ='/Users/javier/Desktop/hydrogen/seed_size/final_results/event_data_trait_BAMM_Species_tree_50_300mill.txt', burnin = 0.3, nsamples = 1000,type='trait')
saveRDS(Div_edata.trait,file='./output/Trait_edata_50shifts_1000samples.RDS')
#get phenotypic rates
analyse_BAMM_output('./output/trees/BAMM_Species_tree.tree',eventfile=Div_edata, traitdata='./output/tables/BAMM_Species_Seeddata.txt',burnin=0.0,mode='trait',name='all_seed')

#this is for the 29,703 species tree
tree<-read.tree('./output/trees/BAMM_Species_tree_noseed.tree')
Div_edata <- getEventData(tree, eventdata="./output/BAMM_results/Zanne_clades/backbone123456.txt", burnin = 0, nsamples = 1000,type='diversification')
saveRDS(Div_edata,file='./output/BAMM_results/Zanne_clades/all_clades_50million_1000samples.RDS')
#get diversification rates
analyse_BAMM_output('./output/trees/BAMM_Species_tree_noseed.tree',eventfile=Div_edata, traitdata='./output/tables/Species_Seeddata.txt',burnin=0.0,mode='diversification',name='all')

#write table with diversification and trait rates 
BAMMdiversification<-read.table('./output/tables/BAMM_diversification_rates_all_seed.txt',header=T,sep='\t',stringsAsFactors = F)
BAMMtrait<-read.table('./output/tables/BAMM_trait_rates_seedrate.txt',header=T,sep='\t',stringsAsFactors = F)
BAMMdivtrait<-merge(BAMMdiversification,BAMMtrait,by='names')
BAMMdivtrait$log.Seed.Weight.y<-NULL
colnames(BAMMdivtrait)[5]<-'log.Seed.Weight'
write.table(BAMMdivtrait,file='./output/tables/BAMMdivtrait_rates_all_seed.txt',sep='\t',row.names=F,quote=F)

#get height phenotypic rates
Height_edata<-readRDS('./output/BAMM_results/Height_edata_trait.RDS')
analyse_BAMM_output(treefile = './output/trees/BAMM_Species_tree_height.tree',Height_edata,traitdata='./output/tables/BAMM_Species_tree_height_table.txt',burnin=0.1,mode='trait',name='heightrate')

#get cvalue phenotypic rates
Cvalue_edata<-readRDS('./output/BAMM_results/Cvalue_edata_trait.RDS')
analyse_BAMM_output(treefile = './output/trees/BAMM_Species_tree_Cvalue.tree',Cvalue_edata,traitdata='./output/tables/BAMM_Species_tree_Cvalue_table.txt',burnin=0.1,mode='trait',name='Cvaluerate')


#get lambda,mu,r and beta scaled trees (warning: takes some hours to run)
Div_edata<-readRDS('./output/BAMM_results/Zanne_clades/all_clades_50million_1000samples.RDS')
Trait_edata<-readRDS('./output/Trait_edata_50shifts_1000samples.RDS')
rate_scaled_trees('./output/trees/BAMM_Species_tree.tree',eventfile = Div_edata,burnin=0,mode='diversification',parameter='speciation')
rate_scaled_trees('./output/trees/BAMM_Species_tree.tree',eventfile = Div_edata,burnin=0,mode='diversification',parameter='extinction')
rate_scaled_trees('./output/trees/BAMM_Species_tree.tree',eventfile = Div_edata,burnin=0,mode='diversification',parameter='ndr')
rate_scaled_trees('./output/trees/BAMM_Species_tree.tree',eventfile = Trait_edata,burnin=0.3,mode='trait',parameter='trait')

#this compares rates from the big dataset (29,703 species) vs the small dataset (13,777 species) and tests correlation
source('./R/compare_big_vs_small_rates.R')

####################SORT THIS OUT############################16/06/17
#######WHERE TO GET THE RATES OF DIVERSIFICATION OF HEIGHT, GENOME SIZE

#run STRAPP, get tables and plots
source('./R/STRAPP_plot.R')
source('./R/plots.R')
#diversification
Div_edata<-readRDS('./output/BAMM_results/Zanne_clades/all_clades_50million_1000samples.RDS')
#subset the event data
traitdata<-read.table('./output/tables/BAMMdivtrait_rates_all_seed.txt',header=T,sep='\t',stringsAsFactors = F)

#have to subsample eventdata object with tree here
Div_edata.sub<-subtreeBAMM(Div_edata,tips=traitdata$names)
saveRDS(Div_edata.sub,file='./output/BAMM_results/Zanne_clades/Div_edata_50shifts_1000samples_subsampled.RDS')
#diversification
STRAPP_plot('./output/trees/BAMM_Species_tree_noseed.tree',eventfile=Div_edata.sub,'./output/tables/BAMMdivtrait_rates_all_seed_1000.txt',burnin=0,'seed')
#trait
STRAPP_plot('./output/trees/BAMM_Species_tree_noseed.tree',eventfile=Div_edata.sub,'./output/tables/BAMMdivtrait_rates_all_seed_1000.txt',burnin=0,'seedrate')
#generate STRAPP correlations through the posterior for plotting (as in Fig.2)
generate_correlations_strapp('./output/trees/BAMM_Species_tree.tree',Div_edata.sub,traitdata = './output/tables/BAMMdivtrait_rates_seed.txt',burnin = 0,name='seed',parameter = 'speciation')
generate_correlations_strapp('./output/trees/BAMM_Species_tree.tree',Div_edata.sub,traitdata = './output/tables/BAMMdivtrait_rates_seed.txt',burnin = 0,name='seed',parameter = 'extinction')
generate_correlations_strapp('./output/trees/BAMM_Species_tree.tree',Div_edata.sub,traitdata = './output/tables/BAMMdivtrait_rates_seed.txt',burnin = 0,name='seed',parameter = 'net diversification')
generate_correlations_strapp('./output/trees/BAMM_Species_tree.tree',Div_edata.sub,traitdata = './output/tables/BAMMdivtrait_rates_seed.txt',burnin = 0,name='seedrate',parameter = 'speciation')
generate_correlations_strapp('./output/trees/BAMM_Species_tree.tree',Div_edata.sub,traitdata = './output/tables/BAMMdivtrait_rates_seed.txt',burnin = 0,name='seedrate',parameter = 'extinction')
generate_correlations_strapp('./output/trees/BAMM_Species_tree.tree',Div_edata.sub,traitdata = './output/tables/BAMMdivtrait_rates_seed.txt',burnin = 0,name='seedrate',parameter = 'net diversification')

#seed mass, C-value, lifecycle and woodiness analysis
source('./R/lifecycle_cvalue_habit_height_analysis.R')
##this creates the dataset
prepare_lifecycle_Cvalue_woodiness_dataset()
prepare_lifecycle_Cvalue_woodiness_height_dataset()
##this compares means of seed weight in annual vs perennials using phylogenetic anova
lifecycle_seed_phylanova(treefile = './output/trees/BAMM_Species_tree.tree',traitfile='./output/tables/seedlifecycleCvaluewoodinessheight_phenorates.txt')

#run analysis and plot strapp graphs
#Figure S5: STRAPP of seed mass, life cycle and genome size
Div_edata<-readRDS('./output/BAMM_results/Zanne_clades/Div_edata_50shifts_1000samples_subsampled.RDS')
run_seedlifecycleCvaluewoodinessheight_phenorate_STRAPP(treefile='./output/trees/BAMM_Species_tree.tree',eventfile = Div_edata,burnin=0,traitfile='./output/tables/kewsid_lifecycle_woodiness_height_TPL.txt',seedratetable = './output/tables/BAMMdivtrait_rates_all_seed_1000.txt',Cvalueratetable = './output/tables/BAMM_trait_rates_Cvaluerate.txt',heightratetable = './output/tables/BAMM_trait_rates_heightrate.txt')
#run Type I error analysis
source('./R/typeI_error_strapp.R')
Div_edata<-readRDS('./output/Div_edata_50shifts_1000samples.RDS')
typeIerror_strapp('./output/trees/BAMM_Species_tree_noseed.tree',eventfile = Div_edata,burnin=0,name = 'seed')


#######################curated dataset of genera analysis,outputs pgls results with epsilon = 0.5
#this runs the analysis for monophyletic genera in Zanne tree with seed size data and with minimum prob of recovering true root = 0.7
#source('./R/genus_level_no_constraints.R')
#this runs the analysis for monophyletic genera with 3 or more sp with seed size data in Zanne tree and with minimum prob of recovering true root = 0.7
#source('./R/genus_level_constrained.R')

#######################run clade based analyses
####this runs the clade based analyses and generates crude versions of plots (should be combined later)
#run this first to prepare the dataset
source('./R/clade_based_analyses_noseed.R')
#get the phylogroups with get_phylogroups (see in clade_based_analyses_noseed.R)
#this runs the clade analyses (with the method-of-moments estimator) across the saved phylogroups from 0 to 20 myr
for (i in seq(from=0,to=18,by=2)){
  run_clades_MS_EB_z0_pgls_savedphylogroups(tree=Qian_dropped,minage=i,maxage=i+2,mincladesize=3,sampling=0.3,ncores=2,table=Qian_Lookup_Angios_TPL_seed)
}
#this runs the clade analyses with congenerics only (with the method-of-moments estimator) across the saved phylogroups from 0 to 20 myr
for (i in seq(from=0,to=18,by=2)){
  run_clades_MSlambda_z0_congenerics_pgls_savedphylogroups(tree=Qian_dropped,minage=i,maxage=i+2,mincladesize=3,sampling=0.3,ncores=2,table=Qian_Lookup_Angios_TPL_seed)
}
#this runs the clade analyses (with RPANDA) across the saved phylogroups from 0 to 20 myr
for (i in seq(from=0,to=18,by=2)){
  run_clades_RPANDA_EB_z0_pgls_savedphylogroups_6models(tree=Qian_dropped,minage=i,maxage=i+2,mincladesize=3,sampling=0.3,ncores=2,table=Qian_Lookup_Angios_TPL_seed)
}

#############################get crude versions of plots (to be edited/combined later)
source('./R/plots.R')
#Figure 1: genus level phylogenetic tree with heatmap of seed mass, seed mass rate and BAMM estimated rates
plot_phylogeny_heatmap('./output/trees/BAMM_Species_tree.tree','./output/tables/BAMMdivtrait_rates_all_seed_1000.txt','maintree_heatmap_all_seed_1000')

#Figure 2 (add the strapp_density plots generated with STRAPP_plot.R): 
plot_strapp_samples_from_file('./output/tables/strapp_all_lambda_seed_1000.txt','lambda','seed','./output/tables/BAMMdivtrait_rates_all_seed_1000.txt')
plot_strapp_samples_from_file('./output/tables/strapp_all_mu_seed_1000.txt','mu','seed','./output/tables/BAMMdivtrait_rates_all_seed_1000.txt')
plot_strapp_samples_from_file('./output/tables/strapp_all_div_seed_1000.txt','div','seed','./output/tables/BAMMdivtrait_rates_all_seed_1000.txt')
plot_strapp_samples_from_file('./output/tables/strapp_all_lambda_seedrate_1000.txt','lambda','seedrate','./output/tables/BAMMdivtrait_rates_all_seed_1000.txt')
plot_strapp_samples_from_file('./output/tables/strapp_all_mu_seedrate_1000.txt.txt','mu','seedrate','./output/tables/BAMMdivtrait_rates_all_seed_1000.txt')
plot_strapp_samples_from_file('./output/tables/strapp_all_div_seedrate_1000.txt.txt','div','seedrate','./output/tables/BAMMdivtrait_rates_all_seed_1000.txt')

##include a clades_size.txt in ./output/tables/
#like this:
#Slice	Nclades
#0_2	131
#2_4	161
#4_6	144
#6_8	148
#8_10	136
#10_12	121
#12_14	95
#14_16	77
#16_18	70
#18_20	69
####
#Figure 3 (correlations RPANDA lambda across time slices)
plot_correlations(folder='./output/clade_analyses/RPANDA_EB_z0/',pattern='*_pgls_RPANDAlambda_results.txt',name = 'RPANDAlambda_correlations',cladesizefile='./output/tables/clades_size.txt')

###########SUPPLEMENTARY INFO

#Figure S1: Species tree with branches scaled by diversification rate and coloured by seed mass evolution rate
plot_colour_ratetree('./output/trees/speciation_transformed.tree','output/trees/beta_transformed.tree','diversification')

#plot prior vs posterior for Fig S2 for speciation/extinction analysis
pdf('./output/plots/div_priorposterior_div_25.pdf',paper='a4r')
plotPriormod(mcmc = './output/BAMM_results/25_div_mcmc_out.txt',expectedNumberOfShifts = 25,burnin = 0,maxNumberShifts = 300)
dev.off()
pdf('./output/plots/div_priorposterior_div_50.pdf',paper='a4r')
plotPriormod(mcmc = './output/BAMM_results/50_div_mcmc_out.txt',expectedNumberOfShifts = 50,burnin = 0,maxNumberShifts = 300)
dev.off()
pdf('./output/plots/div_priorposterior_div_100.pdf',paper='a4r')
plotPriormod(mcmc = './output/BAMM_results/div_100_mcmc_out.txt',expectedNumberOfShifts = 100,burnin = 0,maxNumberShifts = 300)
dev.off()
pdf('./output/plots/div_priorposterior_div_250.pdf',paper='a4r')
plotPriormod(mcmc = './output/BAMM_results/250_div_mcmc_out.txt',expectedNumberOfShifts = 250,burnin = 0,maxNumberShifts = 300)
dev.off()
#plot prior vs posterior for Fig S2 for trait analysis
pdf('./output/plots/trait_priorposterior_trait_25.pdf',paper='a4r')
plotPriormod(mcmc = './output/BAMM_results/mcmc_out_trait_BAMM_Species_tree_25_300mill.txt',expectedNumberOfShifts = 25,burnin = 0.3,maxNumberShifts = 800)
dev.off()
pdf('./output/plots/trait_priorposterior_trait_50.pdf',paper='a4r')
plotPriormod(mcmc = './output/BAMM_results/mcmc_out_trait_BAMM_Species_tree_50_300mill.txt',expectedNumberOfShifts = 50,burnin = 0.2,maxNumberShifts = 800)
dev.off()
pdf('./output/plots/trait_priorposterior_trait_100.pdf',paper='a4r')
plotPriormod(mcmc = './output/BAMM_results/mcmc_out_trait_BAMM_Species_tree_100_300mill.txt',expectedNumberOfShifts = 100,burnin = 0.2,maxNumberShifts = 800)
dev.off()
pdf('./output/plots/trait_priorposterior_trait_250.pdf',paper='a4r')
plotPriormod(mcmc = './output/BAMM_results/mcmc_out_trait_BAMM_Species_tree_250_300mill.txt',expectedNumberOfShifts = 250,burnin = 0.2,maxNumberShifts = 800)
dev.off()
pdf('./output/plots/trait_priorposterior_trait_500.pdf',paper='a4r')
plotPriormod(mcmc = './output/BAMM_results/mcmc_out_trait_BAMM_Species_tree_500_300mill.txt',expectedNumberOfShifts = 500,burnin = 0.3,maxNumberShifts = 800)
dev.off()

##Figure S3 (correlations MS lambda across time slices)
plot_correlations(folder='./output/clade_analyses/MS_EB_z0/',pattern='*_pgls_lambda_results.txt',name = 'MSlambda_correlations',cladesizefile='./output/tables/clades_size.txt')

#Figure S6: Family level tree with sampling information
plot_famtree_sampling('./output/trees/BAMM_Species_tree_noseed.tree','./output/tables/BAMM_Species_FamilySamplingFractions_noseed.txt')

##Figure S9 (correlations MS lambda across time slices congenerics only)
##include a clades_size.txt in ./output/tables/
#like this:
#Slice	Nclades
#0_2	68
#2_4	85
#4_6	70
#6_8	62
#8_10	38
#10_12	38
#12_14	23
#14_16	20
#16_18	17
#18_20	19
####
plot_correlations_congenerics(folder='./output/plots/',pattern='*_pgls_MSlambda_congenerics_results.txt',name = 'MSlambdacongenerics_correlations',cladesizefile='./output/tables/clades_congenerics_size.txt')

#Figure S10 Correlation between speciation rate and rate of seed size evolution in a random sample of the BAMM posterior
###plot correlation between speciation and rate of seed size for a random sample of the posterior 
Div_edata<-readRDS('./output/BAMM_results/Zanne_clades/Div_edata_50shifts_1000samples_subsampled.RDS')
plot_randomposterior_correlation(Div_edata,tablefile='./output/tables/BAMMdivtrait_rates_all_seed_1000.txt')

#Figure S8: Type I error analysis for STRAPP
plot_typeIerror('./output/tables/typeIerror_sims_seed.txt')

#Figure S11: Boxplot of seed mass in strictly annuals and perennials (+PGLS)
source('./R/lifecycle_cvalue_habit_analysis.R')
annual_vs_perennial(traitfile = './output/tables/kewsid_lifecycle_woodiness_TPL.txt')

#Figure S12:Correlation of mean speciation and mean phenotypic rates across all the branches of the angiosperm phylogenetic tree
plot_OLS_speciation_trait('./output/trees/speciation_transformed.tree','./output/trees/beta_transformed.tree')

#crude plots for Figures S13,S14,S17,S18 are generated above (in clade_based_analyses_no_seed)

#Figure S15  (correlations MS NDR across time slices congenerics only)
plot_correlations(folder='./output/clade_analyses/RPANDA_EB_z0/',pattern='*_pgls_ndr_results.txt',name = 'RPANDAndr_correlations',cladesizefile='./output/tables/clades_size.txt')




