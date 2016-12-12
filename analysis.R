#analysis script for seed size analysis
#recreates the analyses in the manuscript and crude versions of the plot
#warning: takes >1 day to run the full analysis

dir.create('./raw_data/')
dir.create('./raw_data/BAMM_results/')
#raw_data folder has to be included with:
#AllKewdata.txt
#QianTree.txt
#genome_size_allKew.txt
#annual_CvalueDB.txt
#perennial_CvalueDB.txt
#/BAMM_results/ with bamm mcmcout and event


#create folder structure
dir.create('./output/')
dir.create('./output/tables/')
dir.create('./output/trees/')
dir.create('./output/plots/')
dir.create('./output/plots/clade_analyses/')
dir.create('./output/BAMM_results/')

#process seed mass data from Kew Seed Information Database
####warning: standardising species names through TPL takes ~2hours
####outputs a series of tables with mean seed mass per genus (and stdev and Cov)
source('./R/KewSID_Taxonlookup.R')

#process Zanne phylogeny, overlap with seed mass dataset and get genus level dataset
####warning: standardising species names through TPL and running phyndr takes ~1 hour
####outputs genus level phylogenetic tree of genera with seed mass data
source('./R/PhylogeneticAnalysis_Phyndr.R')

#get BAMM input files and settings
source('./R/BAMM_input_generate.R')

#run BAMM with input,priors and tree file generated above
#BAMM results files have to go in /raw_data/BAMM_results

#check BAMM run convergence (burnin set to 0.2), outputs to plot folder
source('./R/analyse_BAMM_convergence.R')
pdf('./output/plots/diversification_mcmc.pdf',paper='a4r')
par(mfrow=c(2,2))
analyse_BAMM_convergence('./raw_data/BAMM_results/mcmc_out_divPhyndr25_50m1cluster.txt',0.2)
dev.off()
pdf('./output/plots/traitevolution_mcmc.pdf',paper='a4r')
par(mfrow=c(2,2))
analyse_BAMM_convergence('./raw_data/BAMM_results/mcmc_out_traitPhyndr25_2106_2.txt',0.2)
dev.off()
#plot prior vs posterior
pdf('./output/plots/trait_priorposterior_div_25.pdf',paper='a4r')
plotPriormod('./raw_data/BAMM_results/mcmc_out_divPhyndr25_50m1cluster.txt',25,0.2)
dev.off()
pdf('./output/plots/trait_priorposterior_div_25.pdf',paper='a4r')
plotPriormod('./raw_data/BAMM_results/mcmc_out_traitPhyndr25_2106_2.txt',25,0.2)
dev.off()





#analyse BAMM diversification and trait evolution
source('./R/BAMM_analysis.R')
source('./R/plots.R')
analyse_BAMM_output('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./raw_data/BAMM_results/event_data_divPhyndr25_50m1cluster.txt',traitdata='./output/tables/BAMMPhyndr_GenusSeed_trait_data.txt',burnin=0.2,mode='diversification',name='seed')
analyse_BAMM_output('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./raw_data/BAMM_results/event_data_traitPhyndr25_2106_2.txt', traitdata='./output/tables/BAMMPhyndr_GenusSeed_trait_data.txt',burnin=0.2,mode='trait',name='seed')

#write table with diversification and trait rates 
BAMMdiversification<-read.table('./output/tables/BAMM_diversification_rates_seed.txt',header=T,sep='\t')
BAMMtrait<-read.table('./output/tables/BAMM_trait_rates_seed.txt',header=T,sep='\t')
BAMMdivtrait<-merge(BAMMdiversification,BAMMtrait,by='names')
BAMMdivtrait$log.Seed.Weight.y<-NULL
colnames(BAMMdivtrait)[5]<-'log.Seed.Weight'
write.table(BAMMdivtrait,file='./output/tables/BAMMdivtrait_rates_seed.txt',sep='\t',row.names=F,quote=F)

#get lambda,mu,r and beta scaled trees (warning: takes some hours to run)
rate_scaled_trees('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./raw_data/BAMM_results/event_data_divPhyndr25_50m1cluster.txt',burnin=0.2,mode='diversification',parameter='speciation')
rate_scaled_trees('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./raw_data/BAMM_results/event_data_divPhyndr25_50m1cluster.txt',burnin=0.2,mode='diversification',parameter='extinction')
rate_scaled_trees('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./raw_data/BAMM_results/event_data_divPhyndr25_50m1cluster.txt',burnin=0.2,mode='diversification',parameter='ndr')
rate_scaled_trees('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./raw_data/BAMM_results/event_data_traitPhyndr25_2106_2.txt',burnin=0.2,mode='trait',parameter='trait')

#run STRAPP, get tables and plots
source('./R/STRAPP_plot.R')
#diversification
STRAPP_plot('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./raw_data/BAMM_results/event_data_divPhyndr25_50m1cluster.txt','./output/tables/BAMMdivtrait_rates_seed.txt',0.2,'seed')
#trait
STRAPP_plot('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./raw_data/BAMM_results/event_data_divPhyndr25_50m1cluster.txt','./output/tables/BAMMdivtrait_rates_seed.txt',0.2,'seedrate')
#generate STRAPP correlations through the posterior for plotting (as in Fig.2)
generate_correlations_strapp('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./raw_data/BAMM_results/event_data_divPhyndr25_50m1cluster.txt',traitdata = './output/tables/BAMMdivtrait_rates_seed.txt',burnin = 0.2,name='seed',parameter = 'speciation')
generate_correlations_strapp('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./raw_data/BAMM_results/event_data_divPhyndr25_50m1cluster.txt',traitdata = './output/tables/BAMMdivtrait_rates_seed.txt',burnin = 0.2,name='seed',parameter = 'extinction')
generate_correlations_strapp('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./raw_data/BAMM_results/event_data_divPhyndr25_50m1cluster.txt',traitdata = './output/tables/BAMMdivtrait_rates_seed.txt',burnin = 0.2,name='seed',parameter = 'net diversification')
generate_correlations_strapp('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./raw_data/BAMM_results/event_data_divPhyndr25_50m1cluster.txt',traitdata = './output/tables/BAMMdivtrait_rates_seed.txt',burnin = 0.2,name='seedrate',parameter = 'speciation')
generate_correlations_strapp('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./raw_data/BAMM_results/event_data_divPhyndr25_50m1cluster.txt',traitdata = './output/tables/BAMMdivtrait_rates_seed.txt',burnin = 0.2,name='seedrate',parameter = 'extinction')
generate_correlations_strapp('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./raw_data/BAMM_results/event_data_divPhyndr25_50m1cluster.txt',traitdata = './output/tables/BAMMdivtrait_rates_seed.txt',burnin = 0.2,name='seedrate',parameter = 'net diversification')

#run coefficient of variation analysis
source('./R/cov_analysis.R')
mean_cov_seed_mass('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./output/tables/BAMMPhyndr_GenusSeed_cov_trait_data.txt')
#diversification
STRAPP_plot('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./raw_data/BAMM_results/event_data_divPhyndr25_50m1cluster.txt','./output/tables/BAMMdivtrait_rates_seed.txt',0.2,'seed')
#generate STRAPP correlations through the posterior for plotting (as in Fig.2)
generate_correlations_strapp('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./raw_data/BAMM_results/event_data_divPhyndr25_50m1cluster.txt',traitdata = './output/tables/BAMMPhyndr_GenusSeed_cov_trait_data.txt',burnin = 0.2,name='seedcov',parameter = 'speciation')
generate_correlations_strapp('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./raw_data/BAMM_results/event_data_divPhyndr25_50m1cluster.txt',traitdata = './output/tables/BAMMPhyndr_GenusSeed_cov_trait_data.txt',burnin = 0.2,name='seedcov',parameter = 'extinction')
generate_correlations_strapp('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./raw_data/BAMM_results/event_data_divPhyndr25_50m1cluster.txt',traitdata = './output/tables/BAMMPhyndr_GenusSeed_cov_trait_data.txt',burnin = 0.2,name='seedcov',parameter = 'net diversification')

#seed mass, C-value and life-cycle analysis
source('./R/lifecycleCvalue_analysis.R')
#build Cvalue dataset
build_Cvalue_dataset('./raw_data/genome_size_allKew.txt')
#build life cycle dataset
build_lifecycle_dataset('./raw_data/annual_CvalueDB.txt','./raw_data/perennial_CvalueDB.txt')
#run analysis and plot strapp graphs
#Figure S5: STRAPP of seed mass, life cycle and genome size
run_seedlifecycleCvalue_STRAPP('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./raw_data/BAMM_results/event_data_divPhyndr25_50m1cluster.txt',0.2,'./output/tables/BAMMPhyndr_GenusSeed_trait_data.txt','./output/tables/KewCvalue_GenusLevel.txt','./output/tables/KewLifeCycle_GenusLevel.txt')

#run Type I error analysis
source('./R/typeI_error_strapp.R')
typeIerror_strapp('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./raw_data/BAMM_results/event_data_divPhyndr25_50m1cluster.txt',burnin=0.2,name = 'seed')

#######################curated dataset of genera analysis,outputs pgls results with epsilon = 0.5
#this runs the analysis for monophyletic genera in Zanne tree with seed size data and with minimum prob of recovering true root = 0.7
source('./R/genus_level_no_constraints.R')
#this runs the analysis for monophyletic genera with 3 or more sp with seed size data in Zanne tree and with minimum prob of recovering true root = 0.7
source('./R/genus_level_constrained.R')

#######################run clade based analyses
####this runs the clade based analyses and generates crude versions of plots (should be combined later)
#run this first to prepare the dataset
source('./R/clade_based_analyses.R')

#this runs the clade analysis for congeneric species only
for (i in seq(from=0,to=18,by=2)){
  run_clades_MSlambda_congenerics_pgls(minage=i,maxage=i+2,mincladesize=3,sampling=0.3,ncores=2,table=Kew.TaxLookUp.Qiandropped)
}
#this runs the clade analysis with Magallon & Sanderson method
for (i in seq(from=0,to=18,by=2)){
  run_clades_MSlambda_pgls(minage=i,maxage=i+2,mincladesize=3,sampling=0.3,ncores=2,table=Kew.TaxLookUp.Qiandropped)
}
#this runs the clade analysis with RPANDA method
for (i in seq(from=0,to=18,by=2)){
  run_clades_RPANDA_pgls(minage=i,maxage=i+2,mincladesize=3,sampling=0.3,ncores=2,table=Kew.TaxLookUp.Qiandropped)
}

####################get crude versions of plots (to be edited/combined later)
source('./R/plots.R')
#Figure 1: genus level phylogenetic tree with heatmap of seed mass, seed mass rate and BAMM estimated rates
plot_phylogeny_heatmap('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./output/tables/BAMMdivtrait_rates_seed.txt','maintree_heatmap')

#Figure 2 (add the strapp_density plots generated with STRAPP_plot.R): 
plot_strapp_samples_from_file('./output/tables/strapp_all_lambda_seed.txt','lambda','seed','./output/tables/BAMMdivtrait_rates_seed.txt')
plot_strapp_samples_from_file('./output/tables/strapp_all_mu_seed.txt','mu','seed','./output/tables/BAMMdivtrait_rates_seed.txt')
plot_strapp_samples_from_file('./output/tables/strapp_all_div_seed.txt','div','seed','./output/tables/BAMMdivtrait_rates_seed.txt')
plot_strapp_samples_from_file('./output/tables/strapp_all_lambda_seedrate.txt','lambda','seedrate','./output/tables/BAMMdivtrait_rates_seed.txt')
plot_strapp_samples_from_file('./output/tables/strapp_all_mu_seedrate.txt','mu','seedrate','./output/tables/BAMMdivtrait_rates_seed.txt')
plot_strapp_samples_from_file('./output/tables/strapp_all_div_seedrate.txt','div','seedrate','./output/tables/BAMMdivtrait_rates_seed.txt')

#Figure 3 (correlations across time slices)
plot_correlations(folder='./output/plots/clade_analyses/',pattern='*_pgls_MSlambda_results.txt','MSlambda_correlations',cladesizefile='./output/tables/clades_size.txt')
plot_correlations(folder='./output/plots/clade_analyses/',pattern='*_pgls_MScongenerics_results.txt','MScongenerics_correlations',cladesizefile='./output/tables/clades_congenerics_size.txt')
plot_correlations(folder='./output/plots/clade_analyses/',pattern='*_pgls_lambda_results.txt','RPANDAlambda_correlations',cladesizefile='./output/tables/clades_size.txt')

#Figure S1: posterior for number of rate shifts
plot_number_shifts_BAMM('./raw_data/BAMM_results/mcmc_out_divPhyndr25_50m1cluster.txt',0.2,'diversification')
plot_number_shifts_BAMM('./raw_data/BAMM_results/mcmc_out_traitPhyndr25_2.txt',0.2,'trait')

#Figure S2: Genus level tree with branches scaled by diversification rate and coloured by seed mass evolution rate
plot_colour_ratetree('./output/trees/ndr_transformed.tree','output/trees/trait_transformed.tree','diversification')

#Figure S3: mean genus seed mass vs coefficient of variation
plot_mean_cov_seed_mass('./output/tables/BAMMPhyndr_GenusSeed_cov_trait_data.txt')

#Figure S4: STRAPP with seed mass coefficient of variation
plot_strapp_samples_from_file('./output/tables/strapp_all_lambda_seedcov.txt','lambda','seedcov','./output/tables/BAMMdivtrait_rates_seed.txt')
plot_strapp_samples_from_file('./output/tables/strapp_all_mu_seedcov.txt','mu','seedcov','./output/tables/BAMMdivtrait_rates_seed.txt')
plot_strapp_samples_from_file('./output/tables/strapp_all_div_seedcov.txt','div','seedcov','./output/tables/BAMMdivtrait_rates_seed.txt')


#Figure S7: Family level tree with sampling information
plot_famtree_sampling('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./output/tables/Phyndr_BAMM_Genus_SamplingFractions.txt')

#Figure S8: Type I error analysis for STRAPP
plot_typeIerror('./output/tables/typeIerror_sims_1803.txt')

#Figure S9: Boxplot of seed mass in strictly annuals and perennials (+PGLS)
source('./R/lifecycleCvalue_analysis.R')
annual_vs_perennial('./output/trees/QianSeedPhyndr_nooutliers_GenusTree.tree','./output/tables/BAMMPhyndr_GenusSeed_trait_data.txt','./output/tables/KewCvalue_GenusLevel.txt','./output/tables/KewLifeCycle_GenusLevel.txt')







