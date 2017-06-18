#this compares big(29,703 sp) vs small (13,777 species) dataset rates 
newtable<-read.table('./output/tables/BAMM_diversification_rates_all_seed_1000.txt',header=T,sep='\t',stringsAsFactors = F)
oldtable<-read.table('./output/tables/BAMM_diversification_rates_seed.txt',header=T,sep='\t',stringsAsFactors = F)
colnames(newtable)<-paste(colnames(newtable),'_new',sep='')
colnames(newtable)[1]<-'names'
merge<-merge(oldtable,newtable,by='names',all.x=TRUE)
cor.test(merge$mean.lambda.rates,merge$mean.lambda.rates_new)
pdf('./output/plots/comparison_big_vs_small_rates.pdf')
plot(merge$mean.lambda.rates,merge$mean.lambda.rates_new,log='xy',main='rate comparison of BAMM datasets',xlab='λ rate - 13,577 species', ylab='λ rates - 29,703 species',las=1,cex.axis=.7,pch=21,bg='grey',cex=.9)
dev.off()
