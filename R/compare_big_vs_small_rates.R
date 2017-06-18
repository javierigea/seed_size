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
