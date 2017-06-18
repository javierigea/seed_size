#simulate traits from Rabosky Huang 2015
#treefile: tree to run type I error test
#name: name to append to outputfiles
library(ape)
library(BAMMtools)
typeIerror_strapp<-function(treefile,eventfile,burnin,name){
  #function to simulate traits matrix (Rabosky and Huang 2015)
  simulateTraitsMatrix <- function(phy, sigma2 = 1, rootstate=0, reps=100, verbose=T){  
    states <- matrix(0, nrow=nrow(phy$edge), ncol=reps);  
    root <- length(phy$tip.label) + 1;
    isRootEdge <- phy$edge[,1] == root;
    rootstate <- rep(rootstate, reps);  
    states[isRootEdge, ][1,] <- rootstate + rnorm(reps, sd=sqrt(phy$edge.length[isRootEdge][1] * sigma2));
    states[isRootEdge, ][2,] <- rootstate + rnorm(reps, sd=sqrt(phy$edge.length[isRootEdge][2] * sigma2));  
    for (i in 1:nrow(phy$edge)){
      if (verbose){
        cat(i, '\n');
      }
      node <- phy$edge[i,2];
      par <- phy$edge[i,1];
      if (par != root){
        par.state <- states[phy$edge[,2] == par, ];
        states[i,] <- par.state + rnorm(reps, sd=sqrt(phy$edge.length[i] * sigma2));
      }
      
    }
    rownames(states) <- phy$edge[,2];
    phy$states <- states[as.character(1:length(phy$tip.label)), ];  
    return(phy)
  }
  #############################
  
  tree<-read.tree(treefile)
  if (class(eventfile)!='bammdata'){
    Div_edata <- getEventData(tree, eventfile, burnin = burnin, nsamples = 1000,type='diversification')
  }else if (class(eventfile)=='bammdata'){
    Div_edata<-eventfile
  }
  
  tree<-drop.tip(tree,setdiff(tree$tip.label,Div_edata$tip.label))
  typeItree<-simulateTraitsMatrix(tree,reps = 1000)
  #perform permutation on each of the simulated trait data
  p1<-sapply(1:dim(typeItree$states)[2],function(i){
    cat(i, '\n')
    lt<-typeItree$states[,i]
    names(lt)<-typeItree$tip.label
    #one tailed test
    l<-traitDependentBAMM(Div_edata,lt,1000,return.full = F,logrates = T,two.tailed = T)
    l$p.value
  })
  #calculate the type I error rate.
  error<-sum(p1<=0.05)/length(p1)
  simfile<-paste('./output/tables/typeIerror_sims_',name,'.txt',sep='')
  write(as.character(p1),simfile,sep='')
  #h<-hist(p1,breaks=20)
  #h$density = h$counts/sum(h$counts)
  #plot(h,freq=FALSE,ylim=c(0,1),xlim=c(0,1),yaxs='i',xaxs='i',las=1,ylab='Frequency',xlab='p-value')
  
}
