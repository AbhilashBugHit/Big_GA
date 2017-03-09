Lookup<-function(set1,synlethNA){
  # take two genes and look them up in the epistatic interaction
  # matrix and retrieve their scores
  r1<-match(set1,rownames(synlethNA))
  c1<-match(set1,colnames(synlethNA))
  subsyn<-synlethNA[r1,c1]
  rownames(subsyn)<-NULL
  colnames(subsyn)<-NULL
  epi_score<-unlist(subsyn)
  epi_score<-epi_score[which(!is.nan(epi_score))]
  #Edit 18-1-2016 removing NA's while Lookup because downstream Fitness_fx is giving NA's
  epi_score<-epi_score[which(!is.na(epi_score))]
  return(epi_score)
}

Init_Pop<-function(members=vector(),chrom_size=8,numbers=20){
  # Generate Initial Population
  populus<-vector(mode="list",length=numbers)
  for(i in 1:length(populus))
  {
    populus[[i]]<-sample(members,size=chrom_size,replace=TRUE) 
  }
  return(populus)
}

Fitness_fx<-function(list_element,synlethNA){
  
  # Calculating SNR combinatorial fitness scores for a set of genes
  #create all sets of 2 for a set of genes
  #
  all2<-combn(list_element,2)
  Median_PWEScores<-apply(all2,2,FUN=function(x){return(median(Lookup(x,synlethNA)))})
  SNR<-mean(Median_PWEScores,na.rm=TRUE)/var(Median_PWEScores,na.rm=TRUE)
  
  # As long as you are not removing the negative epistasis scores, you will be fine
  # with estimating on an average (in a distributional sense)
  # whether the deletion is alleviating or aggravating.
  
  # Alternate fitness function that attempts to minimize co-efficient of variation
  # and also prevent errors from shrinking to near zero by apply negative log
  # transform, this should hopefully prevent the variance from exploding as well.
  #SNR<- (-1)* sd(Median_PWEScores,na.rm=TRUE)/mean(Median_PWEScores,na.rm=TRUE)
  return(SNR)
}

Fitness_score<-function(populus,synlethNA){
  # Calculate Fitness Score of population using lookup from synlethNA
  names(populus)<-unlist(lapply(X=populus,FUN=function(x){return(Fitness_fx(x,synlethNA))}))
  return(populus)
}

Selection<-function(populus,PopCount=NULL){
  # Roulette wheel selection
  # select any one solution randomly, generate a random number, if the random number is lesser
  # than the cumulative sum value of fitness from the population as you start adding up 
  
  print("Selection")
  if(length(PopCount)==0){PopCount=length(populus)}
  scores<-as.numeric(names(populus))
  if(length(scores)==0){stop("Population Not Scored Yet! Use Fitness_score to do so.")}
  #! when you do THIS you disregard negative pairs of epistasis scores.
  scores[which(scores<0)]=0
  summatScore<-sum(scores,na.rm = TRUE)
  RandNo<-runif(PopCount,0,summatScore) # Roulette Spin
  RunTot<-cumsum(scores) # Cumulative Spin for area
  randGtIx<-sapply(RandNo,function(x){which(RunTot>x)})
  # Clean up the entries where none of the fitness values were greater than the random number
  randGtIx<-randGtIx[which(unlist(lapply(randGtIx,function(x){length(x)}))!=0)]
  # unlist the list of indices where the population member was greater than random number
  # and select only the first index which is the first instance where the cumulative sum
  # became greater than the random number that was generated
  FirstIx<-unlist(lapply(randGtIx,function(x){x[[1]]}))
  newpop<-populus[FirstIx]
  # reselect the members of the population with the set of first indices
  return(newpop)
}

Splice<-function(x,synlethNA){
  #select any location between the 2nd and the 2nd last to cut and cross over
  spliceSite<-sample(c(2:(length(x[[1]])-1)),1)
  xtemp<-c(x[[1]][1:spliceSite],x[[2]][(spliceSite[1]+1):length(x[[2]])])
  x[[2]]<-c(x[[2]][1:spliceSite],x[[1]][(spliceSite[1]+1):length(x[[1]])])
  x[[1]]<-xtemp
  x<-Fitness_score(x,synlethNA)
  return(x)
}

SpliceR<-function(populus,synlethNA,proportion=0.5){
  # SpliceR - you give it a population it randomly pairs them and splices the pairs
  # sends back a rescored population for better or for worse (Scores)
  print("Splicing")
  randix<-matrix(sample(c(1:length(populus))),,2)
  PopProp<-round(proportion*dim(randix)[1]) # Selecting a proportion of the population
  randix<-randix[1:PopProp,] # subsetting the matrix to reflect the proportion
  RandSplice<-apply(randix,1,function(x){Splice(populus[x],synlethNA)})
  Dotterpop<-list()
  for(i in 1:length(RandSplice))
  {
    Dotterpop<-c(Dotterpop,RandSplice[[i]])
  } # there was a problem with unlisting which is why I used a loop to unlist manually
  return(Dotterpop)
}

Immig<-function(populus,synlethNA,proportion=0.01){
  
  #Immigration - concatenate/add some randomly selected new members to the population assuming that
  # previously in the random selection performed during population intialization, they
  # didn't make it to the list and they might be fit solutions and you never know until
  # you try it out.
  
  print("Immigration")
  isect<-intersect(rownames(synlethNA),colnames(synlethNA)) # lol no need to remove Ts and Damp
  imigpop<-Init_Pop(members=isect,chrom_size=max(unlist(lapply(populus,length))),numbers=round((proportion*length(populus))))
  imigpop<-Fitness_score(imigpop,synlethNA)
  populus<-c(populus,imigpop)
  return(populus)
}

Remove_dupes<-function(populus){
  # Remove chromosomes with duplicate genes, I mean seriously, you are late!
  # 16-09-2015, just use apply, unique and remove indices which have depleted the length
  isect<-intersect(rownames(synlethNA),colnames(synlethNA))
  orig_len<-length(populus)
  chromlen<-unlist(lapply(populus,length))
  uqitpop<-lapply(populus,unique)
  uniqlen<-unlist(lapply(uqitpop,length))
  lendiff<-chromlen-uniqlen
  dupIx<-which(lendiff!=0)
  message(paste(length(dupIx),"duplicates removed"))
  if(length(dupIx)!=0)
  {
    populus<-populus[-c(dupIx)]
    populus<-c(populus,Init_Pop(members=isect,chrom_size=chromlen,numbers=(orig_len-length(populus))))
    populus<-Fitness_score(populus,synlethNA)
  }
  return(populus)
}

SingleGArun<-function(ChromNum=1000,iter=30,synlethNA,ChromLen=8,Crossover=0.5,ImmigRatio=0.01){
  
  #Argument List
  # ChromNum is the number of chromosomes in the population, runtime complexity changes signficantly
  # iter, the number of iterations that you believe should be until convergence occurs
  # ChromLen, the length of the chromosome/ no. of genes it contains.
  # Crossover is the proportion of the population to be considered for crossover
  # ImmigRatio is the proportion of population that must be added to whole population by Immigration
  
  #-------------- So data pre-processing specific to only my dataset-------------------#
  gene_names<-c(rownames(synlethNA),colnames(synlethNA)) # Take common set of all genes in experiment
  gnames_uq<-unique(gene_names) # Take the unique set out of all those genes
  damp_genes<-grep("damp",x=gnames_uq) # grep out the damp genes
  tsq_genes<-grep("tsq",gnames_uq) # grep out the ts alleles
  all_genes<-gnames_uq[-c(damp_genes,tsq_genes)] #all genes in exp with tsq and damp alleles removed
  isect<-intersect(colnames(synlethNA),rownames(synlethNA)) # all common genes between rows and columns
  #------------------------------------------------------------------------------------# 
  #iter=30
  #ChromNum=1000
  #ChromLen=8
  #-----------Initial data declarations and variable storage and initializations-------#
  scoremax<-vector(mode="list",length=iter) # To store all the maximum scores
  initpop<-Init_Pop(members=isect,chrom_size=ChromLen,numbers=ChromNum) # init the first random draw
  initpop<-Remove_dupes(initpop) # Remove duplicate genes from draw event
  initpop<-Fitness_score(initpop,synlethNA) # Score the population
  SelPop<-Selection(initpop) # Apply positive selection and select a population to iterate over
  Elitist<-vector(mode="list",iter) # To store the elite members of the population
  WholePop<-vector(mode="list")
  #-----------------------------------------------------------------------------------#
  
  #  For the first random draw of genes into SelPop (awful and chance games)
  #  enter iterator retrieve fitness scores for the members of SelPop
  #  From SelPop pick one member with the highest score store into Elite list with index=iter
  #  Splice SelPop to SelpopSplice, take (50%) of population, intermix chromosomes
  #  Concatenate spliced population SelPopSplicee and SelPop to TransGenPop(just in case FUBAR)
  #  Selection of good stuff from the mixed population into TransGenPopSel
  #  Mild programming inconsistency send whole population to Immig which returns immigrant + original population as TransGenImmig
  #  Do a Selection on TransGenImmig and return the population to SelPop to iter RINSE and REPEAT
  for(i in 1:iter){
    scoremax[[i]]<-as.numeric(names(SelPop)) # The score store
    WholePop[[i]]<-SelPop # Storing ALL the random gene selections
    Elitist[[i]]<-SelPop[[which(scoremax[[i]]==max(scoremax[[i]]))[1]]] # Store the first of the maximum scoring in the Elite population list
    SelpopSplice<-SpliceR(SelPop,synlethNA,proportion = Crossover) #SpliceR takes some fraction of the population and internally splices them randomly
    SelpopSplice<-Remove_dupes(SelpopSplice) # In case splicing generates duplicates; remove
    TransGenPop<-c(SelPop,SelpopSplice) # Add the spliced population to the original
    TransGenPop<-Remove_dupes(TransGenPop) # remove duplicates if any (unlikely)
    TransGenPopSel<-Selection(TransGenPop,PopCount=ChromNum) # apply a round of selection
    TransGenImmig<-Immig(TransGenPopSel,synlethNA,proportion = ImmigRatio) # Immigrate new sampled members from population
    TransGenImmig<-Remove_dupes(TransGenImmig) # remove duplicates
    SelPop<-Selection(TransGenImmig,PopCount=ChromNum) # another round of selection.
    print(i)
  }
  
  # Data to be returned by function
  # It is a list of structure.
  # SelPop which is the complete population by the end of the run
  # Elitist, which is a list of the best chromosome at each iteration
  # IterScoreHist is the maximum score obtained at each iteration.
  SGARun<-list()
  SGARun$SelPop<-WholePop # Whole Pop is a list of lists(chromosomes)
  SGARun$Elitist<-Elitist
  SGARun$IterScoreHist<-scoremax
  
  return(SGARun)
}
