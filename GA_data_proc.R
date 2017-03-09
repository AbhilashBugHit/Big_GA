MultiRunChomper<-function(MultiSimList,kya=NULL){
  if(length(kya)==0){message("Elite, All or IterScore")}
  if(kya=="Elite")
  {return(EliteOnly<-lapply(MultiSimList,function(x){return(x$Elitist)}))}
  
  if(kya=="All")
  {return(SelPopOnly<-lapply(MultiSimList,function(x){return(x$SelPop)}))}
  
  if(kya=="IterScore")
  {return(IterScore<-lapply(MultiSimList,function(x){return(x$IterScoreHist)}))}
}

BootPolish<-function(Bootno=100,FUN){
  # Will do the analysis done by multichomper iteratively on the bootstrap files
  Output<-vector(mode="list",length=Bootno)
  
  cust_FUN=match.fun(FUN)
  pb<-txtProgressBar(1,Bootno,style = 3)
  for(i in 1:Bootno)
  {
    setTxtProgressBar(pb,i)
    load(paste("Boot",i,".rda",sep=""))
    Output[[i]]<-cust_FUN(Run30x30)
  }
  return(Output)
}

Tabler<-function(x){
  AllSampled<-MultiRunChomper(x,"All")
  tab_revisit<-table(unlist(AllSampled))
  #hist(log10(tab_revisit))
  #plot(density(log10(tab_revisit)))
  return(tab_revisit)
}