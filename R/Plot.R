
stratifiedQQplot = function(p1,p2){
  p1 = p1[p1!=0]
  p2 = p2[p2!=0]
  
  dat = NULL
  for(cutoff in c(1,0.01,0.001,0.0001)){
    p = p1[p2<=cutoff]
    x = -log10(seq(from = 0,to = 1,length.out = length(p)+1))[-1]
    print(max(x))
    y = sort(-log10(p),decreasing = T)
    dat = rbind(dat,data.frame(x=x,y=y,cutoff))
  }
  dat$cutoff = factor(dat$cutoff)
  p = ggplot(dat,aes(x=x,y=y,fill=cutoff,colour=cutoff)) +
    geom_line() +
    geom_abline(intercept=0,slope=1) +
    labs(title = "Stratified Q-Q Plot") +
    labs(x = "Nominal -log p conditional") +
    labs(y = "Nominal -log p")  
}

stratifiedTDRplot = function(p1,p2){
  
  dat = NULL
  for(cutoff in c(1,0.1,0.01,0.001)){
    p = p1[p2<=cutoff]
    y = 1 - unlist(lapply(p/ecdf(p)(p),FUN = function(X)min(X,1)))
    x = -log10(p)
    dat = rbind(dat,data.frame(x=x,y=y,cutoff))
  }
  dat$cutoff = factor(dat$cutoff)
  p = ggplot(dat,aes(x=x,y=y,fill=cutoff,colour=cutoff)) +
    geom_line() +
    labs(title = "Stratified True Discovery Rate Plot") +
    labs(x = "Nominal -log p") +
    labs(y = "TDR: 1 - FDR") 
}

stratifiedQQForGenomeControlplot = function(p){
  p = sort(p)
  y = -log10(GenomicControl(p))
  x = -log10(seq(from = 0,to = 1,length.out = length(p)+1)[-1])
  dat = data.frame(x=x,y=y,type="Adjusted")
  y = -log10(p)
  dat = rbind(dat,data.frame(x=x,y=y,type="raw"))
  p = ggplot(dat) + 
    geom_line(aes(x=x,y=y,colour=type)) + 
    geom_abline(intercept=0,slope=1)
  p
}


