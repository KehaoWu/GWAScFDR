library(ggplot2)

stratifiedQQplot = function(p1,p2,xlab="Nominal -log p conditional",ylab="Nominal -log p"){
  library(ggplot2)
  p1 = p1[p1!=0]
  p2 = p2[p2!=0]
  
  dat = NULL
  for(cutoff in c(1,0.1,0.01,0.001,0.0001)){
    p = p1[p2<=cutoff]
    x = -log10(seq(from = 0,to = 1,length.out = length(p)+1))[-1]
    print(max(x))
    y = sort(-log10(p),decreasing = T)
    dat = rbind(dat,data.frame(x=x,y=y,cutoff))
  }
  dat$cutoff = factor(dat$cutoff)
  p = ggplot(dat,aes(x=x,y=y,fill=cutoff,colour=cutoff)) +
    geom_line(size=1.2) +
    geom_abline(intercept=0,slope=1) +
    labs(title = "Stratified Q-Q Plot") +
    labs(x = xlab) +
    labs(y = ylab) +
    ylim(0,log10(length(p1)))
}

stratifiedTDRplot = function(p1,p2){
  library(ggplot2)
  
  dat = NULL
  for(cutoff in c(1,0.1,0.01,0.001)){
    p = p1[p2<=cutoff]
    y = p * length(p) / (order(order(p)))
    y = 1- ifelse(y>=1,1,y)
    x = -log10(p)
    dat = rbind(dat,data.frame(x=x,y=y,cutoff))
  }
  dat$cutoff = factor(dat$cutoff)
  p = ggplot(dat,aes(x=x,y=y,fill=cutoff,colour=cutoff)) +
    geom_line() +
    labs(title = "Stratified True Discovery Rate Plot") +
    labs(x = "Nominal -log p") +
    labs(y = "TDR: 1 - FDR")
  p
}

stratifiedQQForGenomeControlplot = function(p){
  library(ggplot2)
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

manhattanPlot = function(pvalue,bp,chr,ylab=expression(Conditional -log[10](FDR))){
  pvalue = -log10(pvalue)
  chrLabel = 1:22
  bpMidVec <- vector(length=length(chrLabel))
  maxbp = 0
  bmin = vector(length=length(chrLabel))
  bmax = vector(length=length(chrLabel))
  for(i in chrLabel){
    print(i)
    bp[chr==i] = bp[chr==i] + maxbp
    bmin[i] = min(bp[chr==i])
    bmax[i] = max(bp[chr==i])
    bpMidVec[i] <- ((max(bp[chr==i]) - min(bp[chr==i]))/2) + min(bp[chr==i])
    maxbp = max(bp[chr==i])
  }
  chr = factor(chr)
  print(chrLabel)
  print(bpMidVec)
  p = ggplot() + 
    geom_rect(data = data.frame(bmin,bmax,alpha=0.01),
              aes(xmin=bmin,xmax=bmax,alpha=alpha),
              ymin=0,
              ymax=Inf,
              fill = rep(c("grey90","white"),11),
              size = 0
    ) + 
    geom_point(data = data.frame(P=pvalue,BP=bp,CHR=chr),aes(y=P,x=BP,colour=CHR),alpha=0.8) +
    scale_x_continuous(labels=as.character(chrLabel), breaks=bpMidVec) +
    geom_hline(y=-log10(0.05), linetype=1, col='red', lwd=1) +
    scale_color_manual(values=rep(c('orange1', 'grey20'), 11)) +
    theme_bw() +
    theme(
      panel.grid=element_blank()
    ) +
    ylim(0,1.1*max(pvalue)) +
    xlab("Chromosomal Location") +
    ylab(ylab) +
    theme(legend.position='none')
  p
}

