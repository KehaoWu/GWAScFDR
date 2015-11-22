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

manhattanPlot = function(pvalue,bp,chr,gene=NULL,
                         ylab=expression(-log[10](P)),
                         cutoffline=T,
                         chrLabel = 1:22){
  if(grepl(pattern = "chr",x = chr[1])){
    chr = as.numeric(gsub(pattern = "chr",replacement = "",x = chr))
  }
  bp = bp[pvalue!=0]
  chr = chr[pvalue!=0]
  pvalue = pvalue[pvalue!=0]
  pvalue = pvalue[!is.na(chr)]
  bp = bp[!is.na(chr)]
  chr = chr[!is.na(chr)]
  pvalue = -log10(pvalue)
  pvalue = ifelse(pvalue<0,0,pvalue)
  print(max(pvalue))
  bpMidVec <- vector(length=length(chrLabel))
  maxbp = 0
  bmin = vector(length=length(chrLabel))
  bmax = vector(length=length(chrLabel))
  for(i in chrLabel){
    bp[chr==i] = bp[chr==i] + maxbp
    bmin[i] = min(bp[chr==i])
    bmax[i] = max(bp[chr==i])
    bpMidVec[i] <- ((max(bp[chr==i]) - min(bp[chr==i]))/2) + min(bp[chr==i])
    maxbp = max(bp[chr==i])
  }
  chr = factor(chr)
  yLabel = round(c(0,1,-log10(0.05),2:(max(pvalue)[1])),digits = 2)
  p = ggplot() + 
    geom_rect(data = data.frame(bmin,bmax,alpha=0.01),
              aes(xmin=bmin,xmax=bmax,alpha=alpha),
              ymin=0,
              ymax=Inf,
              fill = rep(c("grey90","white"),11),
              size = 0
    ) + 
    geom_point(data = data.frame(P=pvalue,BP=bp,CHR=chr),
               aes(y=P,x=BP,colour=CHR),
               alpha=0.8) +
    ylim(0,1.3*max(pvalue)) +
    scale_x_continuous(labels=as.character(chrLabel), breaks=bpMidVec) +
    scale_y_continuous(labels=as.character(yLabel), breaks=yLabel) +
    scale_color_manual(values=rep(c('orange1', 'grey20'), 11)) +
    theme_bw() +
    theme(
      panel.grid=element_blank()
    ) +
    xlab("Chromosomal Location") +
    ylab(ylab) +
    theme(legend.position='none')
  if (cutoffline){
    p = p + geom_hline(y=-log10(0.05), linetype=1, col='red', lwd=1) 
  }
  if (!is.null(gene)){
    
    x = bp[pvalue>=-log10(0.05)]
    gene = gene[pvalue>=-log10(0.05)]
    y = pvalue[pvalue>=-log10(0.05)]
    p = p + 
      geom_text(data=data.frame(y=y,x=x,gene=gene),
                      aes(y=y,x=x,label=gene),
                      hjust=0) +
      geom_point(data=data.frame(y=y,x=x),
                aes(y=y,x=x,size=2),colour="red")
  }
  p
}

cFDRDotPlot = function(p1,p2,cFDR){
  x = -log10(p1)
  y = -log10(p2)
  z = -log10(unlist(lapply(cFDR,FUN = function(x)min(x,1))))
  p = ggplot(data=data.frame(x,y,z)) +
    geom_point(aes(x=x,y=y,color=z)) +
    ylim(0,7) + 
    xlim(0,7) +
    scale_colour_gradientn(colours=c("white","yellow","orange","tomato","red"),
                           values=rescale(c(0,1,2,4,max(z))), 
                           space = "Lab")
  plot(p)
}

QQplot = function(x,y = NULL){
  fs = factor(y)
  maxValue = xx = yy = NULL
  x = -log10(x)
  for(f in levels(fs)){
    maxValue = c(max(x[fs==f]))
    yyy = x[fs==f]
    yy = c(yy,sort(yyy,decreasing = T))
    xx = c(xx,-log10((1:sum(fs==f))/sum(fs==f)))
  }
  p = ggplot(data=data.frame(
      x = xx, y = yy, SNPs = fs
    )) +
    geom_line(aes(x=x,y=y,colour=SNPs)) +
    theme_bw() +
    theme(
      panel.grid=element_blank()
    ) +
    xlab(paste("Empirical",expression(-log10(q)))) +
    ylab(paste("Nominal",expression(-log10(p)))) +
    geom_abline(intercept = 0.8*min(maxValue),slope=0,linetype=2,colour="red")
  p
}



