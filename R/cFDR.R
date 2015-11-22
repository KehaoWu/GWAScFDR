cFDR = function(p1,p2) {
  p_vals_ij = data.frame(p1,p2)
  f_i = rep(0,dim(p_vals_ij)[1])
  
  if (F) {
    # Compute cutoff c_i
    p_i = p_vals_ij[,1]; p_j = p_vals_ij[,2];
    ord=order(p_i); cdf=(1:length(ord))/length(ord) # Quantiles of ordered list of p values
    p_i = p_i[ord];
    fdr1 = p_i/cdf;
    c_i = p_i[min(which(fdr1>0.1))]
    
    r1 = which(log(p_vals_ij[,1])+log(p_vals_ij[,2])<0.5*log(c_i)) # values not in triangle (p_i,p_j<1; -log10(p_i) + -log10(p_j) <-log10(c_i))
    r2 = which(p_vals_ij[,1]<(10^-0.75)) # 
    #r3 = which(p_vals_ij[,1]>(10^-12)) # SNPs with values lower than this 
    r = intersect(r1,r2) # intersect(intersect(r1,r2),r3) # Values with ambiguous FDR
    
    total = 1:length(f_i);
    f_i[setdiff(total,r1)]=1; # approximately
    f_i[setdiff(total,r2)]=1; # approximately
    #f_i[setdiff(total,r3)]=0;
    
  } else {
    r = 1:dim(p_vals_ij)[1]
  }
  
  p_i = p_vals_ij[,1]; # Remaining points
  p_j = p_vals_ij[,2];
  
  
  # Compute empirical FDR at remaining points
  
  cdf = rep(0,length(p_i))
  
  oj = order(p_j)
  p_i1 = p_i[oj]; p_j1 = p_j[oj] #
  or = order(oj)[r] # index of snps in R in p_i1/p_j1
  cr = rep(0,length(r))
  pb = txtProgressBar(min = 0,max = length(r),style = 3)
  
  for (i in 1:length(r)) {
    setTxtProgressBar(pb = pb,value = i)
    cdf[r[i]] = length(which(p_i1[1:or[i]]<=p_i1[or[i]]))/or[i]
  }
  cat("\n")
  # Equivalent technique
  
  #for (i in 1:length(r)) {
  # cdf[r[i]] = length(which(p_i<=p_i[r[i]] & p_j<=p_j[r[i]]))/length(which(p_j<=p_j[r[i]]))
  #}
  
  f_i[r]=p1[r]/cdf[r] # Computation of FDR
  
  return(f_i)
}

GenomicControl = function(p,pintergenic){
  z = qnorm(1 - p / 2)
  zintergenic = qnorm(1 - pintergenic / 2)
  lambda = median(zintergenic ^ 2 ) / 0.456
  print(lambda)
  zad = z/sqrt(lambda)
  2 * pnorm(zad,lower.tail = F)
}


