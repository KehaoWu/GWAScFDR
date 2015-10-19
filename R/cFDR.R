cFDR = function(p_i,p_j) {

  # Compute empirical FDR at remaining points
  
  cdf = rep(0,length(p_i))
  
  oj = order(p_j)
  p_i1 = p_i[oj]; p_j1 = p_j[oj] #
  or = order(oj) # index of snps in R in p_i1/p_j1

  pb = txtProgressBar(min = 1,max = length(p_i1),style=3)
  for (i in 1:length(p1)) {
    cdf[i] = length(which(p_i1[1:or[i]]<=p_i1[or[i]]))/or[i]
    setTxtProgressBar(pb = pb,value = i)
  }
  
  
  f_i=p1/cdf
  
  return(f_i)
}

GenomicControl = function(p){
  z = qnorm(1 - p / 2)
  lambda = median(z ^ 2 ) / 0.456
  print(lambda)
  zad = sqrt(z^2/lambda)
  2 * pnorm(zad,lower.tail = F)
}
