cFDR = function(i,p1,p2){

  if(length(p1)!=length(p2))
    return(0)
  
  N = length(p1)
  cat(i,"/",N,"\n")
  pset = p1[p2<=p2[i]]
  index = match(p1[i],pset)
  r = order(order(pset))
  res = p1[i] / (r[index] / N)
  return(res)
}

GenomicControl = function(p){
  z = qnorm(1 - p / 2)
  lambda = median(z ^ 2 ) / 0.456
  print(lambda)
  zad = sqrt(z^2/lambda)
  2 * pnorm(zad,lower.tail = F)
}
