# GWAS_cFDR

## Installation
  
```R  
  install.packages("devtools")
  devtools::install_github("KehaoWu/GWAScFDR")
```
  
## Demo

```R
  p1 = runif(10000,0,7)
  p2 = runif(10000,0,7)
  cFDR(1,p1,p2)   #Compute the 1st cFDR
  stratifiedQQForGenomeControlplot(p1)
  stratifiedQQplot(p1,p2)
  stratifiedTDRplot(p1,p2)
```
