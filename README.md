# GWAS_cFDR

## Installation
  
```R  
install.packages("devtools")
devtools::install_github("KehaoWu/GWAScFDR")
```
  
## Demo

```R
p1 = runif(10000,0,1)
p2 = runif(10000,0,1)
cFDR(p1,p2)
p = stratifiedQQplot(p1,p2)
plot(p)
```

```R
pvalue = runif(100)
bp = sample(10000:20000,100)
chr = sample(1:22,100,replace=T)
gene = as.character(1:100)
p = manhattanPlot(pvalue = pvalue,bp = bp,chr = chr,gene = gene)
plot(p)
```


