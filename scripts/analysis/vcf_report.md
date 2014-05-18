VCF Compare Script
==================


```r
dir()
```

```
##  [1] "analysis.Rproj"    "compare_mmp.R"     "compare_vcf.R"    
##  [4] "figure"            "gtcheck.R"         "plot_stats.py"    
##  [7] "read_depth.R"      "sync_vcf.local.sh" "vcf_report.html"  
## [10] "vcf_report.md"     "vcf_report.Rmd"
```

```r
getwd()
```

```
## [1] "/Users/daniel/Documents/git/Variant-Caller-Pipeline/scripts/analysis"
```






```r
print(args)
```

```
## [1] "../../data/vcf/03_RET1b.txt.Q40.vcf.gz"
## [2] "../../data/vcf/03_RET1.txt.Q40.vcf.gz" 
## [3] "../../data/vcf/03_RET1a.txt.Q40.vcf.gz"
```











## Individual VCF Results
![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-41.png) ![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-42.png) ![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-43.png) ![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-44.png) 


![plot of chunk PSC](figure/PSC1.png) ![plot of chunk PSC](figure/PSC2.png) 


## Ind. Sample Concordance #

![plot of chunk ind_conc](figure/ind_conc.png) 


## Pairwise Concordance

[[1]]
![plot of chunk pairwise_con](figure/pairwise_con1.png) 
[[2]]
![plot of chunk pairwise_con](figure/pairwise_con2.png) 

