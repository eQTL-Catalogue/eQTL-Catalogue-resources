---
title: "eQTL Catalogue coloc.susie example"
output: 
  html_document: 
    keep_md: yes
date: "2024-06-21"
---



## Perform colocalisation with coloc.susie and lbf_variable files from the QTL Catalogue

Load required pacakges:

```r
library("coloc")
```

```
## This is coloc version 5.2.1
```

```r
library("readr")
library("dplyr")
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

Extract BPI eQTL and pQTL LBF variables from eQTL Catalogue

```bash
#gunzip -c QTD000021.lbf_variable.txt.gz | head -n1 > QTD000021_ENSG00000101425.tsv && gunzip -c QTD000021.lbf_variable.txt.gz | grep ENSG00000101425 >> QTD000021_ENSG00000101425.tsv
#gunzip -c QTD000026.lbf_variable.txt.gz | head -n1 > QTD000026_ENSG00000101425.tsv && gunzip -c QTD000026.lbf_variable.txt.gz | grep ENSG00000101425 >> QTD000026_ENSG00000101425.tsv
#gunzip -c QTD000584.lbf_variable.txt.gz | head -n1 > QTD000584_BPI.tsv && gunzip -c QTD000584.lbf_variable.txt.gz | grep "BPI.4126.22.1..1" >> QTD000584_BPI.tsv
```

Import lbf variable values into R

```r
#BLUEPRINT monocytes
mono_lbf = readr::read_tsv("data/QTD000026_ENSG00000101425.tsv.gz", show_col_types = FALSE)

#BLUEPRINT neutrophils
neutrophil_lbf = readr::read_tsv("data/QTD000021_ENSG00000101425.tsv.gz", show_col_types = FALSE)

#INTERVAL plasma pQTLs
interval_lbf = readr::read_tsv("data/QTD000584_BPI.tsv.gz", show_col_types = FALSE)
```

Convert LBF data frames into matrixes suitable for the coloc.bf_bf() method


```r
#BLUEPRINT monocytes
mono_mat = as.matrix(dplyr::select(mono_lbf, lbf_variable1:lbf_variable10))
row.names(mono_mat) = mono_lbf$variant
mono_mat = t(mono_mat)

#BLUEPRINT neutrophils
neutro_mat = as.matrix(dplyr::select(neutrophil_lbf, lbf_variable1:lbf_variable10))
row.names(neutro_mat) = neutrophil_lbf$variant
neutro_mat = t(neutro_mat)

#INTERVAL plasma pQTLs
interval_mat = as.matrix(dplyr::select(interval_lbf, lbf_variable1:lbf_variable10))
row.names(interval_mat) = interval_lbf$variant
interval_mat = t(interval_mat)
```

### Perform all pairwise colocalisations

BLUEPRINT neutrophils eQTL vs INTERVAL plasma pQTL

```r
neutro_coloc = coloc::coloc.bf_bf(neutro_mat, interval_mat)
dplyr::as_tibble(neutro_coloc$summary) %>% dplyr::filter(PP.H4.abf > 0.8)
```

```
## # A tibble: 2 × 10
##   nsnps hit1       hit2  PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf  idx1
##   <int> <chr>      <chr>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl> <int>
## 1  6011 chr20_382… chr2…  8.04e-32  1.77e-10  2.49e-23    0.0528     0.947     1
## 2  6011 chr20_383… chr2…  1.67e-48  8.58e-45  2.57e- 6    0.0112     0.989     2
## # ℹ 1 more variable: idx2 <int>
```


BLUEPRINT monocytes eQTL vs INTERVAL plasma pQTL

```r
mono_coloc = coloc::coloc.bf_bf(mono_mat, interval_mat)
dplyr::as_tibble(mono_coloc$summary) %>% dplyr::filter(PP.H4.abf > 0.8)
```

```
## # A tibble: 4 × 10
##   nsnps hit1       hit2  PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf  idx1
##   <int> <chr>      <chr>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl> <int>
## 1  6015 chr20_382… chr2…  6.55e-15  7.10e-11  2.02e- 6    0.0200     0.980     2
## 2  6015 chr20_383… chr2…  1.79e-52  4.07e-44  2.75e-10    0.0609     0.939     1
## 3  6015 chr20_382… chr2…  6.73e-44  7.30e-40  1.48e- 5    0.159      0.841     2
## 4  6015 chr20_383… chr2…  2.49e-40  2.42e-40  5.48e- 2    0.0515     0.894     3
## # ℹ 1 more variable: idx2 <int>
```


BLUEPRINT monocytes eQTL vs BLUEPRINT neutrophils eQTL

```r
eQTL_coloc = coloc::coloc.bf_bf(neutro_mat, mono_mat)
dplyr::as_tibble(eQTL_coloc$summary) %>% dplyr::filter(PP.H4.abf > 0.8)
```

```
## # A tibble: 3 × 10
##   nsnps hit1       hit2  PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf  idx1
##   <int> <chr>      <chr>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl> <int>
## 1  6249 chr20_383… chr2…  9.31e-14  4.78e-10  2.12e- 5    0.107      0.893     2
## 2  6249 chr20_382… chr2…  1.62e-27  3.56e- 6  1.75e-23    0.0367     0.963     1
## 3  6249 chr20_383… chr2…  4.17e- 5  6.82e- 2  4.11e- 5    0.0655     0.866     3
## # ℹ 1 more variable: idx2 <int>
```

