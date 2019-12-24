
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HTSet

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/zhuchcn/HTSet.svg?branch=master)](https://travis-ci.org/zhuchcn/HTSet)
[![codecov](https://codecov.io/gh/zhuchcn/HTSet/branch/master/graph/badge.svg)](https://codecov.io/gh/zhuchcn/HTSet)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

Quantitative high through-put experiments often generates data in very
similar structures, including an expression/abundance matrix, a table
for feature variables, and a table for sample metadata. The feature data
often holds feature characteristics such as gene ID, protein name,
metabolite annotation, or pathway ID. The sample metadata contains
variables for each sample, such as genotype, phenotype, or geographic
information. Thus in the **HTSet** package, the core S4 class `HTSet`
was created as a general data container for any quantitative high
through-put experiment data, including RNA-seq, microbiome, proteomics,
and metabolomics study.

## Installation

``` r
devtools::install_github("zhuchcn/HTSet")
```
