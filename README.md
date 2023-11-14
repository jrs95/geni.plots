# GENI plots <img src='man/figures/geni.png' align="right" height="139"/>

The GENI plots R package is designed to create a series of plots for visualising results from genome-wide association studies (GWAS) and phenome-wide association studies (PheWAS).  

## Installation
```
install.packages("remotes")
remotes::install_github("jrs95/geni.plots", build_vignettes = TRUE)
```

## Functions
* `fig_manhattan`: creates a Manhattan plot for genomic markers from across the genome, e.g. results from genome-wide association studies.  
* `fig_phewas`: creates a plot visualising results from phenome-wide association studies (PheWAS).  
* `fig_qq`: creates a quantile-quantile (QQ) plot.  
* `fig_region`: creates a regional plot, i.e. a scatter graph of genomic markers associations (e.g. log10(p-values)) with a gene bar underneath.  
* `fig_stack_region`: creates a stacked regional association plot.  

## Example
See vignette for examples: https://jrs95.github.io/geni.plots/  

## Citation
Please cite this R package using the link: https://github.com/jrs95/geni.plots  
