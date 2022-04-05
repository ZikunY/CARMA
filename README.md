# CARMA
## Novel Bayesian model for fine-mapping in meta-analysis studies


We propose a novel Bayesian model, CARMA (CAusal Robust mapping method in Meta-Analysis studies), for fine-mapping in order to identify putative causal variants at GWAS loci. The main features of CARMA are

> Replacing the Normal-Gamma prior family for the effect size distribution by a heavy-tail Cauchy distribution to achieve higher power
>
> Modeling jointly summary statistics and high-dimensional functional annotations
>
> Introducing a novel Bayesian hypothesis testing approach to account for discrepancies between summary statistics and LD from external reference panels

### Vignettes
Vignettes can be found at [https://github.com/ZikunY/CARMA/blob/master/CARMA_demo.pdf].

### Package manual
The manual can be found at [https://github.com/ZikunY/CARMA/blob/master/CARMA_1.0.pdf].

### Installation
```r
devtools::install_github("ZikunY/CARMA")
```

### Citation

> Yang, Z.,  Wang C, Khan A, Vardarajan B,  Mayeux R,  Kiryluk D, Ionita-Laza I, 
> "CARMA: Novel Bayesian model for fine-mapping in meta-analysis studies"
