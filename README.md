# CARMA
## Novel Bayesian model for fine-mapping in meta-analysis studies


We propose a novel Bayesian model, CARMA (CAusal Robust mapping method in Meta-Analysis studies), for fine-mapping in order to identify putative causal variants at GWAS loci. The main features of CARMA are

> Modeling jointly summary statistics and high-dimensional functional annotations
>
> Introducing a novel Bayesian hypothesis testing approach to account for discrepancies between summary statistics and LD from external reference panels
>
> flexible specification of the prior distribution of effect sizes including a heavy-tail Cauchy distribution that is more adaptive to the possibly different signal-to-noise ratios across loci

### Vignettes
Vignettes can be found at [https://github.com/ZikunY/CARMA/blob/master/CARMA_demo.pdf].

### Package manual
The manual can be found at [https://github.com/ZikunY/CARMA/blob/master/CARMA.pdf].

### Installation
```r
devtools::install_github("ZikunY/CARMA")
```
### Intel MLK library 
CARMA requires the Intel MLK library, therefore please specify the pathway for the Intel MLK library. The Intel MLK library can be installed by the following lines

```bash
conda install -c anaconda mkl
LD_PRELOAD=pathto/anaconda3/lib/libmkl_rt.so R
```
Important Note

The pathway for the Intel MLK library and related libraries may differ depending on the operating system. Users should verify the correct pathway on their system and modify the LD_PRELOAD accordingly.

### Citation

> Yang, Z.,  Wang C, ， Khan A, Vardarajan B,  Mayeux R,  Kiryluk D, Ionita-Laza I, 

> CARMA: Novel Bayesian model for fine-mapping in meta-analysis studies
