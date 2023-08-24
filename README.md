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

### Implementation on M1 Mac
If users encounter issues while installing CARMA on an M1 Mac, it may be due to the failure to locate the GSL library or other libraries during the installation process. The problem could be solved by modifying the ‘PKG_CPPFLAGS’ environment variable in the `/.R/Makevars file’. This file should be created if it doesn’t exist. Open the `/.R/Makevars’ file in a text editor, and add the following line to set the PKG_CPPFLAGS environment variable:

```bash
PKG_CPPFLAGS=-I/opt/homebrew/Cellar/gsl/2.7.1/include/
```
The pathway above corresponds to the GSL header files’ location for GSL version 2.7.1. Depending on the user’s system, this pathway may vary.

Setting PKG_LIBS

Next, add the necessary lines to the ~/.R/Makevars file to specify the GSL library’s location and other required libraries:
```bash
PKG_LIBS=-L/opt/homebrew/lib -lgsl -lgslcblas -lm
```
The above line tells the linker (ld) to search for the GSL library (libgsl) and GSL CBLAS (libgslcblas) in the Homebrew installation directory (/opt/homebrew/lib). The -lm flag links the math library.


### Citation

> Yang, Z.,  Wang C, ， Khan A, Vardarajan B,  Mayeux R,  Kiryluk D, Ionita-Laza I, 

> CARMA: Novel Bayesian model for fine-mapping in meta-analysis studies
