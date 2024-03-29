
# R package for "Identifying Disease-Associated Biomarker Network Features Through Conditional Graphical Model"

<img src="https://img.shields.io/badge/Study%20Status-Results%20Available-yellow.svg" alt="Study Status: Results Available"> 

Biomarkers are often organized into networks, in which the strengths of network connections vary across subjects depending on subject‐specific covariates (eg, genetic variants). Variation of network connections, as subject‐specific feature variables, has been found to predict disease clinical outcome. In this work, we develop a two‐stage method to estimate biomarker networks that account for heterogeneity among subjects and evaluate network's association with disease clinical outcome. 

- Title: **Identifying Disease-Associated Biomarker Network Features Through Conditional Graphical Model**
<br/> [Manuscript](https://onlinelibrary.wiley.com/doi/10.1111/biom.13201): Biometrics 76 (3): 995-1006

- Authors: **Shanghong Xie<sup>a</sup> (sx2168@columbia.edu), Xiang Li<sup>b</sup>,  Peter McColgan<sup>c</sup>,  Rachael I. Scahill<sup>c</sup>,  Donglin Zeng<sup>d</sup>, and Yuanjia Wang<sup>a,e</sup>**

- Affiliations: 
   + 1. **Department of Biostatistics, Mailman School of Public Health, Columbia University, New York**
   + 2. **Statistics and Decision Sciences, Janssen Research & Development, LLC, Raritan, New Jersey**
   + 3. **Huntington’s Disease Centre, Department of Neurodegenerative Disease, UCL Institute of Neurology, London, UK**
   + 4. **Department of Biostatistics, University of North Carolina, Chapel Hill, North Carolina**
   + 5. **Department of Psychiatry, Columbia University Medical Center, New York**
  
<p align="center">
<img src="https://github.com/shanghongxie/Covariate-adjusted-network/blob/master/outline.png" width="500" height="600">
</p>



## Setup Requirements
- R
- Install Rcpp and RcppEigen packages

## Code Instructions

The code for the proposed methodology is included in **cNetworkC.cpp** and **cNetworkR.R**. Source the two files to implement the method. The arguments are described inside the code. 

sourceCpp("cNetworkC.cpp")

source("cNetworkR.R") 

### Main functions

- **LmNetwork_1st**: 1st stage of the proposed approach to estimate subject-specific networks.
 
- **EnetLm**: 2nd stage of the proposed approach to assess the assocation between subject-specific network connections with the clinical outcome.

### Examples
- **SimGenerate.R** includes the code to generate simulation data.

- **SimSample.R** provides an example of simulation study.
