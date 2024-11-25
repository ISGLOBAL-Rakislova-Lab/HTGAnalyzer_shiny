<p align="right">
  <img src="https://github.com/ISGLOBAL-Rakislova-Lab/HTGAnalyzer_shiny/blob/main/www/HTGAnalyzer_logo.png" alt="HTGAnalyzer Logo" width="150">
</p>

# INSTALLATION INSTRUCTIONS
The HTGAnalyzer package uses the renv package to ensure that all users have the same package versions as used during development.
It is essential to have the HTGAnalyzer package installed to perform a complete analysis. If you are not familiar with R, 
we have provided a user-friendly script for installation. Follow these steps:

# 1. Install the R program
Download and install the R software from the following link: https://cran.r-project.org/bin/windows/base/
Make sure to install version "4.4.2" of R for compatibility.

# 2. Install necessary packages
Run the following commands to install the renv and shiny packages if they are not already installed:
```{r}
install.packages("renv")
install.packages("shiny")
```

# 3. Restore the package environment
Use the renv package to ensure all dependencies are installed in their correct versions.
The renv.lock file from the GitHub repository will guide the installation.
```{r}
library("renv")
renv::restore(lockfile = "https://github.com/ISGLOBAL-Rakislova-Lab/HTGAnalyzer_shiny/blob/main/renv.lock")
```
# SHINY APP FOR QUALITY CONTROL (QC)
If you only need to perform Quality Control (QC) of the transcriptomic data, you can use the online Shiny app:
* https://isglobal-rakislova-lab.shinyapps.io/htganalyzer_shiny/

# COMPLETE ANALYSIS WITH HTGAnalyzer PACKAGE
The HTGAnalyzer package is an R package designed for the analysis of HTG EdgeSeq data and RNA sequencing results.
It includes tools for:
- Quality Control (QC)
- Differential gene expression analysis
- Tumor microenvironment profiling
- Survival analysis
- Gene set enrichment analysis (GSEA)

# USAGE INSTRUCTIONS
Once you have R installed and the necessary packages restored, you can run the Shiny app locally by using the following commands:

# 1. Load the shiny library
```{r}
library(shiny)
```

# 2. Run the Shiny app directly from GitHub
```{r}
runGitHub(repo = "HTGAnalyzer_shiny", username = "ISGLOBAL-Rakislova-Lab")
```
