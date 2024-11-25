<p align="right">
  <img src="https://github.com/ISGLOBAL-Rakislova-Lab/HTGAnalyzer_shiny/blob/main/www/HTGAnalyzer_logo.png" alt="HTGAnalyzer Logo" width="150">
</p>

# Introduction

This repository contains the Shiny app for the **HTGAnalyzer** package. If you are looking for the R package, you can find it in the following link: [HTGAnalyzer R Package](https://github.com/ISGLOBAL-Rakislova-Lab/HTGAnalyzer/blob/main/README.md).

The **HTGAnalyzer** package is designed to facilitate the analysis of HTG EdgeSeq and RNA sequencing data. It includes the **HTG_auto** function, which automates several key analyses, such as quality control, differential gene expression analysis, tumor microenvironment profiling, survival analysis, and gene set enrichment analysis.

For those who need more control over the analyses or wish to perform the full package analysis, you can install and use the complete **HTGAnalyzer** R package.

On the other hand, if you only need to recognize outliers, perform quality control (QC), or conduct statistical analysis, you can use the Shiny app for **HTGAnalyzer** available at: [HTGAnalyzer Shiny App](https://isglobal-rakislova-lab.shinyapps.io/htganalyzer_shiny/).

Feel free to use the Shiny app for easier interaction with the QC process and statistical analysis or install the full **HTGAnalyzer** package if you need to perform more comprehensive analyses.

# INSTALLATION INSTRUCTIONS
The HTGAnalyzer package uses the renv package to ensure that all users have the same package versions as used during development.
It is essential to have the HTGAnalyzer package installed to perform a complete analysis. If you are not familiar with R, 
we have provided a user-friendly script for installation. Follow these steps:

# 1. Install the R program
Download and install R software from the following link: [R version 4.4.2](https://cran.r-project.org/bin/windows/base/)
Make sure to install version "4.4.2" of R for compatibility.
- This setup works with both R and RStudio, so if you have both programs installed on your system, you can use either.

We recommend using **RStudio** for a better visual interface, which enhances the user experience, but the code will work fine in **R** as well.

To install RStudio, visit the following link: [Download RStudio](https://posit.co/download/rstudio-desktop/)
This page also offers the option to install R along with RStudio, which is convenient for most users.


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
# Activate library
library("renv")

# Download the file renv.lock
if (!file.exists("renv.lock")) {
  download.file(
    url = "https://raw.githubusercontent.com/ISGLOBAL-Rakislova-Lab/HTGAnalyzer_shiny/main/renv.lock",
    destfile = "renv.lock"
  )
}

# Install the packages
renv::restore(lockfile = "renv.lock")
## SELECT OPTION 1.
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
We recommend that, in addition to using the Shiny app, you also keep an eye on the **R console** (or **RStudio console** if you're using RStudio). The console will provide complementary information, such as details on which processes are being executed, the steps completed, and any errors that may occur. This is helpful for troubleshooting and understanding the flow of the analysis.


<div style="display: flex; justify-content: center; gap: 20px;">
    <img src="https://github.com/user-attachments/assets/25dbac67-84eb-4c58-af88-b7e67fdaec33" alt="Image 1" width="200"/>
    <img src="https://github.com/user-attachments/assets/ecf90cb3-9d11-46f7-8b63-cc5c3596902d" alt="ISGlobal Logo" width="200"/>
    <img src="https://github.com/user-attachments/assets/e2680f9a-38e4-4966-bb66-741d2cf58391" alt="Image 2" width="200"/>
</div>

