<meta name="google-site-verification" content="rp1Onw6xf90nsr5dOkzDoYzC_iFLqBVtHRbH5riFcec" />
<p align="right">
  <img src="https://github.com/ISGLOBAL-Rakislova-Lab/HTGAnalyzer_shiny/blob/main/www/HTGAnalyzer_logo.png" alt="HTGAnalyzer Logo" width="150">
</p>

# Introduction

This repository contains the Shiny app for the **HTGAnalyzer** package. If you are looking for the R package, you can find it in the following link: [HTGAnalyzer R Package](https://github.com/ISGLOBAL-Rakislova-Lab/HTGAnalyzer/blob/main/README.md). The **HTGAnalyzer Shiny app** is designed to facilitate the quality control and transcriptomic analysis of HTG EdgeSeq and RNA sequencing data with minimal bioinformatical knowledge. It uses the **HTG_auto** function from the **HTGAnalyzer package**, which automates: 
* Quality control (QC) 
* Differential gene expression analysis (DEA)
* Tumor microenvironment profiling (TME)
* Survival analysis 
* Gene set enrichment analysis (GSEA)

**NOTE:** If you only need to recognize outliers, perform quality control (QC), or conduct statistical analysis, you can also use the Shiny app for **HTGAnalyzer** available at: [HTGAnalyzer Shiny App](http://isglobal-rakislova-lab.shinyapps.io/htg_shinny_app-main).


# INSTALLATION INSTRUCTIONS
The HTGAnalyzer package uses the renv package to ensure that all users have the same package versions as used during development.
It is essential to have the HTGAnalyzer package installed to perform a complete analysis. If you are not familiar with R, 
we have provided a user-friendly script for installation. Follow these steps:

## 1. Install R
Download and install the R software from the following link:  [R version 4.4.2](https://cran.r-project.org/bin/windows/base/)

* Ensure that you install version **"4.4.2"** of R for compatibility.
* This setup is compatible with both **R** and **RStudio**. If you have both installed, you can use either program to run the code.
* We recommend using RStudio for its enhanced visual interface, which improves the user experience. However, the code will work perfectly fine in R as well.

To download and install **RStudio**, use the following link:[Download RStudio](https://posit.co/download/rstudio-desktop/)

This page also offers an option to install R alongside RStudio, which is convenient for most users.

## 2: Install necessary R packages
Instructions for installing the HTGAnalyzer package and the Shiny app can be found here:
* [HTGAnalyzer Installation Guide](https://github.com/ISGLOBAL-Rakislova-Lab/HTGAnalyzer/tree/main)
  
Once you have R installed and the necessary packages restored, you can start using the Shiny app locally.

# USAGE INSTRUCTIONS
After successfully  installed, copy and paste the following command in R to run the Shiny app directly from Github.
```{r}
library(shiny)
runGitHub(repo = "HTGAnalyzer_shiny", username = "ISGLOBAL-Rakislova-Lab")
```
We recommend that, in addition to using the Shiny app, you also keep an eye on the **R console** (or **RStudio console** if you're using RStudio). The console will provide complementary information, such as details on which processes are being executed, the steps completed, and any errors that may occur. This is helpful for troubleshooting and understanding the flow of the analysis.

HTGAnalyzer is a registered digital work under intellectual property law.  
> © 2024 FUNDACIÓN PRIVADA INSTITUTO DE SALUD GLOBAL BARCELONA, FUNDACIÓ DE RECERCA CLÍNIC BARCELONA – INSTITUT D’INVESTIGACIONS BIOMÈDIQUES AUGUST PI I SUNYER, HOSPITAL CLÍNIC BARCELONA. All rights reserved. Registered Coloriuris 20/05/2025.


<div style="display: flex; justify-content: center; gap: 20px;">
    <img src="https://github.com/user-attachments/assets/74a92a2b-6007-4d37-b40c-d53fdf3ba0ef" alt="Image 1" width="200"/>
    <img src="https://github.com/user-attachments/assets/ecf90cb3-9d11-46f7-8b63-cc5c3596902d" alt="ISGlobal Logo" width="200"/>
    <img src="https://github.com/user-attachments/assets/e2680f9a-38e4-4966-bb66-741d2cf58391" alt="Image 2" width="200"/>
</div>
