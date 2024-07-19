# Imaging Mass Cytometry on GBM Tumour Samples

This repository contains R code for carrying out imaging mass cytometry on glioblastoma (GBM) tumor samples. The project uses `renv` for managing package dependencies to ensure reproducibility.

## Table of Contents

-   [Introduction](#introduction)
-   [Installation](#installation)
-   [Usage](#usage)

## Introduction {#introduction}

Glioblastoma (GBM) is an aggressive type of cancer that occurs in the brain or spinal cord. Imaging mass cytometry is a powerful technique that combines mass spectrometry and microscopy to analyse the spatial distribution of proteins and other molecules within tissue sections.

This project provides R scripts and data for performing imaging mass cytometry analysis on GBM tumour samples. The analysis includes data preprocessing, visualization, and statistical analysis.

## Installation {#installation}

### Prerequisites

-   R (version 4.4.1 or higher)
-   RStudio (recommended)

### Clone the Repository

Clone this repository to your local machine using the following command:

```         
sh git clone <https://github.com/yourusername/gbm-imaging-mass-cytometry.git> 

cd gbm-imaging-mass-cytometry
```

This project uses renv to manage R package dependencies. To initialize the renv environment and install the required packages, run the following commands in R:

```         
install.packages("renv")

renv::restore()
```

This will install all the necessary packages specified in the renv.lock file.

## Usage {#usage}

### Running the Analysis

Open the R project file (GBM_IMC_Analysis.Rproj) in RStudio. The main analysis scripts are located in the `process` directory. The order of the scripts are denoted by the prefix number in the file name. Run the scripts in numerical order to perform the analysis. Moreover, the downstream analysis scripts are located in the `Downstream` sub-directory.
