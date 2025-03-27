# Shiny for Olink Data Analysis
(Under active development)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14287630.svg)](https://doi.org/10.5281/zenodo.14287630)
[![GitHub Invite Collaborators](https://img.shields.io/badge/Invite-Collaborators-blue?style=for-the-badge&logo=github)](https://github.com/JD2112/ShinyWGCNA/settings/access)

This Shiny app provides a user-friendly interface for performing quick analysis of Olink data using [OlinkAnalyze package](https://github.com/Olink-Proteomics/OlinkRPackage). It allows for data upload, parameter customization, visualization of results, and exporting of analysis outputs.

## Data Upload and Merging

NOTE:
- needs to change the separator from `;` to `,`.
- should not be any space in the column name.
- should not any special character or Swedish characters (Ö,Ä,Å) in the column name.

In the sidebar, you'll find three file upload buttons:
  - Upload NPX Data (CSV) - needs to have **SampleID** column
  - Upload Key File (CSV) - needs to have **SampleID** and **SUBJID** columns
  - Upload Variables File (CSV) - needs to have **SUBJID** columns

## How to run locally

### Using R
```
R # laod R >= 4.3.1
library(shiny)
shiny::runApp("app.R")
```

### Using Docker
`docker run --rm -p 3838:3838 jd21/shinyolink:v1.2.1`

### Online webserver

[https://olinkwrapper.serve.scilifelab.se/](https://olinkwrapper.serve.scilifelab.se/)

## Highlights

1. Data upload and merging
2. Data preview
3. Descriptive statistics
4. Normality tests for individual proteins
5. PCA plot with customizable options
6. T-Test with downloadable results
7. ANOVA with customizable number of covariates and downloadable results
8. Volcano plot
9. Violin plot for selected proteins


## Credits
- Main Author: 
    - Jyotirmoy Das ([@JD2112](https://github.com/JD2112))

- Collaborators: 
    - Debojyoti Das ([@BioDebojyoti](https://github.com/BioDebojyoti))

## References
See OlinkAnalyze package [https://github.com/Olink-Proteomics/OlinkRPackage](https://github.com/Olink-Proteomics/OlinkRPackage)

## Citation

Jyotirmoy Das. (2024). JD2112/ShinyOlink: ShinyOlink (v1.0). Zenodo. https://doi.org/10.5281/zenodo.14287630

## Acknowledgement

I would like to acknowledge the **Core Facility, Faculty of Medicine and Health Sciences, Linköping University, Linköping, Sweden** and **Clinical Genomics Linköping, Science for Life Laboratory, Sweden** for their support to build this application and SciLifeLab data center to host the application on their server.

