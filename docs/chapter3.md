# 3. A Note on the App's Technical Stack
The Olink Wrapper App is built with a powerful combination of open-source and modern technologies to ensure it is both interactive and reproducible.

## 3.1 Shiny

**Shiny** is an open-source R package that provides a robust framework for building interactive web applications using R. It allows researchers and data scientists to create powerful applications with a live backend powered by R's analytical capabilities, without needing extensive web development knowledge. The app's user interface, plots, and statistical computations are all handled by the Shiny framework.

## 3.2 R Packages

The app relies on a variety of R packages to perform its functions. While the full list can be found in the DESCRIPTION file of the GitHub repository, key packages likely include:

- `shiny`: The core package for the web application framework.
- `dplyr`: For efficient data manipulation and cleaning.
- `ggplot2`: For generating high-quality, customizable data visualizations.
- `plotly`: For creating interactive plots that allow users to explore data.
- `olink`: For specific functions related to Olink data handling and analysis.
- `readxl`: For reading data from Excel files.

## 3.3 Docker

Docker is a platform that uses containerization to package an application and all its dependencies (including the R environment, Shiny, and all required packages) into a single, isolated "container." This container can be run on any system that has Docker installed. The use of Docker provides several key benefits:

- **Reproducibility:** The app will run exactly the same way every time, regardless of the host system's configuration.
- **Simplified Deployment:** Deployment is as simple as building and running the Docker image, eliminating complex dependency management issues.
- **Isolation:** The app runs in its own environment, preventing conflicts with other applications or system configurations.

By using Docker, the app ensures that everyone from the developer to an end-user has a consistent and reliable experience, making the analysis and results fully reproducible.