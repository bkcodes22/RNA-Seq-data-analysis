# RNA-Seq Analysis App

This is an application built using Python, R, Flask, and Shiny for performing RNA-Seq analysis. The app integrates the following components:

- **Python**: Used to call the STAR alignment tool for aligning RNA-Seq reads.
- **Flask**: Provides a frontend for interacting with the app.
- **R Shiny App**: Used to process the counts file and generate Differentially Expressed Genes (DEG) results.

## Features

- **STAR Alignment**: Aligns RNA-Seq data using the STAR tool through Python.
- **Differential Expression Analysis**: Performs differential gene expression analysis on RNA-Seq counts using R and generates the DEG results in a Shiny app.
- **User-friendly Interface**: The frontend, built with Flask, allows users to easily input their data and view the results.

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/repositoryname.git

2. Install required Python packages:

   ```bash
   pip install -r requirements.txt

3. Install required R packages for Shiny and DEG analysis:
   ```R
   install.packages("shiny")
   install.packages("edgeR")
