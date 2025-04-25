# RNA-Seq Analysis App

This is an application built using Python, R, Flask, and Shiny for performing RNA-Seq analysis. The app integrates the following components:

- **Python**: Used to call the STAR alignment tool for aligning RNA-Seq reads and generate counts file.
- **Flask**: Provides a frontend for interacting with the app.
- **R Shiny App**: Used to process the counts file and generate Differentially Expressed Genes (DEG) results.

## Features

- **STAR Alignment**: Aligns RNA-Seq data using the STAR tool through Python.
- **Differential Expression Analysis**: Performs differential gene expression analysis on RNA-Seq counts using R and generates the DEG results in a Shiny app.
- **User-friendly Interface**: The frontend, built with Flask, allows users to easily input their data and view the results.

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/bkcodes22/RNA-Seq-data-analysis.git

2. Install required Python packages:

   ```bash
   pip install -r requirements.txt

3. Install required R packages for Shiny and DEG analysis:
   ```R
   install.packages("shiny")
   install.packages("edgeR")
   install.packages("shinyjs")
   install.packages("DT")
   install.packages("ggplot2")
   install.packages("FactoMineR")
   install.packages("RColorBrewer")
   install.packages("limma")

## Usage

1. Run the Flask app
   ```
   python main.py
   
2. Select the type of analysis you want to do (From Fastq files or directly through raw counts)

3. Give the path for the STAR executer, Input file / Directory and the Output Directory

4. BAM Files, Count files along with QC stats will be generated in the output directory

5. Run R shiny app, to generate DEG Results.

6. Our pipeline inbuilt uses TMM Normalisation and DEG analysis with limma package and allows you to save the results and plots

7. All the outputs will be generated in the designated directory
