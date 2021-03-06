# Predictors and Microbiology of Co-Infection in Patients with COVID-19: Living Rapid Review Update and Meta-Regression

## Purpose of this repository

This repository contains all code and data necessary to reproduce the analyses presented in the manuscript "[Predictors and microbiology of respiratory and bloodstream bacterial infection in patients with COVID-19: living rapid review update and meta-regression](https://doi.org/10.1016/j.cmi.2021.11.008)" by Langford et al (doi: [10.1016/j.cmi.2021.11.008](https://doi.org/10.1016/j.cmi.2021.11.008)). This includes every table and figure except for figure 1, figure 4, supplementary table 3 and supplementary table 4. Supplementary table 1 comprises a subset of the complete dataset found in the directory `data`. Note that the included scripts may produce figures and tables beyond what are directly presented in the manuscript.

For more details on the methodological approach used in this meta-analysis, please see the manuscript [Meta-analysis of Proportions Using Generalized Linear
Mixed Models](https://doi.org/10.1097/EDE.0000000000001232) by Lin & Chu (2020) and the sample code they provide.

## Requirements

All code is written in the programming language [R](https://www.r-project.org/). It is mostly easily run using the IDE [RStudio](https://rstudio.com/). An .Rproj file is included with this repository for your convenience.

The R packages required to reproduce the tables and figures at listed at the top of their respective scripts. They must be installed using `install.packages` prior to running the script.

## Reproductive the tables

Run `tables.R`.

## Reproducing the figures

Run `figures.R`.
