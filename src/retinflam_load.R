#! /usr/bin/Rscript
#
# created on Wed Jun 14 12:03:36 2023
# Finn Rabe, <finn dot rabe at bli dot uzh dot ch>
#-----------------------------------------------------------------------
#
# common libraries
libs <- c(
    "tidyr",
    "dplyr",
    "devtools",
    "htmltools",
    "png",
    "grid",
    "gridExtra",
    "represearch",
    "ggplot2",
    "tidyverse",
    "cowplot",
    "knitr",
    "rmarkdown",
    "papaja",
    "kableExtra",
    "tidyr",
    "lme4",
    "robustlmm",
    "lmerTest",
    "lmtest",
    "sjPlot",
    "webshot2",
    "reticulate",
    "memisc",
    "MASS",
    "DiagrammeR",
    "DiagrammeRsvg",
    "magrittr",
    "rsvg",
    "sfsmisc",
    "FDRestimation",
    "sensemakr",
    "xtable",
    "magick",
    "car",
    "geomtextpath",
    "insight",
    "parameters",
    "psycho",
    "boot",
    "boot.pval",
    "Gmisc",
    "visreg",
    "rticles",
    "mice",
    "modelsummary"
)

if (!require("pacman")) install.packages("pacman")
library("pacman")
pacman::p_load(char = libs)
