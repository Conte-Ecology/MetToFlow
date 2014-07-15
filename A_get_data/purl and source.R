library(knitr)
source("C:/R/util.alr.r")

rm(list=ls())

setwd("c:/alr/models/mettoflow/a_get_data")

#script 0
purl(input="0_general_functions.Rmd",output="0_general_functions.R",documentation=0)
source("0_general_functions.R")
file.remove("0_general_functions.R")

