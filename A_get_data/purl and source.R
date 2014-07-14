library(knitr)
source("C:/R/util.alr.r")

rm(list=ls())

#script I
purl(input="C:/ALR/Models/CTRflows3/get_flow_data/daily_seasonal/I_get_gages_flow_agg.Rmd",
     output="C:/ALR/Models/CTRflows3/get_flow_data/daily_seasonal/I_get_gages_flow_agg.R",
     documentation=0)
source("C:/ALR/Models/CTRflows3/get_flow_data/daily_seasonal/I_get_gages_flow_agg.R")
# file.remove("C:/ALR/Models/CTRflows3/get_flow_data/daily_seasonal/I_get_gages_flow_agg.R")


#script II
purl(input="C:/ALR/Models/CTRflows3/get_flow_data/daily_seasonal/II_plot_gages_dams.Rmd",
     output="C:/ALR/Models/CTRflows3/get_flow_data/daily_seasonal/II_plot_gages_dams.R",
     documentation=0)
source("C:/ALR/Models/CTRflows3/get_flow_data/daily_seasonal/II_plot_gages_dams.R")
# file.remove("C:/ALR/Models/CTRflows3/get_flow_data/daily_seasonal/II_plot_gages_dams.R")

#script III
purl(input="C:/ALR/Models/CTRflows3/get_flow_data/daily_seasonal/III_get_met.Rmd",
     output="C:/ALR/Models/CTRflows3/get_flow_data/daily_seasonal/III_get_met.R",
     documentation=0)
source("C:/ALR/Models/CTRflows3/get_flow_data/daily_seasonal/III_get_met.R")
# file.remove("C:/ALR/Models/CTRflows3/get_flow_data/daily_seasonal/III_get_met.R")

# #script IV
# purl(input="C:/ALR/Models/CTRflows3/get_flow_data/daily_seasonal/IV_standardize.R",
#      output="C:/ALR/Models/CTRflows3/get_flow_data/daily_seasonal/IV_standardize.R",
#      documentation=0)
# source("C:/ALR/Models/CTRflows3/get_flow_data/daily_seasonal/IV_standardize.R")
# file.remove("C:/ALR/Models/CTRflows3/get_flow_data/daily_seasonal/IV_standardize.R")

#script V
purl(input="C:/ALR/Models/CTRflows3/get_flow_data/daily_seasonal/V_general_functions.Rmd",
     output="C:/ALR/Models/CTRflows3/get_flow_data/daily_seasonal/V_general_functions.R",
     documentation=0)
source("C:/ALR/Models/CTRflows3/get_flow_data/daily_seasonal/V_general_functions.R")
# file.remove("C:/ALR/Models/CTRflows3/get_flow_data/daily_seasonal/V_general_functions.R")



rm(qdaily,cdaily,date.template,cmatrices,catchments)
gc()
ls.objects()