# Seasonal flow timeseries model<br> from paired flow and meterological/weather record
## 
## Exploratory Data Analysis: comparison of precipitation vs flow quantiles
## ALR July 2014<br>Conte Anadromous Fish Research Center
##   

### Load libraries, functions, etc


```r
library(lme4)
```

```
## Loading required package: Matrix
## Loading required package: Rcpp
```

```r
library(plyr)
library(ggplot2)
library(reshape2)
library(sp)
library(devtools)
```

```
## WARNING: Rtools is required to build R packages, but is not currently installed.
## 
## Please download and install Rtools 3.1 from http://cran.r-project.org/bin/windows/Rtools/ and then run find_rtools().
## 
## Attaching package: 'devtools'
## 
## The following objects are masked from 'package:utils':
## 
##     ?, help
## 
## The following object is masked from 'package:base':
## 
##     system.file
```

```r
library(xtable)


model_dir<-"C:/ALR/Models/MetToFlow"
model_data_dir<-"C:/ALR/Models_processed_data"


#utility functions, load (source) from my saved gist
source_gist("https://gist.github.com/anarosner/ba285306fc0ce9d812a5", sha1="b25a1b73e02cc2b2d2c590f6c0b2c9c9945fa980")
```

```
## Sourcing https://gist.githubusercontent.com/anarosner/ba285306fc0ce9d812a5/raw/48b3efa59d36c7dceacf9c1c0a6fd77ca20bfdb5/util.r
```



```r
setwd(file.path(model_data_dir,"flow_timeseries"))

load("gages_char_spatial.Rdata")
long_site_no<-gages.char.spatial$site_no[gages.char.spatial$qannual>=25 & gages.char.spatial$TNC_DamCount==0 & gages.char.spatial$OnChannelWaterSqKM<0.5]
length(long_site_no)
```

```
## [1] 7
```


```r
# View(t(gages.char.spatial@data[gages.char.spatial$site_no %in% long_site_no,]))
g<-t(gages.char.spatial@data[gages.char.spatial$site_no %in% long_site_no,])
dimnames(g)[[2]]<-g["site_no",]
print(xtable(g),type="html")
```

<!-- html table generated in R 3.0.3 by xtable 1.7-3 package -->
<!-- Mon Jul 14 16:29:47 2014 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> 01064400 </TH> <TH> 01073000 </TH> <TH> 01075800 </TH> <TH> 01105585 </TH> <TH> 01105600 </TH> <TH> 01118300 </TH> <TH> 01434025 </TH>  </TR>
  <TR> <TD align="right"> FEATUREID </TD> <TD> 9311777 </TD> <TD> 5845058 </TD> <TD> 6730101 </TD> <TD> 5863947 </TD> <TD> 5863773 </TD> <TD> 6140826 </TD> <TD> 4147946 </TD> </TR>
  <TR> <TD align="right"> agency_cd </TD> <TD> USGS </TD> <TD> USGS </TD> <TD> USGS </TD> <TD> USGS </TD> <TD> USGS </TD> <TD> USGS </TD> <TD> USGS </TD> </TR>
  <TR> <TD align="right"> site_no </TD> <TD> 01064400 </TD> <TD> 01073000 </TD> <TD> 01075800 </TD> <TD> 01105585 </TD> <TD> 01105600 </TD> <TD> 01118300 </TD> <TD> 01434025 </TD> </TR>
  <TR> <TD align="right"> station_nm </TD> <TD> LUCY BROOK NEAR NORTH CONWAY, NH </TD> <TD> OYSTER RIVER NEAR DURHAM, NH </TD> <TD> STEVENS BROOK NEAR WENTWORTH, NH </TD> <TD> TOWN BROOK AT QUINCY, MA </TD> <TD> OLD SWAMP RIVER NEAR SOUTH WEYMOUTH, MA </TD> <TD> PENDLETON HILL BROOK NEAR CLARKS FALLS, CT. </TD> <TD> BISCUIT BK ABOVE PIGEON BK AT FROST VALLEY NY </TD> </TR>
  <TR> <TD align="right"> dec_lat_va </TD> <TD> 44.07 </TD> <TD> 43.15 </TD> <TD> 43.84 </TD> <TD> 42.25 </TD> <TD> 42.19 </TD> <TD> 41.47 </TD> <TD> 42.00 </TD> </TR>
  <TR> <TD align="right"> dec_long_va </TD> <TD> -71.17 </TD> <TD> -70.97 </TD> <TD> -71.89 </TD> <TD> -71.00 </TD> <TD> -70.94 </TD> <TD> -71.83 </TD> <TD> -74.50 </TD> </TR>
  <TR> <TD align="right"> coord_acy_cd </TD> <TD> S </TD> <TD> S </TD> <TD> S </TD> <TD> S </TD> <TD> U </TD> <TD> H </TD> <TD> 1 </TD> </TR>
  <TR> <TD align="right"> dec_coord_datum_cd </TD> <TD> NAD83 </TD> <TD> NAD83 </TD> <TD> NAD83 </TD> <TD> NAD83 </TD> <TD> NAD83 </TD> <TD> NAD83 </TD> <TD> NAD83 </TD> </TR>
  <TR> <TD align="right"> huc_cd </TD> <TD> 01060002 </TD> <TD> 01060003 </TD> <TD> 01070001 </TD> <TD> 01090001 </TD> <TD> 01090001 </TD> <TD> 01090005 </TD> <TD> 02040104 </TD> </TR>
  <TR> <TD align="right"> drain_area_va </TD> <TD>  4.68 </TD> <TD> 12.10 </TD> <TD>  2.94 </TD> <TD>  4.11 </TD> <TD>  4.50 </TD> <TD>  4.02 </TD> <TD>  3.72 </TD> </TR>
  <TR> <TD align="right"> sv_begin_date </TD> <TD> 1987-10-09 </TD> <TD> 1936-03-19 </TD> <TD> 1987-10-08 </TD> <TD> 1982-12-16 </TD> <TD> 1967-04-18 </TD> <TD> 1958-07-11 </TD> <TD> 1983-02-16 </TD> </TR>
  <TR> <TD align="right"> sv_end_date </TD> <TD> 1992-05-27 </TD> <TD> 2014-06-02 </TD> <TD> 1998-10-01 </TD> <TD> 2014-06-02 </TD> <TD> 2014-06-05 </TD> <TD> 2014-06-25 </TD> <TD> 2014-05-20 </TD> </TR>
  <TR> <TD align="right"> sv_count_nu </TD> <TD>  42 </TD> <TD> 314 </TD> <TD>  98 </TD> <TD> 127 </TD> <TD> 268 </TD> <TD> 482 </TD> <TD> 413 </TD> </TR>
  <TR> <TD align="right"> da_sqkm </TD> <TD> 12.121 </TD> <TD> 31.339 </TD> <TD>  7.615 </TD> <TD> 10.645 </TD> <TD> 11.655 </TD> <TD> 10.412 </TD> <TD>  9.635 </TD> </TR>
  <TR> <TD align="right"> qseasonal </TD> <TD> 112 </TD> <TD> 249 </TD> <TD> 141 </TD> <TD> 105 </TD> <TD> 179 </TD> <TD> 209 </TD> <TD> 111 </TD> </TR>
  <TR> <TD align="right"> qannual </TD> <TD> 27 </TD> <TD> 63 </TD> <TD> 35 </TD> <TD> 27 </TD> <TD> 45 </TD> <TD> 53 </TD> <TD> 28 </TD> </TR>
  <TR> <TD align="right"> Forest </TD> <TD>  97.799 </TD> <TD>  73.528 </TD> <TD>  99.118 </TD> <TD>   8.895 </TD> <TD>  33.360 </TD> <TD>  83.289 </TD> <TD> 100.000 </TD> </TR>
  <TR> <TD align="right"> Herbacious </TD> <TD> 0.6441 </TD> <TD> 5.7393 </TD> <TD> 0.3083 </TD> <TD> 4.6497 </TD> <TD> 0.2096 </TD> <TD> 4.9337 </TD> <TD> 0.0000 </TD> </TR>
  <TR> <TD align="right"> Agriculture </TD> <TD> 0.62107 </TD> <TD> 5.88559 </TD> <TD> 0.33586 </TD> <TD> 0.06358 </TD> <TD> 0.08576 </TD> <TD> 6.15487 </TD> <TD> 0.00000 </TD> </TR>
  <TR> <TD align="right"> Developed </TD> <TD>  0.9354 </TD> <TD> 11.8234 </TD> <TD>  0.2382 </TD> <TD> 79.8569 </TD> <TD> 66.3444 </TD> <TD>  5.3472 </TD> <TD>  0.0000 </TD> </TR>
  <TR> <TD align="right"> DevelopedNotOpen </TD> <TD>  0.1150 </TD> <TD>  3.6475 </TD> <TD>  0.0000 </TD> <TD> 73.7455 </TD> <TD> 50.1991 </TD> <TD>  0.7673 </TD> <TD>  0.0000 </TD> </TR>
  <TR> <TD align="right"> Impervious </TD> <TD>  0.08174 </TD> <TD>  2.16175 </TD> <TD>  0.01724 </TD> <TD> 48.61790 </TD> <TD> 27.67411 </TD> <TD>  0.70356 </TD> <TD>  0.00000 </TD> </TR>
  <TR> <TD align="right"> AnnualTmaxC </TD> <TD> 11.77 </TD> <TD> 14.33 </TD> <TD> 11.39 </TD> <TD> 15.39 </TD> <TD> 15.47 </TD> <TD> 15.44 </TD> <TD> 10.40 </TD> </TR>
  <TR> <TD align="right"> AnnualTminC </TD> <TD> -0.6727 </TD> <TD>  2.4447 </TD> <TD> -0.1976 </TD> <TD>  5.2691 </TD> <TD>  4.8915 </TD> <TD>  4.5252 </TD> <TD>  0.1656 </TD> </TR>
  <TR> <TD align="right"> AnnualPrcpMM </TD> <TD> 1401 </TD> <TD> 1161 </TD> <TD> 1282 </TD> <TD> 1278 </TD> <TD> 1256 </TD> <TD> 1216 </TD> <TD> 1524 </TD> </TR>
  <TR> <TD align="right"> SummerPrcpMM </TD> <TD> 369.1 </TD> <TD> 292.1 </TD> <TD> 344.6 </TD> <TD> 299.2 </TD> <TD> 289.3 </TD> <TD> 292.3 </TD> <TD> 398.3 </TD> </TR>
  <TR> <TD align="right"> WinterPrcpMM </TD> <TD> 313.9 </TD> <TD> 247.2 </TD> <TD> 274.6 </TD> <TD> 314.1 </TD> <TD> 312.6 </TD> <TD> 283.0 </TD> <TD> 319.4 </TD> </TR>
  <TR> <TD align="right"> DrainageClass </TD> <TD> 1.593 </TD> <TD> 3.347 </TD> <TD> 3.105 </TD> <TD> 3.276 </TD> <TD> 3.919 </TD> <TD> 3.570 </TD> <TD> 2.488 </TD> </TR>
  <TR> <TD align="right"> HydrologicGroupAB </TD> <TD> 87.5000 </TD> <TD> 49.5472 </TD> <TD> 76.6733 </TD> <TD> 45.7568 </TD> <TD> 44.2069 </TD> <TD> 68.9397 </TD> <TD>  0.1498 </TD> </TR>
  <TR> <TD align="right"> HydrologicGroupCD </TD> <TD> 12.50 </TD> <TD> 47.43 </TD> <TD> 22.52 </TD> <TD> 54.24 </TD> <TD> 55.79 </TD> <TD> 31.06 </TD> <TD> 98.35 </TD> </TR>
  <TR> <TD align="right"> SurficialCoarseC </TD> <TD> 48.000 </TD> <TD>  7.842 </TD> <TD> 13.059 </TD> <TD> 15.538 </TD> <TD> 26.931 </TD> <TD>  4.692 </TD> <TD>  0.000 </TD> </TR>
  <TR> <TD align="right"> PercentSandy </TD> <TD>  4.1290 </TD> <TD>  0.2672 </TD> <TD>  9.1447 </TD> <TD> 11.0852 </TD> <TD> 17.1089 </TD> <TD>  0.4179 </TD> <TD>  0.0000 </TD> </TR>
  <TR> <TD align="right"> ReachElevationM </TD> <TD> 295.79000 </TD> <TD>  48.70422 </TD> <TD> 282.19044 </TD> <TD>   0.01083 </TD> <TD>  33.64731 </TD> <TD>  56.35676 </TD> <TD> 727.91603 </TD> </TR>
  <TR> <TD align="right"> BasinElevationM </TD> <TD> 431.60 </TD> <TD>  58.71 </TD> <TD> 486.11 </TD> <TD>  24.05 </TD> <TD>  42.92 </TD> <TD>  90.72 </TD> <TD> 865.67 </TD> </TR>
  <TR> <TD align="right"> ReachSlopePCNT </TD> <TD> 5.2180 </TD> <TD> 0.8770 </TD> <TD> 7.6185 </TD> <TD> 0.0113 </TD> <TD> 0.5419 </TD> <TD> 1.2008 </TD> <TD> 4.9676 </TD> </TR>
  <TR> <TD align="right"> BasinSlopePCNT </TD> <TD> 25.189 </TD> <TD>  4.996 </TD> <TD> 22.682 </TD> <TD>  4.987 </TD> <TD>  3.474 </TD> <TD>  6.912 </TD> <TD> 25.826 </TD> </TR>
  <TR> <TD align="right"> TotDASqKM </TD> <TD> 14.296 </TD> <TD> 32.043 </TD> <TD>  7.866 </TD> <TD> 15.985 </TD> <TD> 11.826 </TD> <TD> 13.701 </TD> <TD>  9.918 </TD> </TR>
  <TR> <TD align="right"> TNC_DamCount </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> </TR>
  <TR> <TD align="right"> deg_barr_1 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> </TR>
  <TR> <TD align="right"> deg_barr_2 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> </TR>
  <TR> <TD align="right"> deg_barr_3 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> </TR>
  <TR> <TD align="right"> deg_barr_4 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> </TR>
  <TR> <TD align="right"> deg_barr_6 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> </TR>
  <TR> <TD align="right"> deg_barr_7 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> </TR>
  <TR> <TD align="right"> OnChannelWaterSqKM </TD> <TD> 0.000000 </TD> <TD> 0.445174 </TD> <TD> 0.000000 </TD> <TD> 0.164880 </TD> <TD> 0.005342 </TD> <TD> 0.033938 </TD> <TD> 0.000000 </TD> </TR>
  <TR> <TD align="right"> OnChannelWetlandSqKM </TD> <TD> 0.1151 </TD> <TD> 2.3047 </TD> <TD> 0.0000 </TD> <TD> 0.1046 </TD> <TD> 0.7245 </TD> <TD> 0.7986 </TD> <TD> 0.0000 </TD> </TR>
  <TR> <TD align="right"> OffChannelWaterSqKM </TD> <TD> 0.000000 </TD> <TD> 0.048012 </TD> <TD> 0.000000 </TD> <TD> 0.036367 </TD> <TD> 0.007902 </TD> <TD> 0.049986 </TD> <TD> 0.000000 </TD> </TR>
  <TR> <TD align="right"> OffChannelWetlandSqKM </TD> <TD> 0.03532 </TD> <TD> 0.90935 </TD> <TD> 0.00272 </TD> <TD> 0.17506 </TD> <TD> 0.91673 </TD> <TD> 0.30900 </TD> <TD> 0.00000 </TD> </TR>
  <TR> <TD align="right"> large_barriers </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> </TR>
  <TR> <TD align="right"> small_barriers </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> <TD> 0 </TD> </TR>
   </TABLE>


```r
load("dseasonal.Rdata")
```

```
## Warning: cannot open compressed file 'dseasonal.Rdata', probable reason
## 'No such file or directory'
```

```
## Error: cannot open the connection
```

```r
dseasonal.long<-dseasonal[dseasonal$site_no %in% long_site_no,]
```

```
## Error: object 'dseasonal' not found
```

```r
dseasonal.long<-subset(dseasonal.long,subset=!is.na(precip.e) & !is.na(precip.e.lag1) & !is.na(precip.e.lag2) )
```

```
## Error: object 'dseasonal.long' not found
```

```r
# dseasonal.long[is.na(dseasonal.long$precip.e.lag2),]
# dim(dseasonal.long)



# 1 = Complete barrier to all fish (12+ feet)
# 2 = Small dam barrier (1-12 feet)
# 3 = Partial breach
# 4 = Barrier with fish ladder
# 5 = Unlikely barrier - fully breached, weir, under 1ft dam (also COND=NO DAM or COND=REM or COND=DEL)
# 6 = Unknown, assumed full barriers
# 7 = Locks
```


```r
basic<-lmer(log(flow) ~ 
              log(da_sqkm) +
               precip.e +  precip.e.lag1 + precip.e.lag2 +    #precip.e.lag3 + 
               (1|year)+(1|site_no),
          na.action=NULL, #to ensure we attached fitted values w/ correct site_no
         data=dseasonal.long) 
```

```
## Error: 'data' not found, and some variables missing from formula
## environment
```

```r
season.names<-c("winter","spring","summer","fall")
# basic.season<-list()
weights<-as.data.frame(matrix(nrow=4,ncol=3))

for (i in 1:length(season.names)) {
     temp<-update(basic, subset=season==season.names[i])
     weights[i,]<-fixef(temp)[c("precip.e", "precip.e.lag1", "precip.e.lag2")]/
          sum(fixef(temp)[c("precip.e", "precip.e.lag1", "precip.e.lag2")])
#      basic.season[season.names[i]]<-temp
}
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'update': Error: object 'basic' not found
```

```r
row.names(weights)<-season.names
names(weights)<-c("precip.e", "precip.e.lag1", "precip.e.lag2")

weights
```

```
##        precip.e precip.e.lag1 precip.e.lag2
## winter       NA            NA            NA
## spring       NA            NA            NA
## summer       NA            NA            NA
## fall         NA            NA            NA
```

### compare precip vs flow quantiles from sites w long records

```r
# 4 sites, 4 seasons

rm(quant.melt)
for (i in long_site_no) {
#      i<-long_site_no[2]
     print(i)
     for (j in 1:length(season.names)) {
          print(j)
          print(season.names[j])
          temp<-dseasonal.long[dseasonal.long$site_no == i &
                                    dseasonal.long$season==season.names[j],
                            c("precip.e","precip.e.lag1","precip.e.lag2","flow","year")]
          temp.weighted<-temp[,1]*weights[j,1] + temp[,2]*weights[j,2] + temp[,3]*weights[j,3]
          
          cdf.q<-ecdf(temp$flow)
          cdf.p<-ecdf(temp$precip.e)
          cdf.p.weighted<-ecdf(temp.weighted)
          
          temp.q<-data.frame(year=temp$year,
                             flow=cdf.q(temp$flow),
                             precip=cdf.p(temp$precip.e),
                             precip.weighted=cdf.p.weighted(temp.weighted),
                             site=i,
                             season=season.names[j],
                             stringsAsFactors=F)

          temp.q$high<-apply(temp.q[,c("flow","precip")],MARGIN=1,max)
          temp.q$low<-apply(temp.q[,c("flow","precip")],MARGIN=1,min)
          temp.q$high.weighted<-apply(temp.q[,c("flow","precip.weighted")],MARGIN=1,max)
          temp.q$low.weighted<-apply(temp.q[,c("flow","precip.weighted")],MARGIN=1,min)
#           temp.q$under<-(temp.q$flow-temp.q$precip)>0
          head(temp.q)
     
          temp.melt<-melt(temp.q,id.vars=c("year","site", "season","high","low","high.weighted","low.weighted"))
          {     
          if(!exists("quant.melt"))
               quant.melt<-temp.melt
          else{
#                print((nrow(quant.melt)+1))
#                print((nrow(quant.melt)+nrow(temp.melt)))
               quant.melt[(nrow(quant.melt)+1):
                               (nrow(quant.melt)+nrow(temp.melt)),]<-temp.melt
          }
          }
     }
}
```

```
## [1] "01064400"
## [1] 1
## [1] "winter"
```

```
## Error: object 'dseasonal.long' not found
```

```r
head(quant.melt)
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'head': Error: object 'quant.melt' not found
```

```r
quant.melt$variable<-(capitalize(as.character(quant.melt$variable)))
```

```
## Error: object 'quant.melt' not found
```

```r
quant.melt$season<-capitalize(quant.melt$season)
```

```
## Error: object 'quant.melt' not found
```

```r
# names(quant.melt)[names(quant.melt)=="value"]<-"quantile"

# quant.melt$year.factor<-as.factor(quant.melt$year)
unique(quant.melt$variable)
```

```
## Error: object 'quant.melt' not found
```

    

```r
gg.quant<-ggplot(data=subset(quant.melt,variable %in% c("Flow","Precip" )),
                 aes(x=year,y=value,color=variable,lty=variable)) +
     theme_bw()+scale_colour_hue(l=40)+ylab("quantile")+theme(legend.title=element_blank())
```

```
## Error: object 'quant.melt' not found
```

```r
limits <- aes(ymax = high, ymin=low)

#      i<-capitalize(season.names[1])
for (i in capitalize(season.names)) {  
     temp<-gg.quant +  
          geom_line(subset=.(season==i),colour="grey50") + 
          geom_errorbar(limits,colour="grey30",lwd=1.1,subset=.(season==i)) +
          geom_point(pch=15,cex=3,subset=.(season==i)) +
          facet_wrap(~site,nrow=2) + ggtitle(i)
     print(temp)
     
     print(
          temp %+% subset(quant.melt,variable %in% c("Flow","Precip.weighted")) +
              ggtitle(paste("Weighted quantiles:",i))
           )
     
#      setwd(file.path(model_dir,"B_eda"))
#      ggsave(temp,filename=paste0("quantile_",i,".png"),width=16,height=8)
#      ggsave(temp,filename=paste0("quantile_weighted_",i,".png"),width=16,height=8)
}
```

```
## Error: object 'gg.quant' not found
```




