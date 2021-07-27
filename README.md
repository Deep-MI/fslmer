
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fslmer

<!-- badges: start -->

<!-- badges: end -->

The goal of fslmer is to …

## Installation

You can install the released version of fslmer from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("fslmer")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Deep-MI/fslmer")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(fslmer)
## basic example code
```

<!-- # lme-r -->

<!-- ## Overview -->

<!-- This project is a port of Freesurfer's Matlab-based LME tools to the R programming language. -->

<!-- Please refer to the original documentation at https://surfer.nmr.mgh.harvard.edu/fswiki/LongitudinalStatistics and https://surfer.nmr.mgh.harvard.edu/fswiki/LinearMixedEffectsModels for an overview and futher information about the software. For the original code, see the original repository at https://github.com/NeuroStats/lme. -->

<!-- If you use these tools in your analysis please cite: -->

<!-- - Bernal-Rusiel J.L., Greve D.N., Reuter M., Fischl B., Sabuncu M.R., 2012. Statistical Analysis of Longitudinal Neuroimage Data with Linear Mixed Effects Models, NeuroImage 66, 249-260, 2012, https://dx.doi.org/10.1016%2Fj.neuroimage.2012.10.065 -->

<!-- - Bernal-Rusiel J.L., Greve D.N., Reuter M., Fischl B., Sabuncu M.R., 2013. Spatiotemporal Linear Mixed Effects Modeling for the Mass-univariate Analysis of Longitudinal Neuroimage Data, NeuroImage 81, 358–370, 2013, https://doi.org/10.1016/j.neuroimage.2013.05.049 -->

<!-- - Reuter M., Schmansky N.J., Rosas H.D., Fischl B, 2012.Within-Subject Template Estimation for Unbiased Longitudinal Image Analysis, NeuroImage 61, 1402-1418, 2012, http://dx.doi.org/10.1016/j.neuroimage.2012.02.084 -->

<!-- ## Installing the program -->

<!-- TODO -->

<!-- ## Running the program -->

<!-- There are two types of analyses that can be done: univariate and mass-univariate.  -->

<!-- ### Univariate analyses -->

<!-- Usually you should already have a longitudinal Qdec table containing the ids, subject ids, time, and demographics. Using Freesurfer's `asegstats2table` and `aparcstats2table` tools, you can directly grab the stats values from the longitudinally processed data in your subjects directory: -->

<!-- ``` -->

<!-- asegstats2table --qdec-long qdec.table.dat -t PATH_TO_QDEC_TABLE/aseg.long.table  -->

<!-- ``` -->

<!-- Load these tables into R: -->

<!-- ```R -->

<!-- aseg <- read.table("PATH_TO_DATA/aseg.long.table", header=True) -->

<!-- qdec <- read.table("PATH_TO_QDEC_TABLE/qdec.table.dat", header=True ) -->

<!-- ``` -->

<!-- For computational efficiency reasons, these tools require the data ordered according to time for each individual (that is, your design matrix needs to have all the repeated assessments for the first subject, then all for the second and so on). Assume that your longitudinal Qdec table contains four columns: "fsid" (image ID), "fsid-base" (subject ID), "time" (time from baseline) and "group"(a binary variable) then you can use the following code: -->

<!-- ```R -->

<!-- # Create a set of auxiliary variables: -->

<!-- # qdec ... ordered assessment variables -->

<!-- # aseg ... ordered aseg table -->

<!-- # Y    ... ordered data  -->

<!-- # ni   ... a matrix with the number of repeated measures for each subject -->

<!-- # Load auxiliary libraries -->

<!-- library(plyr) -->

<!-- # Convert categorical strings into numeric values -->

<!-- qdec["group"] <- factor(qdec$group, labels=c("CN"="1", "AD"="0")) -->

<!-- # Add columns fsid.base and time to the aseg dataframe -->

<!-- aseg["fsid.base"] <- qdec$fsid.base -->

<!-- aseg["time"] <- qdec$time -->

<!-- # Sorts aseg and qdec dataframes first by fsid.base and then by time -->

<!-- aseg <-aseg[order(fsid.base,time),] -->

<!-- qdec <- qdec[order(fsid.base,time),] -->

<!-- # Extract the structure of interest from the aseg dataframe -->

<!-- Y <- aseg$Left.Hippocampus -->

<!-- # Count the number of occurances for each fsid, delete column `x`, and convert  -->

<!-- # to a matrix -->

<!-- ni <- count(qdec$fsid.base)          -->

<!-- ni$x <- NULL -->

<!-- ni <- as.matrix(ni) -->

<!-- ``` -->

<!-- Once you have your ordered data, you need to build your design matrix. As ab example, a simple linear model containing a group by time interaction can be obtained with the following design matrix: -->

<!-- ```R -->

<!-- X <- model.matrix(~time*group, qdec) -->

<!-- ``` -->

<!-- The contrast `[0 0 0 1]` can be used to test the interaction between `group` and `time`, which indicates diverging slopes of thickness changes across time for the two groups:  -->

<!-- ```R -->

<!-- C <- matrix(c(0, 0, 0, 1), 1) -->

<!-- ``` -->

<!-- Estimate the model by using this function: -->

<!-- ```R -->

<!-- stats <- lme_fit_FS(X, Zcols, y, ni) -->

<!-- ``` -->

<!-- Conduct an F-test: -->

<!-- ```R -->

<!-- F_C <- lme_F(stats, C) -->

<!-- ``` -->

<!-- The p-value of the test is in `F_C$pval`. -->

<!-- ## Mass-univariate analyses -->

<!-- # read data -->

<!-- Yin<-lme_openmgh("~/Work/lme/A_data/orig/mass_univariate/lh.50sMCI_vs_50cMCI_long_thickness_sm10.mgh") -->

<!-- Xin<-read.csv("~/Work/lme/A_data/eval/sMCI50_vs_cMCI50.csv") -->

<!-- nin<-read.csv("~/Work/lme/A_data/eval/sMCI50_vs_cMCI50_idx.csv") -->

<!-- lh.sphere<-lme_readsurf("~/Work/freesurfer/subjects/fsaverage/surf/lh.sphere") -->

<!-- lh.cortex<-lme_readlabel("~/Work/freesurfer/subjects/fsaverage/label/lh.cortex.label")[,1] -->

<!-- # prepare data -->

<!-- Y<-t(as.matrix(drop(Yin$x))) -->

<!-- X<-as.matrix(Xin) -->

<!-- Zcols<-c(1,2) -->

<!-- ni<-as.matrix(nin) -->

<!-- maskvtx<-sort(lh.cortex+1) -->

<!-- prs<-12 -->

<!-- maskvtx2<-maskvtx -->

<!-- maskvtx2<-maskvtx2[-which(colSums(Y[,maskvtx2])==0)] -->

<!-- # algorithm with region-growing -->

<!-- timeFitInit<-system.time(outFitInit<-lme_mass_fit_init(X,Zcols,Y,ni,maskvtx=maskvtx2,numcore=prs)) -->

<!-- #outFitEMInit<-lme_mass_fit_EMinit(X,Zcols,Y,ni,maskvtx=maskvtx,numcore=prs) -->

<!-- timeRgGrow<-system.time(outRgGrow<-lme_mass_RgGrow(lh.sphere,outFitInit$Re0,outFitInit$Theta0,maskvtx=maskvtx2,nst=2,prc=95)) -->

<!-- timeFitRgw<-system.time(outFitRgw<-lme_mass_fit_Rgw(X,Zcols,Y,ni,outFitInit$Theta0,outRgGrow$Regions,lh.sphere,prs=prs)) -->

<!-- save(file="~/Work/lme/D_analyses/D2_mass-univariate-rgw/run22/lme_mass_rgw_result_run21.RData",list=c("maskvtx","maskvtx2","timeFitInit","timeFitRgw","timeRgGrow","outFitInit","outRgGrow","outFitRgw","ni","prs","X","Xin","Zcols")) -->
