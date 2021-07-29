
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fslmer

<!-- badges: start -->

<!-- badges: end -->

## Overview

This project is a port of Freesurfer’s Matlab-based LME tools to the R
programming language.

Please refer to the original documentation at
<https://surfer.nmr.mgh.harvard.edu/fswiki/LongitudinalStatistics> and
<https://surfer.nmr.mgh.harvard.edu/fswiki/LinearMixedEffectsModels> for
an overview and futher information about the software. For the original
code, see the original repository at
<https://github.com/NeuroStats/lme>.

If you use these tools in your analysis please cite:

  - Bernal-Rusiel J.L., Greve D.N., Reuter M., Fischl B., Sabuncu M.R.,
    2012. Statistical Analysis of Longitudinal Neuroimage Data with
    Linear Mixed Effects Models, NeuroImage 66, 249-260, 2012,
    <https://dx.doi.org/10.1016%2Fj.neuroimage.2012.10.065>

  - Bernal-Rusiel J.L., Greve D.N., Reuter M., Fischl B., Sabuncu M.R.,
    2013. Spatiotemporal Linear Mixed Effects Modeling for the
    Mass-univariate Analysis of Longitudinal Neuroimage Data, NeuroImage
    81, 358–370, 2013,
    <https://doi.org/10.1016/j.neuroimage.2013.05.049>

  - Reuter M., Schmansky N.J., Rosas H.D., Fischl B, 2012.Within-Subject
    Template Estimation for Unbiased Longitudinal Image Analysis,
    NeuroImage 61, 1402-1418, 2012,
    <http://dx.doi.org/10.1016/j.neuroimage.2012.02.084>

## Installation

You can install fslmer from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Deep-MI/fslmer")
```

## Running the program

There are two types of analyses that can be done: univariate and
mass-univariate. The univariate analysis is, in principle, similar to
classical mixed-effects analyses as implemented in the `lme4` or `nlme`
packages. The mass-univariate analysis is specifically tailored for the
surface-based vertex data from a preceding run of FreeSurfer’s
longitudinal analysis pipeline.

### Univariate analyses

  - Loading the data

Usually you should already have a longitudinal Qdec table, which
contains the subject IDs, image/observation/session IDs, time, and
demographics (see
<https://surfer.nmr.mgh.harvard.edu/fswiki/LongitudinalStatistics> for a
description of the format). Using Freesurfer’s `asegstats2table` and/or
`aparcstats2table` tools, you can extract the volume and/or thickness
estimates from the longitudinally processed data in your subjects
directory:

    asegstats2table --qdec-long PATH_TO_QDEC_TABLE/qdec.table.dat -t PATH_TO_ANALYSIS_DIRECTORY/aseg.long.table

  - Preparing the data

Your longitudinal Qdec table needs to contain the “fsid-base” (subject
ID) columns and “fsid” (image/observation/session ID), and should
contain an additional column indicating time from baseline or the
sequence of observations. Additional columns such as diagnoses,
demographics, or covariates may also be present. For the current
example, we assume that the timing info is stored in a numerical
variable called ‘time’, indicating years from baseline, and that a
categorical variable ‘DX’ indicates a diagnosis at baseline. You can
used other names and/or additional variables in your analysis, but the
example code has to be adapted accordingly.

After loading the data, some further preparations are needed, including
merging and sorting the data. In particular, the fslmer tools require
the data ordered according to time for each individual (that is, your
design matrix needs to have all the repeated assessments for the first
subject, then all for the second and so on).

``` r
# Load aseg/aparc and qdec tables into R:
aseg <- read.table("PATH_TO_DATA/aseg.long.table", header=True)
qdec <- read.table("PATH_TO_QDEC_TABLE/qdec.table.dat", header=True )

# Read the qdec and aseg tables; note that R replaces '-' (and other characters) 
# in variable names with '.'.
qdec <- read.csv("adni-180-long-fs53-new.qdec")
aseg <- read.csv("adni-180-long-fs53-new-aseg.csv")

# For the aseg table, split the longitudinal ID into fsid and fsid.base.
aseg$fsid <- sub("\\.long\\..*", "", aseg$Measure.volume)
aseg$fsid.base <- sub(".*\\.long\\.", "", aseg$Measure.volume)

# Merge qdec and aseg tables based on fsid and fsid.base variables; create new 
# table 'dat'.
dat <- merge(qdec, aseg, by=c("fsid.base", "fsid"))

# Sort the by 'fsid.base' and then by 'time'. 
dat <- dat[order(dat$fsid.base, dat$Time.From.Baseline), ]

# Extract the structure of interest (here: volume of left hippocampus; can be 
# changed) and store it in a new variable 'Y'. The variable needs to be column 
# vector (hence the 'matrix' command).
Y <- matrix(dat$Left.Hippocampus, ncol=1)

# As an auxiliary variable, create a vector of number of observations per each 
# subject, i.e. count the number of occurances for each fsid.base. 
ni <- matrix(unname(table(dat$fsid.base)), ncol=1)
```

  - Creating the design matrix and contrasts

Once you have your ordered data, you need to build your design matrix.
As an example, a simple linear model containing a group by time
interaction can be obtained with the following design matrix:

``` r
X <- model.matrix(~time*DX, dat)
```

Let us assume that the categorical variable ‘DX’ has two levels,
patients and controls. Then the model matrix ‘X’ will have four columns:
one for the intercept, a time variable, a binary group indicator created
from ‘DX’, and an interaction term between ‘DX’ and ‘time’. The contrast
`[0 0 0 1]`, which is a row vector with four elements, can then be used
to test the interaction between `DX` and `time`, which indicates
diverging slopes of thickness changes across time for the two groups.
Instead of a contrast vector, it also possible to specify contrast
matrices that have more than one row. Note that the number of columns of
the contrast vector or matrix always needs to correspond to the number
of columns in the model matrix.

``` r
C <- matrix(c(0, 0, 0, 1), nrow=1)
```

  - Estimating the model and conducting statistical inference

Estimate the model by using the `lme_fit_FS` function: The arguments
`X`, `Y`, and `ni` have been defined before, and `Zcols` is a vector
that indicates which terms of the model / columns of the design matrix
should be regarded as random effects: assuming that the Intercept is the
first column and time is the second column in `X`, use `Zcols=1` for a
random-intercept model, and `Zcols=c(1, 2)` for a
random-intercept-and-slope model.

``` r
stats <- lme_fit_FS(X, Zcols, Y, ni)
```

The function returns a list representing a model fit, with entries Bhat,
CovBhat and bihat, among others; Bhat contains the beta values, CovBhat
is the error covariance matrix, and bihat the random-effects
coefficients.

Finally, conduct an F-test using the `lme_F` function, which takes the
model fit and the contrast vector / matrix as inputs:

``` r
F_C <- lme_F(stats, C)
```

This will return a list with entries F, pval, sgn, and df: F is the
F-value, pval is the p-value, sgn is the sign of the beta coefficient,
and df are the degrees of freedom.

### Mass-univariate analyses

The mass-univariate analysis is run per hemisphere (lh and rh). Here and
in the following, we only describe the analysis for the left hemisphere,
which needs to done also for the right hemisphere.

Start with FreeSurfer’s `mris_preproc` comamand, which uses the
longitudinal qdec table to automatically find the longitudinally
processed data and assembles it into a single lh.thickness.mgh file.
Note that it is possible to use a different study template than the
standard fsaverage template, and that other measures than thickness can
be used as well (see the help for mris\_preproc).

    mris_preproc --qdec-long PATH_TO_QDEC_TABLE/qdec.table.dat --target fsaverage --hemi lh --meas thickness --out lh.thickness.mgh

The next step is to smooth the data; here we use a 10 mm FWHM kernel.

    mri_surf2surf --hemi lh --s fsaverage --sval lh.thickness.mgh --tval lh.thickness_sm10.mgh --fwhm-trg 10 --cortex --noreshape

  - Loading the data

The resulting file as well as will as a set of associated files will be
read in to R

``` r

# read thickness file
lh.thickness<-lme_openmgh("/PATH/TO/lh.thickness_sm10.mgh")

# read template surface
lh.sphere<-lme_readsurf("/PATH/TO/FREESURFER/DIRECTORY/subjects/fsaverage/surf/lh.sphere")

# read cortical label file
lh.cortex<-lme_readlabel("/PATH/TO/FREESURFER/DIRECTORY/subjects/fsaverage/label/lh.cortex.label")[,1]
```

  - Preparing the data

The design matrix is prepared in a similar way as for the univariate
analysis. We just repeat the code here with minimal comments, see above
for an explanation.

``` r

# read qdec table
dat <- read.csv("PATH_TO_QDEC_TABLE/qdec.table.dat")

# order qdec table based on fsid.base and time
dat <- dat[order(dat$fsid.base, dat$time), ]

# create vector of number of observations per subject
ni <- matrix(unname(table(dat$fsid.base)), ncol=1)

# create matrix from thickness data
Y <- t(drop(lh.thickness$x))

# create vector of cortical labels
maskvtx <- sort(lh.cortex)+1
```

  - Creating the design matrix and contrasts

<!-- end list -->

``` r
# create model matrix
X <- model.matrix(~time*DX, dat)

# create contrasts
C1 <- matrix(c(0, 1, 0, 0), nrow=1)
C2 <- matrix(c(0, 0, 0, 1), nrow=1)
```

  - Estimating the model and conducting statistical inference

<!-- end list -->

``` r
# obtain initial estimates
outFitInit <- lme_mass_fit_init(X=X, Zcols=1, Y=Y, ni=ni, maskvtx=maskvtx, numcore=12)

# run algorithm to identify spatially homogeneous regions
outRgGrow <- lme_mass_RgGrow(lh.sphere, outFitInit$Re0, outFitInit$Theta0, maskvtx=maskvtx, nst=2, prc=95)

# fit model
outFitRgw <- lme_mass_fit_Rgw(X, Zcols, Y, ni, outFitInit$Theta0, outRgGrow$Regions, lh.sphere, prs=12)

# inference
F_C1 <- lme_mass_F(outFitRgw$stats, C1)

# save
vol <- NULL
vol$ndim1 <- length(F_C1$F)
vol$ndim2 <- 1
vol$ndim3 <- 1
vol$nframes <- 1
vol$x <- F_C1$F
lme_savemgh(vol=vol, fname=OUTPUT_FILE_NAME)
```
