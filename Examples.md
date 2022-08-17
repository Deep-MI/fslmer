
<!-- Examples.md is generated from Examples.Rmd. Please edit that file -->

# Examples

<!-- badges: start -->

<!-- badges: end -->

Here, we illustrate some application examples for the `fslmer` tools.

## Predicting individual values from model estimates

A common task after fitting a model is to make predictions given the
original data, or even given hypothetical data. Here we describe how
this can be done. We assume that an analysis as detailed on the main
page has been run already, at least up to the model-fitting stage,
including the estimation / extraction of random effects.

For a random-intercept, random-slope model where the intercept is in
column 1 of the design matrix `X` (denoted as \(X_0\), not used in the
formula below since it’s just a vector of \(1\)s) and the slope variable
(often: time) is in column 2 of the design matrix (denoted as \(X_1\) in
the formula below), the prediction would be computed as
follows:\`

\(\hat{Y} = \beta_0 + X_1 * \beta_1 + ... + X_n * \beta_n + b_0 + X_1 * b_1\)

where \(\beta_0\) … \(\beta_n\) are the population-level coefficients,
and \(b_0\) and \(b_1\) are the random-intercept and random-slope
coefficients, respectively.

The population-level betas in `rfx$Bhat[,i]` are the same as those in
`FitRgw$stats[[i]]$Bhat[1:2]`, where `i` indicates any particular
vertex. In our scenario, they are the model intercept and the
coefficient for the time/age variable.

To predict individual values, the subject-level betas in `rfx$Rfx` need
to be added / subtracted from the population-level betas, and we also
need to take into account the other covariates and their interactions.
All of this is done in the above formula.

In the matrix of random-effects estimates (`rfx$Rfx`), we first
replicate each individual random-effects coefficient as often as there
are observations per individual (based on the `ni` vector), and then
rearrange this matrix to a more convenient format. This gives a
multidimensional array of size `(number of observations, number of
random effects, number of vertices)`, while the original layout of
`rfx$Rfx` was `(number of cases, number of ranom effects * number of
vertices)`. We call the resulting array `gamma3D`, and use it to
introduce the subject-specific deviations from the population mean into
the prediction. In the following example, we do the prediction at the
vertex with the index `1234`.

``` r
# replicate cases to match observations
gamma2D <- apply(rfx$Rfx, MARGIN=2, FUN=rep, times=ni)

# reshape to 3D array
gamma3D <- array(gamma2D, dim=c(nrow(X), length(Zcols), ncol(Y)))

# make predictions
vertexIdx <- 1234
Y_pred <- X %*% FitRgw$stats[[vertexIdx]]$Bhat + gamma3D[, 1, vertexIdx] + X[,2] * gamma3D[, 2, vertexIdx]
```

The predicted and original thickness values should correlate fairly well
(depending on the overall model fit).

``` r
plot(Y[,vertexIdx], Y_pred)
```

We can also compare the original and predicted thickness maps using
FreeView. First, we need to do the prediction for all vertices, and then
create and save the corresponding images. Here, we do this for one
specific observation in one particular case. Again, the resulting maps
should look somewhat similar, but not totally identical.

``` r
# do the prediction
Y_pred <- matrix(0, nrow(Y), ncol(Y))

for (i in c(1:ncol(Y_pred))) {
    print(i)
    if ((!is.null(FitRgw$stats[[i]]$Bhat)) && (!any(is.na(FitRgw$stats[[i]]$Bhat)))) {
        Y_pred[,i] <- X %*% FitRgw$stats[[i]]$Bhat + gamma3D[, 1, i] + X[,2] * gamma3D[, 2, i]
    }
}

# save case 1, observation 1 (=first row of the design matrix) for original thickness data
vol <- NULL
vol$ndim1 <- ncol(Y)
vol$ndim2 <- 1
vol$ndim3 <- 1
vol$nframes <- 1
vol$x <- array(data=Y[1,], dim=c(ncol(Y), 1, 1, 1))
lme_savemgh(vol=vol, fname="Y_1_1b.mgh")
rm(vol)

# save case 1, observation 1 (=first row of the design matrix) for predicted thickness data
vol <- NULL
vol$ndim1 <- ncol(Y_pred)
vol$ndim2 <- 1
vol$ndim3 <- 1
vol$nframes <- 1
vol$x <- array(data=Y_pred[1,], dim=c(ncol(Y_pred), 1, 1, 1))
lme_savemgh(vol=vol, fname="Y_pred_1_1b.mgh")
rm(vol)
```

We have so far predicted the thickness values using existing data, but
it’s also possible to do this for arbitrary data. To do so, we need to
use a custom vector of predictors instead of the design matrix `X`.
Unless we want to model any deviations from the population estimate, I
guess we can neglect the random-effects estimate - because these
represent just deviations from the population estimate (e.g., population
intercept and slope).

Assuming that our design matrix `X` has 6 colums (intercept, time,
covariate1, covariate2, time \(\times\) covariate1, time \(\times\)
covariate 2), we now predict the thickness for a hypothetical case at
time 10, with values 3 and 4 for covariate 1 and covariate 2,
respectively, by creating a new design matrix `X_new`. It has just one
row, since we are predicting for one observation of one case. Note that
we also include the model intercept with a value of 1. Again, we do the
prediction for one vertex (id: 1234).

``` r
X_new <- matrix(c(1, 10, 3, 4, 30, 40), nrow=1, ncol=6)

Y_pred_new <- X_new %*% FitRgw$stats[[vertexIdx]]$Bhat
```

## How to compute prediction intervals

Afer computing the predictions of individual values as described above,
we still need the variability estimates for the predictions. A
complicating thing for mixed models is that for this purpose we need to
take both the between-subjects variability and the within-subjects
residual variance into account. This is typically done via
bootstrapping.

There are some variants how to do the bootstrapping, primarily
parametric and non-parametric, see <https://doi.org/10.1002/pst.1561>
for an an overview. Below is an example how to do a parametric
bootstrap. The basic idea is that we use the estimated parameters of
both the random-effects and the residual distributions to draw random
samples from these distributions. This is in contrast to the
non-parametric bootstrep, where draw random samples directly from the
observed random effects and residuals. In both cases, however, we obtain
variability estimates for . In the following example, these are the beta
weights and the predicted values, but other parameters of interest could
be added as well.

We assume that the `fslmer` package has been loaded and that a
vertex-wise analysis has been conducted at least until the model-fitting
stage. We also assume that the variables as in the above example are
present, in particular `X_new` and `vertexIdx`. Running the analysis
will take a while; we recommend to set the `n.boot` variable, which
determines how many samples are drawn, to a lower value initially to
test if everything works and to get an estimate of the expected run
time.

``` r

# load additional libraries

library(MASS)
library(Matrix)

# determine individual random-effects design matrices and create large block-diagonal matrix

zi <- list()
for (i in c(1:length(ni))) {
  zi[[i]] <- matrix(X[(sum(ni[0:(i-1)])+1):sum(ni[0:i]), Zcols], nrow=ni[i], ncol=length(Zcols))
}
Z <- bdiag(zi)  

# get variance-covariance matrix of the random efects

vcRE <- FitRgw$stats[[vertexIdx]]$Dhat

# compute population-level estimates

pop_resp <- X %*% FitRgw$stats[[vertexIdx]]$Bhat

# initialization and settings for the bootstrap

beta <- NULL
y_pred <- NULL

set.seed(1234)
n.boot <- 10000

# do the bootstrap

for (i in c(1:n.boot)) {
  
  # draw from random effects distribution
  gamma <- mvrnorm(n=length(ni), mu=rep(0, length(diag(vcRE))), Sigma=vcRE)
  
  # draw from residuals distribution
  e <- NULL
  for (j in c(1:length(ni))) {
    e <- c(e, mvrnorm(n=1, mu=rep(0, ni[j]), Sigma=diag(FitRgw$stats[[vertexIdx]]$phisqhat, ni[j], ni[j])))
  }
  
  # compute responses for current residual and random effects values
  boot_resp <- pop_resp + Z %*% matrix(t(gamma), ncol=1) + e
  
  # fit model for new responses
  boot_mod <- lme_fit_FS(X, Zcols, as.matrix(boot_resp), ni)
  
  # make prediction for new data using new coefs
  boot_pred <- X_new %*% boot_mod$Bhat
  
  # store betas and predictions
  beta <- rbind(beta, matrix(boot_mod$Bhat, nrow=1))
  y_pred <- rbind(y_pred, matrix(boot_pred, nrow=1))
  
  # clean up
  rm(gamma, e, boot_resp, boot_pred)
  
}

# summare bootstrap results

apply(beta, 2, sd)

apply(y_pred, 2, sd)
```
