% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lme_fit_FS.R
\name{lme_fit_FS}
\alias{lme_fit_FS}
\title{Fit a univariate linear mixed-effects model using Fisher scoring}
\usage{
lme_fit_FS(X, Zcols, y, ni, e = 10^-3)
}
\arguments{
\item{X}{Design matrix}

\item{Zcols}{Vector of random effect indices in design matrix}

\item{y}{Outcome variable}

\item{ni}{Vector indicating the repeated observations of each subject}

\item{e}{Tolerance (default: 10^-3)}
}
\value{
The function returns stats, a list of statistics, with the following entries:
Bhat, CovBhat, bihat, Covbihat, phisqhat, SIGMA, W, Dhat, X, Zcols, re,
ni, and lreml.
}
\description{
Fit a univariate linear mixed-effects model using Fisher scoring
}
\examples{
\dontrun{stats <- lme_fit_FS(X, Zcols, y, ni)}
}
