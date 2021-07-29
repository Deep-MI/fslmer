#' Mass-univariate linear mixed-effects model fitting
#'
#' @param X Design matrix
#' @param Zcols Vector of random effect indices in design matrix
#' @param Y Outcome variable
#' @param ni Vector indicating the repeated observations of each subject
#' @param maskvtx Indices (one-based) for vertices to include in the analysis (default: all)
#' @param prs (currently ignored)
#' @param e  Tolerance (default: 10^-1)
#' @param Xrows Internal parameter (default: NA)
#' @param numcore Number of cores for parallel computing (default: 1)
#'
#' @return
#' This function returns a list with the following parameters: 
#' Bhat, CovBhat, phisqhat, Dhat, X, Zcols, invEI, Pth, Qthth, lreml.

#' @export
#'
#' @examples
#' stats <- lme_mass_fit_vw(X, Zcols, Y, ni, numcore=1)
    
    
lme_mass_fit_vw<-function(X,Zcols,Y,ni,maskvtx=NA,prs=1,e=10^-1,Xrows=NA,numcore=1) {
    
    #
    
    stats<-NULL
    
    # check if parallel computing is feasible
    
    if (numcore==1) print("No parallel computing enabled (not recommended)",quote=F)
    
    #
    
    nv0<-ncol(Y)

    if (any(is.na(maskvtx)))
    {
        print(paste("Selecting all available vertices (",nv0,").",sep=""),quote=F)
        maskvtx<-c(1:nv0)
    }

    Y<-Y[,maskvtx,drop=F]
    
    p<-ncol(X)
    
    nv<-ncol(Y)
    
    rfcols<-Zcols
    Zcols<-matrix(0,1,p)
    Zcols[,rfcols]<-matrix(1,1,length(rfcols))
    
    if (!is.na(Xrows))
    {
        Xrows<-Xrows[,maskvtx]
    }
    
    stats1<-lme_mass_fit(X,NA,Xrows,Zcols,Y,ni,prs,e,numcore)
    
    for (i in c(1:nv))
    {
        stats[maskvtx[i]]<-stats1[i]
    }        
    
    return(stats)
}
