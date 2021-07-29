#' Estimate F statistics for many vertices
#'
#' @param stats Model fit as returned by \code{^lme_mass_fit_Rgw}
#' @param C Contrast vector
#'
#' @return
#' The function returns a list with entries F, pval, sgn and df for each vertex.
#'
#' @export
#'
#' @examples
#' \dontrun{C <- matrix(c(0, 1, 0, 0, 0, 0), nrow=1)}
#' \dontrun{FitRgw <- lme_mass_fit_Rgw(...)}
#' \dontrun{F_C <- lme_mass_F(FitRgw$stats, C)}

lme_mass_F<-function(stats, C) {

    nv<-length(stats)

    Fval<-NULL
    pval<-NULL
    sgn<-NULL
    df<-matrix(nrow=2,ncol=nv)

    fstats<-NULL

    for (i in c(1:nv)) {

        if (!(any(is.na(stats[[i]]$Bhat)))) {
            Bhat = stats[[i]]$Bhat
            Zcols = stats[[i]]$Zcols
            q = length(Zcols)
            nth = q*(q+1)/2+1
            invEI = stats[[i]]$invEI
            Pth = stats[[i]]$Pth
            Qthth = stats[[i]]$Qthth

            # check matrices: C, CBhat, OM

            # Estimation of the bias in the covariance matrix and computation of the
            # F-statistic and the degrees of freedom of the test.

            Bias = 0
            CBhat = stats[[i]]$CovBhat

            OM = t(C)%*%solve((C%*%CBhat%*%t(C)))%*%C # check this
            A1 = 0
            A2 = 0
            Term1 = OM%*%CBhat
            Term2 = CBhat%*%OM%*%CBhat

            for (k in c(1:nth))
            {
                Pk = drop(Pth[k,,])
                for (j in c(1:nth))
                {
                    Qkj = drop(Qthth[k,j,,])
                    Pj = drop(Pth[j,,])
                    Bias = Bias + invEI[k,j]*(Qkj-Pk%*%CBhat%*%Pj)
                    A1 = A1+invEI[k,j]*sum(diag(Term1%*%Pk%*%CBhat))*sum(diag(Term1%*%Pj%*%CBhat))
                    A2 = A2+invEI[k,j]*sum(diag(Term1%*%Pk%*%Term2%*%Pj%*%CBhat))
                }
            }

            szC = nrow(C)
            Bdf = (A1+6*A2)/(2*szC)
            g = ((szC+1)*A1-(szC+4)*A2)/((szC+2)*A2)
            d = 3*szC+2*(1-g)
            c1 = g/d
            c2 = (szC-g)/d
            c3 = (szC+2-g)/d
            EF = (1-A2/szC)^-1
            VF = (2/szC)*(1+c1*Bdf)/((1-c2*Bdf)^2*(1-c3*Bdf))
            ro = VF/(2*EF^2)
            m = 4+(szC+2)/(szC*ro-1)
            l = m/(EF*(m-2))
            Bias = CBhat%*%Bias%*%CBhat
            CovBhat = CBhat + 2*Bias
            Fstat = l%*%t(Bhat)%*%t(C)%*%solve((C%*%CovBhat%*%t(C)))%*%C%*%Bhat/szC # check this

            if (Fstat<0) {Fstat = 0}

            Fval[i] = Fstat
            pval[i] = max(1-stats::pf(Fstat,szC,m), 1e-30)
            contrast = C%*%Bhat
            sgn[i] = sign(contrast[1])
            df[,i] = matrix(c(szC,m),2)
        }
        else {

            Fval[i] = NA
            pval[i] = NA
            contrast = NA
            sgn[i] = NA
            df[,i] = NA
        }
    }

    fstats$F = Fval
    fstats$pval = pval
    fstats$sgn = sgn
    fstats$df = df

    return(fstats)

}
