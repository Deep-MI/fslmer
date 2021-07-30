#' Estimate F statistics
#'
#' @param stats Output from the \code{lme_fit_FS} function
#' @param C Contrast vector
#'
#' @return
#' This function returns a list with entries F, pval, sgn, and df.
#'
#' @export
#'
#' @examples
#' \dontrun{stats <- lme_fit_FS(...)}
#' \dontrun{C <- matrix(c(0, 1, 0, 0), nrow=1)}
#' \dontrun{F <- lme_F(stats, C)}

lme_F<-function(stats, C) {

	# --------------------------------------------------------------------------

	X = stats$X

	if (ncol(C) != ncol(X))
	    {
	    stop(cat("The number of colums in C must be equal to the number of colums in the design matrix X"));
	    }

	p = ncol(X)
	ni = stats$ni
	Zcols = stats$Zcols
	Z = X[,Zcols,drop=F]
	q = length(Zcols)
	nth = (q * (q + 1) / 2) + 1
	W = stats$W
	V = stats$SIGMA
	CBhat = stats$CovBhat
	L = chol(stats$Dhat)
	phi = sqrt(stats$phisqhat)

	# Computation of Rijs
	Rthth = array(0,dim = c(nth,nth,p,p))
	jk = 0
	for (k in 1:q)
	    {
	    for (j in 1:k)
	        {
	        jk = jk + 1
	        uv = 0
	        for (v in 1:q)
	            {
	            for (u in 1:v)
	                {
	                uv = uv + 1
	                posi = 1;
	                SumR = 0
	                for (i in 1:length(ni))
	                    {
	                    posf = posi + ni[i] - 1
	                    Xi = X[posi:posf,,drop=F];
	                    Zi = Z[posi:posf,,drop=F];
	                    Wi = W[posi:posf,1:ni[i],drop=F]
	                    Ekj = matrix(0,q,q);
	                    Ekj[k,j] = 1;
	                    Euv = matrix(0,q,q);
	                    Euv[u,v] = 1
	                    Ai = Zi %*% Ekj %*% Euv %*% t(Zi)
	                    Ri = t(Xi) %*% Wi %*% (Ai+ t(Ai)) %*% Wi %*% Xi
	                    SumR = SumR + Ri
	                    posi = posf + 1
	                    }
	                Rthth[jk,uv,,] = SumR
	                }
	            }
	        }
	    }

	#Computation of Pis,Qijs and the expected information matrix EI.

	list2env(lme_EI(X,Zcols,W,CBhat,V,L,phi,ni),envir=environment())

	invEI = lme_mldivide(EI,diag(nth),ginv.do=F) # check this

	#Estimation of the bias in the covariance matrix and computation of the
	#F-statistic and the degrees of freedom of the test.

	Bias = 0
    OM = t(C) %*% solve(C %*% CBhat %*% t(C)) %*% C # check this
    A1 = 0;
    A2 = 0
    Term1 = OM %*% CBhat
    Term2 = CBhat %*% OM %*% CBhat
    for (k in 1:nth)
        {
    	Pk = drop(Pth[k,,])
    	for (j in 1:nth)
    	    {
            Qkj = drop(Qthth[k,j,,]);
            Rkj = drop(Rthth[k,j,,]);
            Pj = drop(Pth[j,,])
            Bias = Bias + invEI[k,j] * (Qkj - Pk %*% CBhat %*% Pj - 0.25 * Rkj)
            A1 = A1 + invEI[k,j] %*% sum(diag(Term1 %*% Pk %*% CBhat)) * sum(diag(Term1 %*% Pj %*% CBhat))
            A2 = A2 + invEI[k,j] %*% sum(diag(Term1 %*% Pk %*% Term2 %*% Pj %*% CBhat))
    	    }
        }
    szC = nrow(C)
    Bdf = (A1 + 6 %*% A2) / (2*szC)
    g = ((szC + 1) * A1-(szC + 4) * A2) / ((szC + 2) * A2)
    d = 3 * szC + 2 * (1-g)
    c1 = g / d
    c2 = (szC - g) / d
    c3 = (szC + 2 - g) / d
    EF = (1 - A2 / szC)^-1
    VF = (2 / szC) * (1 + c1 * Bdf) / ((1 - c2 * Bdf)^2 * (1 - c3 * Bdf))
    ro = VF / (2 * EF^2)
    m = 4 + (szC + 2) / (szC * ro -1)
    l = m / (EF * (m-2))
    Bias = CBhat %*% Bias %*% CBhat
    CovBhat = CBhat + 2 * Bias
    Bhat = stats$Bhat

    F = l * t(Bhat) %*% t(C) %*% solve(C %*% CovBhat %*% t(C)) %*% C %*% Bhat /szC # check this
    if (F<0)
        {
    	F = 0
        }

    pval = max(1-stats::pf(F,szC,m))
    sgn = sign(C %*% Bhat)
    df = c(szC,m)

	  F_C<- list("F" = F,"pval" = pval,"sgn" = sgn,"df" = df)
    return(F_C)

}
