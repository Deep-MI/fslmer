lme_mass_F<-function(stats,CM)
{
    if (is.element("Bhat",names(stats))) {stats<-list(stats)}
    
    if (is.element("C",names(CM))) {CM<-list(CM)}
    
    nv<-length(stats)
    
    pval<-NULL
    sgn<-NULL
    df<-matrix(nrow=2,ncol=nv)
    
    fstats<-NULL
    
    for (i in c(1:nv))
        {
        Bhat = stats[[i]]$Bhat
        C = CM[[i]]$C
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

        F[i] = Fstat
        pval[i] = max(1-pf(Fstat,szC,m), 1e-30)
        contrast = C%*%Bhat
        sgn[i] = sign(contrast[1])
        df[,i] = matrix(c(szC,m),2)

    }

    fstats$F = F
    fstats$pval = pval
    fstats$sgn = sgn
    fstats$df = df

    return(fstats)

}
