lme_RgFSfit<-function(X,Zcols,Y,ni,Theta0,Dist,model='exp',e=0.1)
{
    # -------------------------------------------------------------------------
    # check libraries
    
    st<-require(MASS)
    
    if (!st) stop("MASS package is required for this analysis, please install.")
    
    #
    
    # -------------------------------------------------------------------------
    # Auxiliary functions

    reshapeM<-function(Y,m,ni,nv)
    {
        
        posi<-1
        posiM<-1
        RM<-matrix(0,sum(ni)*nv,1)
        for (i in c(1:m))
        {
            posf<-posi+ni[i]-1
            posfM<-posiM+ni[i]*nv-1
            RM[c(posiM:posfM)]<-matrix(Y[posi:posf,,drop=F],ni[i]*nv,1)
            posi<-posf+1
            posiM<-posfM+1
        }
        
        return(RM)
        
    }
    
    sptDer<-function(Dist,a,b,model)
    {
        
        if (model=='gauss')
        {
            Distsq<-Dist*Dist
            G <- exp(-b*b*Dist-a*a*Distsq)
            tryCatch(
            {
                chol(G)
            },
            error=function(e)
            {
                outEig <- eigen(G)
                G <- outEig$vectors%*%diag(sapply(outEig$values,max,1e-5))%*%t(outEig$vectors)
            }
            )
            GDa<-(-2*a*Distsq*G)
            GDb<-(-2*b*Dist*G)
        }
        else
        {
            if (model=='exp')
            {
                G <- exp(-a*a*Dist)
                tryCatch(
                {
                    chol(G)
                },
                error=function(e)
                {
                    outEig <- eigen(G)
                    G <- outEig$vectors%*%diag(sapply(outEig$values,max,1e-5))%*%t(outEig$vectors)
                }
                )
                GDa<-(-2*a*Dist*G)
                GDb<-NA
            }
            else
            {
                if (model=='sph')
                {
                    C<-a*a*Dist
                    G<-C
                    Ind<-(C<=1)
                    C3<-C[Ind]^3
                    G[Ind]<-1-0.5*(3*C[Ind]-C3)
                    G[C>1]<-0
                    tryCatch(
                    {
                        chol(G)
                    },
                    error=function(e)
                    {
                        outEig <- eigen(G)
                        G <- outEig$vectors%*%diag(sapply(outEig$values,max,1e-5))%*%t(outEig$vectors)
                    }
                    )
                    GDa<-G
                    GDa[Ind]<-(-3*(a*Dist[Ind]-C3/a))
                    GDb<-NA
                }
                else
                {
                    stop('Valid spatial parametric models are exp, gauss, or sph')
                }
                
            }
            
        }
        
        # Output
        
        out<-NULL
        out$G<-G
        out$GDa<-GDa
        out$GDb<-GDb
        
        return(out)
        
    }
    
    # -------------------------------------------------------------------------
    # Main function
    
    Z<-X[,Zcols,drop=F]

    m<-length(ni)
    n<-sum(ni)
    p<-ncol(X)
    q<-length(Zcols)
    qsq<-q*q
    nth<-(qsq+q)/2+1
    ind<-c(rep(FALSE,qsq),rep(TRUE,2))
    
    for (k in c(1:q))
    {
        for (j in c(1:k))
        {
            ind[(k-1)*q+j] <- TRUE
        }
    }
    
    W<-matrix(0,n,max(ni))
    SIGMA<-W
    
    # Starting values
    
    mTheta0<-rowMeans(Theta0)
    D<-matrix(0,qsq,1)
    D[ind[1:(length(ind)-2)]]<-mTheta0[1:(length(mTheta0)-1)]
    D<-matrix(D,q,q)
    triuD<-D; triuD[lower.tri(triuD,1)]<-0
    D <- D + t(triuD); rm(triuD) 
    phisq<-mTheta0[length(mTheta0)]
    L <- chol(D)  
    phi<-sqrt(phisq)
    nv<-nrow(Dist)
    b<-NULL
    
    if (model=='gauss')
    {
        if (nv>2)
        {
            ind<-c(ind,TRUE)
            b<-0.01
        }
        else
        {
            model<-'exp'
            Dist<-Dist*Dist
        }
    }
    a<-0.05
    theta<-c(as.vector(L),phi,a,b)
    
    # Iterations
    
    nit<-35
    lr<-matrix(0,1,nit)
    Y<-reshapeM(Y,m,ni,nv)
    r<-Y
    
    it<- 0
    tf<-TRUE
    while (tf)
    {
        it<-it+1
        # Computation of W = SIGMA^-1 and H.
        posi<-1; H<-0
        scInvD<-solve(D)*phisq # check this

        for (i in c(1:m))
        {
            posf<-posi+ni[i]-1
            Zi<-Z[posi:posf,,drop=F]
            Wi<-(diag(ni[i])-lme_mrdivide(Zi,t(Zi)%*%Zi+scInvD,ginv.do=F)%*%t(Zi))/phisq
            SIGMA[c(posi:posf),c(1:ni[i])]<-Zi%*%D%*%t(Zi)+diag(ni[i])*phisq
            W[c(posi:posf),c(1:ni[i])]<-Wi
            Xi<-X[posi:posf,,drop=F]
            H<-H+t(Xi)%*%Wi%*%Xi
            posi<-posf+1
        }
        
        #invH <- lme_mldivide(H,diag(p),ginv.do=F)
        invH <- solve(H) # check this
        
        # Estimation
        
        posi<-1; posiY<-1; Bhat<-0
        for (i in c(1:m))
        {
            posf<-posi+ni[i]-1
            posfY<-posiY+ni[i]*nv-1
            Xi<-X[c(posi:posf),,drop=F]
            Wi<-W[c(posi:posf),c(1:ni[i])]
            Bhat <- Bhat + kronecker(diag(nv),invH%*%t(Xi)%*%Wi)%*%Y[posiY:posfY,drop=F] # KD
            posi<-posf+1
            posiY<-posfY+1
        }
            
        # spatial derivative
        outSptDer<-sptDer(Dist,a,b,model)
        G<-outSptDer$G
        GDa<-outSptDer$GDa
        GDb<-outSptDer$GDb
            
        # log-likelihood and residuals
        posi<-1
        posiY<-1
        lreml<-0
        term<-0
            
        for (i in c(1:m))
        {
            posf<-posi+ni[i]-1
            posfY<-posiY+ni[i]*nv-1
            Xi<-X[c(posi:posf),,drop=F]
            Wi<-W[c(posi:posf),c(1:ni[i]),drop=F]
            SIGMAi<-SIGMA[c(posi:posf),c(1:ni[i]),drop=F]
            ri <- Y[posiY:posfY,drop=F]-kronecker(diag(nv),Xi)%*%Bhat
            r[c(posiY:posfY)]<-ri
            term <- term + t(ri)%*%lme_mldivide(kronecker(G,SIGMAi),ri,ginv.do=FALSE) # check this
            lreml<-lreml+log(det(Wi)) 
            posi<-posf+1
            posiY<-posfY+1
        }

        lreml = 0.5*((n-p)*log(1/det(G)) + nv*(lreml-log(det(H))) - term) 

        # Fisher scoring
        
        svdG<-svd(G)
        if ((svdG$d[1]/svdG$d[length(svdG$d)])<1e+10) 
        {
            invG <- lme_mldivide(G,diag(nv),ginv.do=F) 
            gr<-lme_mass_RgGradient(X,Zcols,W,invH,L,phi,r,ni,invG,GDa,GDb)
            outRgEI<-lme_mass_RgEI(X,Zcols,W,invH,L,phi,ni,invG,GDa,GDb)
            EI<-outRgEI$EI
            Pth<-outRgEI$Pth
            Qthth<-outRgEI$Qthth
        }
        else
        {
            gr<-lme_mass_RgGradient1(X,Zcols,SIGMA,W,invH,L,phi,r,ni,G,GDa,GDb)
            outRgEI1<-lme_mass_RgEI1(X,Zcols,W,invH,L,phi,ni,G,GDa,GDb)
            EI<-outRgEI1$EI
            Pth<-outRgEI1$Pth
            Qthth<-outRgEI1$Qthth
        }

        svdEI<-svd(EI)
        if ((svdEI$d[1]/svdEI$d[length(svdEI$d)])<1e+10) 
        {
            invEI <- lme_mldivide(EI,diag(nrow(EI)),ginv.do=F)
        }
        else 
        {
            outEig <- eigen(EI)
            invEI <- outEig$vectors%*%diag(1/sapply(outEig$values,max,1e-5))%*%t(outEig$vectors) # check this
        }
        theta[ind]<-theta[ind]+invEI%*%gr
        gnorm<-norm(gr,type="2") 
        lr[it]<-lreml

        # Termination
        
        if ((it==nit) || (gnorm<e))
        {
            tf<-FALSE
            
            if (gnorm<e)
            {
                st<-matrix(1,nv,1)
            }
            else
            {
                st<-matrix(0,nv,1)
            }
            
            invEI<-invEI[c(1:nth),c(1:nth)]
            
            stats<-NULL
            for (i in c(1:nv))
            {
                stats[[i]] <- list(Bhat=matrix(0,p,1),CovBhat=invH,phisqhat=phisq,Dhat=D,Zcols=Zcols,invEI=invEI,Pth=Pth,Qthth=Qthth,lreml=lr[1:it])
                stats[[i]]$Bhat<-Bhat[((i-1)*p+1):(i*p)] 
            }
         }
        else
        {
            L <- matrix(theta[1:qsq],q,q)
            D<-t(L)%*%L
            phi<-theta[qsq+1]
            phisq<-phi*phi
            a<-theta[qsq+2]
            if (model=='gauss') b<-theta[length(theta)]
        }

    }
    
    # Output
    
    out<-NULL
    out$stats<-stats
    out$st<-st
    out$a<-a
    out$b<-b
    
    return(out)

}
