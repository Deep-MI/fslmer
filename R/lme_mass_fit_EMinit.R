lme_mass_fit_EMinit<-function(X,Zcols,Y,ni,maskvtx=NA,nit=5,numcore=1)
{
    
    #
    
    nv0<-ncol(Y)
    
    if (is.na(maskvtx)) maskvtx<-c(1:nv0)
    
    Y<-Y[,maskvtx,drop=F]
    nv<-ncol(Y)
    
    m<-length(ni)
    n<-sum(ni)
    
    p<-ncol(X)
    q<-length(Zcols)
    
    nimax<-max(ni)

    nth<-q*(q+1)/2+1
    nDth<-q*q
    
    ind<-matrix(FALSE,nDth,1)
    
    for (k in c(1:q))
    {
        for (j in c(1:k))
        {
            ind[(k-1)*q+j]<-TRUE
        }
    }
    
    outInit<-lme_mass_fit_init(X,Zcols,Y,ni,maskvtx=NA,numcore)
    Theta1<-outInit$Theta0
    Re1<-outInit$Re0
    
    #
    
    print('Starting Expectation Maximization iterations ...',quote=F)
    
    for (j in c(1:nv))
    {
        
        # Iterations
        W<-matrix(0,n,nimax)
        ths<-Theta1[,j,drop=F]
        D<-matrix(0,nDth,1)
        D[ind]<-ths[1:(length(ths)-1)]
        D<-matrix(D,q,q)
        Dup<-D
        Dup[lower.tri(Dup,diag=T)]<-0
        D<-D+t(Dup)
        phisq<-ths[length(ths)]

        for (k in c(1:nit))
        {

            # Computation of W = SIGMA^-1 and H.
            posi<-1; H<-0; Term<-0
            y<-Y[,j,drop=F]
            scInvD<-lme_mldivide(D,diag(q),ginv.do=F)*phisq # check this
        
            for (i in c(1:m))
            {
                posf<-posi+ni[i]-1
                Zi<-X[c(posi:posf),Zcols,drop=F]
                Wi<-(diag(ni[i])-lme_mrdivide(Zi,(t(Zi)%*%Zi+scInvD),ginv.do=F)%*%t(Zi))/phisq # check this
                W[c(posi:posf),c(1:ni[i])]<-Wi
                Xi<-X[c(posi:posf),,drop=F]
                Ti<-t(Xi)%*%Wi
                H<-H+Ti%*%Xi
                Term<-Term+Ti%*%y[c(posi:posf)]
                posi<-posf+1
            }

            invH<-lme_mldivide(H,diag(p),ginv.do=F) # check this
    
            # Estimation
            
            Bhat<-invH%*%Term
            r<-y-X%*%Bhat
            Re1[,j]<-r
            Dhat<-matrix(0,q,q)
            phisqhat<-0
    
            posi<-1
            for (i in c(1:m))
            {
                posf<-posi+ni[i]-1
                Zi<-X[c(posi:posf),Zcols,drop=F]
                Wi<-W[c(posi:posf),c(1:ni[i])]
                Xi<-X[c(posi:posf),,drop=F]
                ri<-r[c(posi:posf)]
                bihat<-D%*%t(Zi)%*%Wi%*%ri
                Pi<-Wi-Wi%*%Xi%*%invH%*%t(Xi)%*%Wi
                phisqhat<-phisqhat+t(ri-Zi%*%bihat)%*%(ri-Zi%*%bihat)+phisq*sum(diag(diag(ni[i])-phisq*Pi))
                Dhat<-Dhat+bihat%*%t(bihat)+D-D%*%t(Zi)%*%Pi%*%Zi%*%D
                posi<-posf+1
            }
        
            D<-Dhat/m
            phisq<-as.vector(phisqhat/n)
        }
    
        Theta1[,j]<-c(as.vector(D)[ind],phisq)
        
    }

    Theta0<-matrix(0,nth,nv0)
    Theta0[,maskvtx]<-Theta1
    Re0<-matrix(0,n,nv0)
    Re0[,maskvtx]<-Re1
    
    # output
    
    out<-NULL
    out$Theta0<-Theta0
    out$Re0<-Re0

    return(out)
}
