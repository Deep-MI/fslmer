lme_FSfit<-function(X,Zcols,y,ni,e=0.1)
{
    stats<-NULL
    st<-(-1)

    # --------------------------------------------------------------------------
    # Initialize variables
    
    nit<-20
    m<-length(ni)
    p<-ncol(X)
    q<-length(Zcols)
    nth<-q*(q+1)/2+1
    ind<-c(rep(FALSE,times=q*q),TRUE)
    for (k in 1:q)
    {
        for (j in 1:k)
        {
            ind[(k-1)*q+j]<-TRUE
        }
    }
    
    n<-sum(ni)
    W<-matrix(0,nrow=n,ncol=max(ni))
    SIGMA<-W
    lr<-rep(0,nrow=1,ncol=nit)

    # --------------------------------------------------------------------------
    # Starting values
    
    list2env(lme_fit_init(X,Zcols,y,ni),envir=environment())
    
    D<-D0
    phisq<-phisq0
    L<-chol(D)
    phi<-sqrt(phisq)
    theta<-c(L,phi)    
    
    # --------------------------------------------------------------------------
    # Iterations
    
    print("Starting Fisher scoring iterations",quote=F)
    
    tf<-TRUE
    it<-0
    while (tf)
    {
        it<-it+1
        
        # Computation of W = SIGMA^-1 and H
        posi<-1
        H<-0
        Term<-0
        scInvD<-lme_mldivide(D,diag(q),ginv.do=F)*phisq
        for (i in 1:m)
        {
            posf<-posi+ni[i]-1
            Zi<-Z[posi:posf,,drop=F]
            Wi<-(diag(ni[i])-lme_mrdivide(Zi,(t(Zi)%*%Zi+scInvD),ginv.do=F)%*%t(Zi))/phisq 
            W[posi:posf,1:ni[i]]<-Wi
            SIGMA[posi:posf,1:ni[i]]<-Zi %*% D %*% t(Zi)+diag(ni[i])*phisq
            Xi<-X[posi:posf,,drop=F]
            Ti<-t(Xi) %*% Wi
            H<-H+(Ti %*% Xi)
            Term<-Term+(Ti %*% y[posi:posf,,drop=F])
            posi<-posf+1
        }
        
        invH<-lme_mldivide(H,diag(p),ginv.do=F)
        
        # Estimation
        Bhat<-invH %*% Term
        r<-y-X %*% Bhat
        posi<-1
        lreml<-0
        
        for (i in 1:m)
        {
            posf<-posi+ni[i]-1
            Wi<-W[posi:posf,1:ni[i],drop=F]
            ri<-r[posi:posf,,drop=F]
            lreml<-lreml+log(det(Wi))-t(ri) %*% Wi %*% ri
            posi<-posf+1
        }
        
        lreml<-0.5*(lreml-log(det(H)))
        
        list2env(lme_Gradient(X,Zcols,W,invH,L,phi,r,ni),envir=environment())
        gr<-as.matrix(gr)
        
        list2env(lme_EI(X,Zcols,W,invH,SIGMA,L,phi,ni),envir=environment())
        
        invEI<-lme_mldivide(EI,diag(nth),ginv.do=F)
        
        theta[ind]<-theta[ind]+invEI%*%gr
        
        eps<-norm(gr,type="2")
        
        lr[it]<-lreml
        
        # Termination
        
        if ((it == nit) | (eps<e))
        {
            tf<-FALSE
            stats<-list("Bhat" = Bhat, "CovBhat" = invH, "phisqhat" = phisq,"Dhat" = D, "X" = X, "Zcols" = Zcols, "invEI" = invEI, "Pth" = Pth, "Qthth" = Qthth, "lreml" = lr[1:it])
            if (it == nit)
            {
                st<-0 
            }      
            else
            {
                st<-1
            }
        }
        else
        {
            L<-matrix(theta[1:(length(theta)-1)],c(q,q))
            D<-t(L) %*% L
            phi<-theta[length(theta)]
            phisq<-phi*phi;                   
        }
        
    }
        
    stats$st<-st
    return(stats)

}
