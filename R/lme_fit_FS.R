lme_fit_FS<-function(X,Zcols,y,ni,e=10^-3)
    {
    # --------------------------------------------------------------------------
    # Input checks
    
    if (!exists("X")) stop("")
    if (!exists("Zcols")) stop("")
    if (!exists("y")) stop("")
    if (!exists("ni")) stop("")
    
    if (!is.matrix(X)) stop("X must be a matrix")
    if (!is.matrix(y) || ncol(y)>1) stop("y must be a mx1 matrix")
    if (!is.matrix(ni) || ncol(ni)>1) stop("y must be a sx1 matrix")
    
    if (nrow(X)!=nrow(y)) stop("X and y must have the same number of rows")
    if (nrow(X)!=sum(ni)) stop("X and sum(ni) do not match")
        
    if (max(Zcols)>ncol(X)) stop("max(Zcols) may not exceed ncol(X)")
        
    # --------------------------------------------------------------------------
    # Initialize variables
  
    # number of iterations
    nit<-50
    
    # program status
    st<-1
    
    # descriptor variables
    m<-length(ni)
    n<-sum(ni)
    p<-ncol(X)
    q<-length(Zcols)
    
    # indicator vector for theta: q^2+1 fields
    ind<-c(rep(FALSE,times=q*q),TRUE)
    for (k in 1:q)
        {
        for (j in 1:k)
            {
            ind[(k-1)*q+j]<-TRUE
            }
        }
    
    # Estimated covariance of the subject especific random coefficients
    Cbihat<-matrix(0,q*m,q)
    
    # Inverses of the estimated marginal covariance matrices for each subject 
    # stacked in W
    W<-matrix(0,nrow=n,ncol=max(ni))
    
    # Estimated marginal covariance matrices for each subject stacked in SIGMA
    SIGMA<-matrix(0,nrow=n,ncol=max(ni))
  
    # --------------------------------------------------------------------------
    # Starting values
  
    list2env(lme_fit_init(X,Zcols,y,ni),env=environment())
    
    # Estimated random effects covariance matrix and its Cholesky decomposition
    D<-D0
    L<-chol(D)
    
    # Estimated within-subject variability
    phisq<-phisq0
    phi<-sqrt(phisq)
    
    # 
    theta<-c(L,phi)
  
    # --------------------------------------------------------------------------
    # Iterations
   
    # status and current iteration 
    tf<-TRUE
    it<-0
  
    while (tf)
        {
        it<-it+1
        
        #Computation of W = SIGMA^-1 and H
        H<-0
        Term<-0
        
        scInvD<-lme_mldivide(D,diag(q)*phisq,ginv.do=F)
        
        posi<-1
        for (i in 1:m)
            {
            posf<-posi+ni[i]-1
            Zi<-Z[posi:posf,,drop=F]
            
            # B = (t(Zi) %*% Zi) + scInvD
            # B = lme_mrdivide(Zi,B)
            # Wi = (diag(ni[i]) - (B %*% t(Zi))) / phisq
            
            Wi<-(diag(ni[i])-(lme_mrdivide(Zi,((t(Zi)%*% Zi)+scInvD),ginv.do=F)%*%t(Zi)))/phisq # check this
            
            W[posi:posf,1:ni[i]]<-Wi
            
            SIGMA[posi:posf,1:ni[i]]<-Zi%*%D%*%t(Zi)+diag(ni[i])*phisq
            Xi<-X[posi:posf,,drop=F]
            Ti<-t(Xi)%*%Wi
            H<-H+(Ti%*%Xi)
            Term<-Term+(Ti%*%y[posi:posf,,drop=F])
            posi<-posf+1
            }
        
        invH<-lme_mldivide(H,diag(p),ginv.do=F)
        
        #Estimation
        Bhat<-invH %*% Term
        r<-y - X %*% Bhat
        
        lreml<-0
        posi<-1
        for (i in 1:m)
            {
            posf<-posi + ni[i] - 1
            Wi<-W[posi:posf,1:ni[i],drop=F]
            ri<-r[posi:posf,,drop=F]
            lreml<-lreml + log(det(Wi)) - t(ri) %*% Wi %*% ri
            posi<-posf+1
            }
 
        # compute gradient; this will return 'gr'
        list2env(lme_Gradient(X,Zcols,W,invH,L,phi,r,ni),env=environment())
        
        # compute expected information; this will return 'EI' (and 'Pth', 'Qthth', but we don't need these)
        list2env(lme_EI(X,Zcols,W,invH,SIGMA,L,phi,ni),env=environment())
    
        # update theta
        theta[ind]<-theta[ind] + lme_mldivide(EI,gr,ginv.do=F) # check this
    
        # restricted log-likelihood
        lreml<-0.5 * (lreml - log(det(H)))
        cat("\nLikelihood at FS iteration", it ,":", round(lreml,digits = 4))
        
        # gradient norm
        eps<-sqrt(sum(gr^2))
        cat("\nGradient norm:",round(eps,digits = 7))

        # termination?
        if ((it == nit) | (eps<e))
            {
            tf<-FALSE
            posi<-1
            Cbihat<-matrix(0,q * m,q)
            bihat<-matrix(0,q,m)
            for (i in 1:m)
                {
                posf<-posi + ni[i] - 1
                Zi<-Z[posi:posf,,drop=F]
                Wi<-W[posi:posf,1:ni[i],drop=F]
                Xi<-X[posi:posf,,drop=F]
                Pi<-Wi - Wi %*% Xi %*% invH %*% t(Xi) %*% Wi
                ri<-r[posi:posf,drop=F]
                bihat[,i]<-D %*% t(Zi) %*% Wi %*% ri
                Cbihat[(((i - 1)* q) + 1):(i * q),]<-D - (D %*% t(Zi) %*% Pi %*% Zi %*% D)
                posi<-posf + 1
                }
      
            stats<-list("Bhat" = Bhat, "CovBhat" = invH, "bihat" = bihat, "Covbihat" = Cbihat, "phisqhat" = phisq, "SIGMA" = SIGMA, "W" = W, "Dhat" = D, "X" = X, "Zcols" = Zcols, "re" = r, "ni" = ni, "lreml" = lreml)
      
            if (rcond(EI) < 10^-11)
                {
                cat("\nMatrix EI may be close to singular. Results may be inaccurate. RCOND =",rcond(EI))
                }
      
            if (it == nit)
                {
                st<-0 
                cat("\nAlgorithm does not converge after",it,"iterations!!!");  
                }
            }  
        else 
            {
            L<-matrix(theta[1:(length(theta) - 1)],c(q,q))
            D<-t(L) %*% L
            phi<-theta[length(theta)]
            phisq<-phi * phi
            }
  
        }   
  
    return(stats)
    }
