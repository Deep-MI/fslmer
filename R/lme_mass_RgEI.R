lme_mass_RgEI<-function(X,Zcols,W,CBhat,L,phi,ni,invG,GDa,GDb)
{
    m<-length(ni)
    n<-sum(ni)
    nv<-nrow(invG)
    q<-length(Zcols)
    nth<- q*(q+1)/2+1
    p<-ncol(X)
    Pth <- array(0,c(nth,p,p))
    Qthth <- array(0,c(nth,nth,p,p))
    Der <- array(0,c(nth,n,max(ni))) 
    
    # Computation of the first order derivatives of the temporal covariance matrix
    
    jk<-0
    for (k in c(1:q))
    {
        for (j in c(1:k))
        {
            jk<-jk+1
            posi<-1
            for (i in c(1:m))
            {
                posf<-posi+ni[i]-1
                Zi<-X[c(posi:posf),Zcols,drop=F]
                Zki<-Zi[,k,drop=F]    
                Mjki<-Zki%*%L[j,,drop=F]%*%t(Zi)
                Mjki<-Mjki+t(Mjki)
                Der[jk,c(posi:posf),c(1:ni[i])]<-Mjki
                posi<-posf+1
            }
        }
    }
                
    posi<-1
    for (i in c(1:m))
    {
        posf<-posi+ni[i]-1
        Der[nth,c(posi:posf),c(1:ni[i])] <- 2*phi*diag((ni[i])) 
        posi<-posf+1
    }
    
    # Computation of Pis,Qijs and the expected information matrix EI
    
    for (j in c(1:nth))
    {
        posi<-1; Pj<-0
        Bj <- drop(Der[j,,])
        for (i in c(1:m))
        {
            posf<-posi+ni[i]-1
            Wi<-W[c(posi:posf),c(1:ni[i]),drop=F]
            Pj<-Pj-t(X[c(posi:posf),,drop=F])%*%Wi%*%Bj[c(posi:posf),c(1:ni[i]),drop=F]%*%Wi%*%X[c(posi:posf),,drop=F]
            posi<-posf+1
        }
        Pth[j,,]<-Pj
    }
    
    EI<-matrix(0,nth+2,nth+2)

    # Expected information among Lijs (including phi)
    
    for (k in c(1:nth))
    {
        Bk <- drop(Der[k,,])
        Pk <- drop(Pth[k,,])
        for (j in c(1:k))
        {
            Bj <- drop(Der[j,,])
            Pj <- drop(Pth[j,,])
            posi<-1; Qkj<-0
            traceBkj<-0
            for (i in c(1:m))
            {
                posf<-posi+ni[i]-1
                Wi<-W[c(posi:posf),c(1:ni[i]),drop=F]
                Bkji<-Wi%*%Bk[c(posi:posf),c(1:ni[i]),drop=F]%*%Wi%*%Bj[c(posi:posf),c(1:ni[i]),drop=F]
                traceBkj <- traceBkj + sum(diag(Bkji))  
                Bkji<-Bkji%*%Wi
                Qkji<-t(X[c(posi:posf),,drop=F])%*%Bkji%*%X[c(posi:posf),,drop=F]
                Qkj<-Qkj+Qkji
                posi<-posf+1
            }
            Qthth[k,j,,]<-Qkj
            Qthth[j,k,,]<-Qkj
            EI[k,j] <- nv*(traceBkj - sum(diag(CBhat%*%(2*Qkj-Pk%*%CBhat%*%Pj))))
            EI[j,k]<-EI[k,j]
        }
    
    }
            
    # Expected information between Lijs (including phi) and a (first spatial parameter)

    Mauxa<-GDa%*%invG
    trMauxa <- sum(diag(Mauxa))
    
    for (k in c(1:nth))
    {
        Bk <- drop(Der[k,,])
        Pk <- drop(Pth[k,,])
        posi<-1
        traceBk<-0
        for (i in c(1:m))
        {
            posf<-posi+ni[i]-1
            Wi<-W[c(posi:posf),c(1:ni[i]),drop=F]
            Bki<-Wi%*%Bk[c(posi:posf),c(1:ni[i]),drop=F]
            traceBk <- traceBk + sum(diag(Bki))
            posi<-posf+1
        }
        EI[k,nth+1] <- trMauxa*(traceBk + sum(diag(CBhat%*%Pk)))
        EI[nth+1,k]<-EI[k,nth+1]
    }
    
    # Expected information between a and a
    
    EI[nth+1,nth+1] <- (n-p)*sum(diag(Mauxa%*%Mauxa))
    
    if (!is.na(GDb))
    {
        # Expected information between Lijs (including phi) and b (second spatial parameter)
        
        Mauxb<-GDb%*%invG
        trMauxb <- sum(diag(Mauxb)) # check this 
        for (k in c(1:nth))
        {
            Bk <- drop(Der[k,,]) # check this
            Pk <- drop(Pth[k,,]) # check this
            posi<-1
            traceBk<-0
            for (i in c(1:m))
            {
                posf<-posi+ni[i]-1
                Wi<-W[c(posi:posf),c(1:ni[i]),drop=F]
                Bki<-Wi%*%Bk[c(posi:posf),c(1:ni[i]),drop=F]
                traceBk <- traceBk + sum(diag(Bki))
                posi<-posf+1
            }
            EI[k,nth+2] <- trMauxb*(traceBk + sum(diag(CBhat%*%Pk)))
            EI[nth+2,k]<-EI[k,nth+2]
        }
        
        # Expected information between b and b
        EI[nth+2,nth+2] <- (n-p)*sum(diag(Mauxb%*%Mauxb))
        
        # Expected information between a and b 
        EI[nth+1,nth+2] <- (n-p)*sum(diag(Mauxa%*%Mauxb))
        EI[nth+2,nth+1]<-EI[nth+1,nth+2]
    }
    else
    {
        EI<-EI[c(1:(nth+1)),c(1:(nth+1))] 
    }
    
    EI<-0.5*EI
    
    # Output
    
    out<-NULL
    out$EI<-EI
    out$Pth<-Pth
    out$Qthth<-Qthth
    
    return(out)

}