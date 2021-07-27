lme_EI<-function(X,Zcols,W,CBhat,SIGMA,L,phi,ni)
{
    # Descriptor variables
    
    m<-length(ni)
    n<-sum(ni)
    p<-ncol(X)
    q<-length(Zcols)
    
    nth<-(q * (q + 1) / 2) + 1
    Pth<-array(0,dim<-c(nth,p,p))
    Qthth<-array(0,dim<-c(nth,nth,p,p))
    Der<-array(0,dim<-c(nth,n,ncol(W)))
    
    # Computation of the first order derivatives of W
    
    jk<-0
    for (k in 1:q)
    {
        
        for (j in 1:k)
        {    
            jk<-jk + 1
            posi<-1
            for (i in 1:m)
            {
                posf<-posi + ni[i] - 1
                Zi<-X[posi:posf,Zcols,drop=F]
                Zki<-Zi[,k,drop=F]
                Mjki<-Zki %*% L[j,] %*% t(Zi)
                Mjki<-Mjki + t(Mjki)
                Wi<-W[posi:posf,1:ni[i],drop=F]
                Der[jk,posi:posf,1:ni[i]]<--Wi %*% (Mjki %*% Wi)
                posi<-posf + 1
            }
        }
    }
    
    posi<-1
    for (i in 1:m)
    {
        posf<-posi + ni[i] - 1
        Wi<-W[posi:posf,1:ni[i],drop=F]
        Der[nth,posi:posf,1:ni[i]]<--2 * phi * Wi %*% Wi
        posi<-posf + 1
    }
    
    # Computation of Pis, Qijs and the expected information matrix EI.
    
    for (j in 1:nth)
    {
        posi<-1
        Pj<-0
        Bj<-drop(Der[j,,])
     
        for (i in 1:m)
        {  
            posf<-posi + ni[i] - 1
            Pj<-Pj + t(X[posi:posf,,drop=F]) %*% Bj[posi:posf,1:ni[i],drop=F] %*% X[posi:posf,,drop=F]
            posi<-posf + 1
        }
        Pth[j,,]<-Pj
    }
    
    EI<-matrix(0,nth,nth)
    
    for (k in 1:nth)
    {
        Bk<-drop(Der[k,,])
        Pk<-drop(Pth[k,,])
        
        for (j in 1:k)                                 
        {	
            Bj<-drop(Der[j,,])
            Pj<-drop(Pth[j,,])
            posi<-1
            Qkj<-0
            traceBkj<-0
            
            for (i in 1:m)
            {
                posf<-posi + ni[i] - 1
                Vi<-SIGMA[posi:posf,1:ni[i],drop=F]
                Bkji<-Bk[posi:posf,1:ni[i],drop=F] %*% Vi %*% Bj[posi:posf,1:ni[i],drop=F]
                traceBkj<-traceBkj + sum(diag(Bkji %*% Vi))
                Qkji<-t(X[posi:posf,,drop=F]) %*% Bkji %*% X[posi:posf,,drop=F]
                Qkj<-Qkj + Qkji
                posi<-posf + 1
            }
            Qthth[k,j,,]<-Qkj
            Qthth[j,k,,]<-Qkj
            EI[k,j]<-.5 * (traceBkj-sum(diag((CBhat %*% (2 * Qkj - Pk %*% CBhat %*% Pj)))))
            EI[j,k]<-EI[k,j]
        }
    }
    return(list(EI=EI,Pth=Pth,Qthth=Qthth))
}