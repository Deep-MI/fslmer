lme_mass_RgGradient<-function(X,Zcols,W,invH,L,phi,re,ni,invG,GDa,GDb)
{

    m<-length(ni)
    q<-length(Zcols)
    Z<-X[,Zcols,drop=F]
    
    nth<-q*(q+1)/2+1
    nv<-nrow(invG)
    
    # Log-likelihoood derivative for L
    
    gr<-matrix(0,nth+1,1)
    jk<-0
    
    for (k in c(1:q))
    {
        for (j in c(1:k))
        {
            jk<-jk+1
            posi<-1; posir<-1
            a1<-0; a2<-0
            M1<-0
            
            for (i in c(1:m))
            {
                posf<-posi+ni[i]-1
                posfr<-posir+ni[i]*nv-1
                Zi<-Z[c(posi:posf),,drop=F]
                Wi<-W[c(posi:posf),c(1:ni[i]),drop=F]
                Xi<-X[c(posi:posf),,drop=F]
                rei<-re[c(posir:posfr)]
                Mjki<-Zi[,k,drop=F]%*%L[j,,drop=F]%*%t(Zi)
                Mjki<-Mjki+t(Mjki)
                Maux<-Mjki%*%Wi
                a1 <- a1 - sum(diag(Maux))
                a2 <- a2 + t(rei)%*%kronecker(invG,Wi%*%Maux)%*%rei
                M1<-M1+t(Xi)%*%Wi%*%Maux%*%Xi
                posi<-posf+1
                posir<-posfr+1
            }
            a3 <- (-sum(diag(invH%*%M1)))
            gr[jk]<-(nv*a1+a2-nv*a3)/2
        }
    }
    
    # Log-likelihoood derivative for phi
    
    posi<-1; posir<-1
    a1<-0; a2<-0
    M1<-0
    
    for (i in c(1:m))
    {
        posf<-posi+ni[i]-1
        posfr<-posir+ni[i]*nv-1;
        Wi<-W[c(posi:posf),c(1:ni[i]),drop=F]
        Xi<-X[c(posi:posf),,drop=F]
        rei<-re[c(posir:posfr)]
        a1 <- a1 - sum(diag(Wi))
        Maux<-Wi%*%Wi
        a2 <- a2 + t(rei)%*%kronecker(invG,Maux)%*%rei
        M1<-M1+t(Xi)%*%Maux%*%Xi
        posi<-posf+1
        posir<-posfr+1
    }
    
    a3 <- (-sum(diag(invH%*%M1)))
    gr[nth]<-phi*(nv*a1+a2-nv*a3)
    
    # Log-likelihoood derivative for a
    
    posi<-1; posir<-1
    a2<-0
    Maux1<-GDa%*%invG
    Maux2<-invG%*%Maux1
    
    for (i in c(1:m))
    {
        posf<-posi+ni[i]-1
        posfr<-posir+ni[i]*nv-1
        Wi<-W[c(posi:posf),c(1:ni[i]),drop=F]
        rei<-re[c(posir:posfr)]
        a2 <- a2 + t(rei)%*%kronecker(Maux2,Wi)%*%rei
        posi<-posf+1
        posir<-posfr+1
    }

    aux <- (-sum(diag(Maux1))) 
    a1<-sum(ni)*aux
    a3<-ncol(X)*aux
    gr[nth+1]<-(a1+a2-a3)/2

    if (!is.na(GDb))
    {
        # Log-likelihoood derivative for b
        posi<-1; posir<-1
        a2<-0
        Maux1<-GDb%*%invG
        Maux2<-invG%*%Maux1
        
        for (i in c(1:m))
        {
            posf<-posi+ni[i]-1
            posfr<-posir+ni[i]*nv-1
            Wi<-W[c(posi:posf),c(1:ni[i]),drop=F]
            rei<-re[c(posir:posfr)]
            a2 <- a2 + t(rei)%*%kronecker(Maux2,Wi)%*%rei
            posi<-posf+1
            posir<-posfr+1
        }

    aux <- (-sum(diag(Maux1)))
    a1<-sum(ni)*aux
    a3<-ncol(X)*aux
    gr(nth+2) <- (a1+a2-a3)/2
    }
    else
    {
        gr<-gr[c(1:(nth+1))]
    }
    
    #
    
    return(gr)
    
} 