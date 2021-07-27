lme_Gradient<-function(X,Zcols,W,invH,L,phi,re,ni)
{
    # descriptor variables
    m<-length(ni)
    q<-length(Zcols)
    
    # random effects
    Z<-X[,Zcols,drop=F]
    
    # 
    nth<-q * (q + 1) / 2 + 1
    gr<-rep(0,nth)
    gr[nth]<-1

    #
    jk<-0
    
    for (k in 1:q)
    { 
        
        for (j in 1:k)
        { 
            jk<-jk + 1
            posi<-1
            a1<-0
            a2<-0
            M1<-0
        
            for (i in 1:m)
            {  
                posf<-posi + ni[i,] - 1
                Zi<-Z[posi:posf,,drop=F]
                Wi<-W[posi:posf,1:ni[i],drop=F]
                Xi<-X[posi:posf,,drop=F]
                rei<-re[posi:posf]
                Mjki<-(Zi[,k] %*% t(L[j,])) %*% t(Zi)
                Mjki<-Mjki + t(Mjki)
                Maux<-Mjki %*% Wi
                a1<-a1 + sum(diag(-Maux))
                a2<-a2 + t(rei) %*% Wi %*% Maux %*% rei
                M1<-M1 + t(Xi) %*% Wi %*% Maux %*% Xi
                posi<-posf + 1
            }
            a3<-sum(diag(-invH %*% M1))
            gr[jk]<-(a1 + a2 - a3) / 2
        }
             
        #
        posi<-1
        a1<-0
        a2<-0
        M1<-0
        
        for (i in 1:m)
        { 
            posf<-posi + ni[i] - 1
            Wi<-W[posi:posf,1:ni[i],drop=F]
            Xi<-X[posi:posf,,drop=F]
            rei<-re[posi:posf,,drop=F]
            a1<-a1 + sum(diag(-Wi))
            Maux<-Wi %*% Wi
            a2<-a2 + t(rei) %*% Maux %*% rei
            M1<-M1 + t(Xi) %*% Maux %*% Xi
            posi<-posf+1
        }
        a3<-sum(diag(-invH %*% M1))
        gr[nth]<-phi * (a1 + a2 - a3)
    }
    
    return(list(gr=gr))
}