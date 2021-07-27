lme_fit_init<-function(X,Zcols,y,ni)
{
    library(MASS)
    
    # --------------------------------------------------------------------------
    # Initialize variables
    
    # descriptor variables
    m<-length(ni)
    n<-sum(ni)
    p<-ncol(X)
    q<-length(Zcols)
    
    # matrix of random effects
    Z<-X[,Zcols,drop=F]
    
    # estimated vector of the population regression parameters
    Bhat<-ginv(X,tol=0)%*%y
    
    # --------------------------------------------------------------------------
    #
    
    # auxiliary variables
    t1<-0
    t2<-0
    t3<-0
    t4<-0
    
    # position indicator
    posi<-1
    
    for(i in 1:m)
    {   
        # select subject-specific Z_i, X_i, y_i, and compute r_i
        posf<-posi + ni[i]-1
        Zi<-Z[posi:posf,,drop=F]
        Xi<-X[posi:posf,,drop=F]
        yi<-y[posi:posf,,drop=F]
        ri<-yi - Xi %*% Bhat
        
        #
        t<-ginv(t(Zi) %*% Zi)
        t1<-t1 + t
        bihat<-t %*% t(Zi) %*% ri
        t2<-t2 + t(yi) %*% yi- t(bihat) %*% t(Zi) %*% ri
        t3<-t3 + t(Xi) %*% yi
        t4<-t4 + bihat %*% t(bihat)
        
        #
        posi<-posf + 1
    }
    
    # 
    phisq0<-(t2 - t(Bhat) %*% t3) / (n - (m - 1) * q - p)
    phisq0<-as.vector(phisq0)

    if (phisq0<=0) phisq0<-1
    
    #
    D0<-(t4 - (phisq0 * t1)) / m
    
    # Handling non-positive definite initial D
    if(min(Re(eigen(D0)$values))<0)
    {
        EV<-(eigen(D0)$vectors)
        EV<-EV * -1
        ED<-eigen(D0)$values
        ED[ED < 0]<-10^-4
        ED<-diag(x=ED,nrow=length(ED),ncol=length(ED))
        D0<-EV %*% ED %*% (solve(EV))
    }
    
    return(list(D0=D0,phisq0=phisq0,Z=Z))
}