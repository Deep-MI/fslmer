lme_mass_fit_init<-function(X,Zcols,Y,ni,maskvtx=NA,numcore=1)
{
   
    if (numcore==1) print("No parallel computing enabled (not recommended)",quote=F)
    
    #
    
    nv0<-ncol(Y)
    
    if (any(is.na(maskvtx))) 
    {
        print("Selecting all vertices, possibly check 'maskvtx' argument",quote=F)
        maskvtx<-c(1:nv0)
    }
    
    Y<-Y[,maskvtx,drop=F]
    nv<-ncol(Y)
    
    m<-length(ni)
    n<-sum(ni)

    p<-ncol(X)
    q<-length(Zcols)

    ind<-rep(FALSE,q*q)

    for (k in c(1:q))
    {
        for (j in c(1:k))
        {
            ind[(k-1)*q+j]<-TRUE
        }
    }

    nth<-q*(q+1)/2+1
    
    Theta1<-matrix(0,nth,nv)
    phisqd<-(n-(m-1)*q-p)
    
    posi<-1
    
    t<-matrix(0,m*q,q)
    t1<-0
    
    for (i in c(1:m))
    {
        posf<-posi+ni[i]-1
        Zi<-X[posi:posf,Zcols,drop=F]
        t2<-MASS::ginv(t(Zi)%*%Zi) 
        t[c(((i-1)*q+1):(i*q)),]<-t2
        t1<-t1+t2
        posi<-posf+1
    }
    
    #
    
    print('Computing initial values ...',quote=F)
    
    sinvX<-MASS::ginv(X)
    Bhat<-sinvX%*%Y
    
    # define function for parallel processing    
    
    compute_theta<-function(b,y,X,Zcols,phisqd,m)
    {
        # Estimation

        t2<-0; t3<-0; t4<-0
        
        posi<-1
        
        for (i in c(1:m))
        {
            posf<-posi+ni[i]-1
            Zi<-X[c(posi:posf),Zcols,drop=F]
            Xi<-X[c(posi:posf),,drop=F]
            yi<-y[c(posi:posf)]
            ri<-yi-Xi%*%b
            bihat<-t[c((i-1)*q+1):(i*q),,drop=F]%*%t(Zi)%*%ri
            t2<-t2+t(yi)%*%yi-t(bihat)%*%t(Zi)%*%ri
            t3<-t3+t(Xi)%*%yi
            t4<-t4+bihat%*%t(bihat)
            posi<-posf+1
        }

        phisq<-as.vector((t2-t(b)%*%t3))/phisqd
        
        if (phisq<=0) phisq<-1
    
        D<-(t4-phisq*t1)/m
        
        D<-tryCatch({chol(D);D},error=function(e){
            eigRes<-eigen(D)
            eigRes$values[eigRes$values<0]<-(10^-4)
            D<-eigRes$vectors%*%diag(eigRes$values,nrow=length(eigRes$values),ncol=length(eigRes$values))%*%solve(eigRes$vectors) # check this
            })
    
        Theta1j<-c(as.vector(D)[ind],phisq)
        
        return(Theta1j)
    
    }
    
    #
    
    if (numcore==1)
    {      
        for (j in c(1:nv))
        {
            Theta1[,j]<-compute_theta(Bhat[,j,drop=F],Y[,j,drop=F],X,Zcols,phisqd,m)
        }
    }
    else
    {
        if (numcore>1)
        {
            Theta1<-simplify2array(parallel::mclapply(c(1:nv),function(x) compute_theta(Bhat[,x,drop=F],Y[,x,drop=F],X,Zcols,phisqd,m),mc.cores=numcore))
        }
    }

    
    #

    Theta0<-matrix(0,nth,nv0)
    Theta0[,maskvtx]<-Theta1
    Re0<-matrix(0,n,nv0)
    Re0[,maskvtx]<-Y-X%*%Bhat

    # output
    
    output<-NULL
    output$Theta0<-Theta0
    output$Re0<-Re0
    output$Bhat0<-Bhat
    
    return(output)
}
