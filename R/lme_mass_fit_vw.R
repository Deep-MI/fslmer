lme_mass_fit_vw<-function(X,Zcols,Y,ni,maskvtx=NA,fname=NA,prs=1,e=10^-1,Xrows=NA,numcore=1)
{
    stats<-NULL
    
    # check if parallel computing is feasible
    
    if (numcore>1 & !exists("mclapply")) 
    {
        st<-require(parallel)
        
        if (!st) stop("Please install the 'parallel' package or set numcore=1")
        
        print("loading 'parallel' package")
    
    }
    
    if (numcore==1) print("No parallel computing enabled (not recommended)",quote=F)
    
    #
    
    nv0<-ncol(Y)

    if (any(is.na(maskvtx)))
    {
        print(paste("Selecting all available vertices (",nv0,").",sep=""),quote=F)
        maskvtx<-c(1:nv0)
    }

    Y<-Y[,maskvtx,drop=F]
    
    p<-ncol(X)
    
    nv<-ncol(Y)
    
    rfcols<-Zcols
    Zcols<-matrix(0,1,p)
    Zcols[,rfcols]<-matrix(1,1,length(rfcols))
    
    if (!is.na(Xrows))
    {
        Xrows<-Xrows[,maskvtx]
    }
    
    stats1<-lme_mass_fit(X,NA,Xrows,Zcols,Y,ni,prs,e,numcore)
    
    for (i in c(1:nv))
    {
        stats[maskvtx[i]]<-stats1[i]
    }        
    
    return(stats)
}