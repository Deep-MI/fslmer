lme_mass_FDR<-function(p,fdr)
{
    p<-sort(abs(p))
    nv<-length(p)
    
    nn<-matrix(c(1:nv),ncol=1)
    
    if (length(which(p<=(fdr*nn/nv)))>0)
    {
        imax<-max(which(p<=(fdr*nn/nv)))
        pthresh<-p[imax]
    }
    else
    {
        pthresh<-min(p)/10
    }
        
    return(pthresh)    
}