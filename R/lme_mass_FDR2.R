lme_mass_FDR2<-function(pval,sgn,maskvtx=NA,rate=0.05,tail=0)
{
    # subfunction

    sidedPval<-function(pval,tail)
    {
        spval<-abs(pval)
        
        if (tail==-1)
        {
            spval[pval<0]<-spval[pval<0]*0.5
            spval[pval>0]<-1-0.5*spval[pval>0]
        }
        if (tail==1)
        {
            spval[pval>0]<-spval[pval>0]*0.5
            spval[pval<0]<-1-0.5*spval[pval<0]
        }
        if (!(tail==0 || tail==-1 || tail==1))
        {
            stop(print("Tail must be -1 for left-sided or 0 for two-sided or 1 for right-sided hypothesis testing"))
        }
        
        return(spval)
    }
    
    #

    nv0<-length(pval)
    if (is.na(maskvtx)) {maskvtx<-c(1:nv0)}
    
    p<-pval[maskvtx]*sgn[maskvtx]
    p[pval[maskvtx]==1]<-1
    nv<-length(p)
    
    # first stage
    
    q0<-rate/(1+rate)
    spval<-sidedPval(p,tail)
    pth0<-lme_mass_FDR(spval,q0)
    
    if (tail==0) 
    {
        # two-sided thresh
        detv0<-maskvtx[spval<=pth0]
    }
    if  (tail==-1)
    {
        # left-sided thresh
        vtx<-maskvtx[p<0]
        detv0<-vtx[spval[p<0]<=pth0]
    }
    if (tail==1)
    {
        # right-sided thresh
        vtx<-maskvtx[p>0]
        detv0<-vtx[spval[p>0]<=pth0]
    }
    
    ndetv0<-length(detv0)
    m0<-nv-ndetv0
    
    # Second stage
    
    if ((ndetv0!=0) && (ndetv0!=nv))
    {                
        # one sided-thresh
        pth<-lme_mass_FDR(spval,q0*nv/m0)
    
        if (tail==0)
        {
            # two-sided thresh
            detvtx<-maskvtx[spval<=pth]
        }
        if  (tail==-1)
        {
            # left-sided thresh
            detvtx<-vtx[spval[p<0]<=pth]
        }
        if (tail==1)
        {            
            # right-sided thresh
            detvtx<-vtx[spval[p>0]<= pth]
        }
    }
    else
    {
        detvtx<-detv0
        pth<-pth0
    }
    
    sided_pval<-matrix(1,1,nv0)
    sided_pval[maskvtx]<-spval
    
    out=list(detvtx=detvtx,sided_pval=sided_pval,pth=pth,m0=m0)
    return(out)
    
}
