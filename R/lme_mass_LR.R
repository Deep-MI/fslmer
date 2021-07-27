lme_mass_LR<-function(statsfull,statsred,q)
{
    # --------------------------------------------------------------------------
    
    if (exists("nargin") == TRUE) 
    {
        if (nargin < 3) 
        {
            stop(cat("Error: Too Few Inputs"))
        } 
    } 
    
    # --------------------------------------------------------------------------
    
    nv<-length(statsfull)
    pval<-matrix(0,nrow=1,ncol=nv)
    
    for (i in c(1:nv))
    {
        lrstat<-lme_LR(statsfull[i]$lreml,statsred[i]$lreml,q)
        pval[i]<-lrstats$pval
    }
    
    return(pval)    
}