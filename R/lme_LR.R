lme_R<-function(lrmlfull,lrmlred,q)
{
    # --------------------------------------------------------------------------
    # preparations
    
    # check inputs
    
    if (exists("nargin") == TRUE) 
    {
        if (nargin < 3) 
        {
            stop(cat("Error: Too Few Inputs"))
        } 
    } 
    
    # --------------------------------------------------------------------------
    
    G<-2*(lrmlfull-lrmlred)
    
    n1<-1000000
    n2<-1000000
    
    lrstats<-NULL
    lrstats$pval<-sum(c(rchisq(n=n1,df=q),rchisq(n=n2,df=(q+1)))>G)/(n1+n2)
    lrstats$G<-G
    lrstats$df<-c(q,q+1)
    
    return(lrstats)
}