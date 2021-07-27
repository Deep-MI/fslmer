lme_mrdivide<-function(A,B,ginv.do=TRUE) 
{
    #
    
    stopifnot(is.numeric(A) || is.complex(A),is.numeric(B) || is.complex(B))
    
    #
    
    if (is.vector(A)) A <- t(A)
    if (is.vector(B)) B <- t(B)
    
    #
    
    if (ncol(A) != ncol(B))
    {
        stop("Matrices 'A' and 'B' must have the same number of columns.")
    }
    
    #
    
    out<-t(lme_mldivide(t(B),t(A),ginv.do=ginv.do))
    
    #
    
    return(out)
    
}