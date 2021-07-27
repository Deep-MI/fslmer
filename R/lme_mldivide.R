lme_mldivide<-function(A,B,ginv.do=TRUE,fallback=TRUE) 
{
    #
    
    stopifnot(is.numeric(A) || is.complex(A),is.numeric(B) || is.complex(B))
    
    #
    
    if (is.vector(A)) A <- as.matrix(A)
    if (is.vector(B)) B <- as.matrix(B)
    
    #
    
    if (nrow(A) != nrow(B))
    {
        stop("Matrices 'A' and 'B' must have the same number of rows.")
    }
    
    #
    
    if (ginv.do) 
    {
        out<-MASS::ginv(t(A) %*% A) %*% t(A) %*% B
    } 
    else 
    {
        st<-0
        tryCatch(
        {
            if (nrow(A)==ncol(A))
            {
                cA<-chol(A)
                invA<-backsolve(cA,forwardsolve(t(cA),diag(nrow(A))))
                out<-invA%*%B
                rm(cA,invA)
                st<-1
            }
            else
            {
                out<-qr.solve(A, B)
                st<-1
            }
        }
        ,error=function(e){}
        )
        if (st==0)
        {
            if (fallback)
            {
                if (nrow(A)==ncol(A))
                {
                    print("chol() failed, falling back to solve()")
                }
                else
                {
                    print("qr.solve() failed, falling back to solve()")
                }
                
                tryCatch(
                {
                    out<-solve(A, B)
                    st<-1
                }
                ,error=function(e){}
                )
                if (st==0)
                {
                    print("solve() failed, falling back to ginv()")
                    tryCatch(
                    {
                        out<-MASS::ginv(t(A) %*% A) %*% t(A) %*% B
                        st<-1
                    }
                    ,error=function(e){print("We have an error here. Aborting.")}
                    )
                }
            }
            else
            {
                print("All failed. Aborting.")  
            }
        }
    }
    
    #
    
    return(out)
    
}
