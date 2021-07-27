lme_mass_AdjMtx<-function(Surf,maskvtx=NA)
{
    # -------------------------------------------------------------------------
    # Auxiliary functions
    
    add<-function(v,last,el)
    {
        if (sum(v==el)==0)
        {
            v[last]<-el
            last<-last+1
        }
    
        out<-NULL
        out$v<-v
        out$last<-last
        
        return(out)
    }
 
    # -------------------------------------------------------------------------
    # Main function

    surftri<-Surf$faces
    cn<-max(table(matrix(surftri,1)))
    
    if (!is.na(maskvtx))
    {
        
        logictri<-surftri
        
        for (i in c(1:3))
        {
    
            logictri[,i]<-is.element(surftri[,i],maskvtx)
        }
        
        surftri<-surftri[rowSums(logictri)==3,,drop=F] 
    }
    
    nv<-ncol(Surf$vertices)
    
    AdjM<-matrix(0,nv,cn)
    last<-matrix(1,nv,1)
    ntri<-nrow(surftri)
    
    if (ntri>0) {
    
        for (i in c(1:ntri))
        {
            tri<-surftri[i,,drop=F]
            
            tmp<-add(AdjM[tri[1],,drop=F],last[tri[1]],tri[2])
            AdjM[tri[1],]<-tmp$v
            last[tri[1]]<-tmp$last
            rm(tmp)
            
            tmp<-add(AdjM[tri[1],,drop=F],last[tri[1]],tri[3])
            
            AdjM[tri[1],]<-tmp$v
            last[tri[1]]<-tmp$last
            rm(tmp)
                
            tmp<-add(AdjM[tri[2],,drop=F],last[tri[2]],tri[1])                
            AdjM[tri[2],]<-tmp$v
            last[tri[2]]<-tmp$last
            rm(tmp)
                
            tmp<-add(AdjM[tri[2],,drop=F],last[tri[2]],tri[3])
            AdjM[tri[2],]<-tmp$v
            last[tri[2]]<-tmp$last
            rm(tmp)
                
            tmp<-add(AdjM[tri[3],,drop=F],last[tri[3]],tri[1])
            AdjM[tri[3],]<-tmp$v
            last[tri[3]]<-tmp$last
            rm(tmp)
                
            tmp<-add(AdjM[tri[3],,drop=F],last[tri[3]],tri[2])
            AdjM[tri[3],]<-tmp$v
            last[tri[3]]<-tmp$last
            rm(tmp)
        }
        
    }
    
    # Output
    
    out<-NULL
    out$AdjM<-AdjM
    out$cn<-cn
        
    return(out)
    
}