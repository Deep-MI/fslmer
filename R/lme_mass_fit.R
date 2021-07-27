lme_mass_fit<-function(X,Xcols=NA,Xrows=NA,Zcols,Y,ni,prs=1,e=10^-1,numcore)
{    
    # Initialize
    
    out<-NULL
    
    # Validation of inputs

    nm<-nrow(Y)
    nv<-ncol(Y)

    p<-ncol(X)
    
    if (is.na(Xcols)) 
    {
        Xcols<-matrix(1,nv,p)
    }

    if (nrow(Zcols)==1)
    {
        rfcols<-which(Zcols==1)
        Zcols<-matrix(0,nv,p)
        Zcols[,rfcols]<-matrix(1,nv,length(rfcols))
    }

    #
    
    if (!is.na(Xrows))
    {
        sID<-matrix(0,nm,1)
        m<-length(ni)
        count<-0
        for (i in c(1:m))
        {
            sID[(count+1):(count+ni[i])]<-i
            count<-count+ni[i]
        }
        
        #
        
        build_the_model<-function(i,Y,X,Xrows,Xcols,Zcols,ni,e)
        {
            # Build the model for this specific location
            rows<-which(Xrows[,i]==1)
            XM<-X[rows,Xcols[i,]==1]
            sIDM<-sID[rows]
			mM<-length(unique(sort(sIDM)))
            niM<-matrix(1,mM,1)
            
            count<-1
            for (j in c(2:length(rows)))
            {
                if (sIDM[j-1]==sIDM[j])
                {
                    niM[count]<-niM[count]+1
                }
                else
                {
                    count<-count+1
                }
            }
            
            rfcols<-which(Zcols[i,]==1)
            
            # Estimation
            stat<-NULL
            stat<-lme_FSfit(XM,rfcols,Y[,i,drop=F],ni,e)
            stat$lreml<-try(stat$lreml[length(stat$lreml)])
            
            # Print
            
			if ((i%%100)==0) {print(paste("Working on vertex ",i))}
            
            # Output
            
            return(stat)
            
        }
        
        #
        
        if (numcore==1)
        {      
            for (i in c(1:nv))
            {
                out[[i]]<-build_the_model(i,Y,X,Xrows,Xcols,Zcols,ni,e)
            }
        }
        else
        {
            if (numcore>1)
            {
                out<-parallel::mclapply(c(1:nv),build_the_model,Y,X,Xrows,Xcols,Zcols,ni,e,mc.cores=numcore)
            }
        }
    
    }
    else
    {
        
        build_the_model<-function(i,Y,X,Xcols,Zcols,ni,e)
        {
            # Build the model for this specific location
            XM<-X[,Xcols[i,]== 1]
            rfcols<-which(Zcols[i,]==1)
            
            # Estimation
            stat<-NULL
            stat<-lme_FSfit(XM,rfcols,Y[,i,drop=F],ni,e)
            stat$lreml<-try(stat$lreml[length(stat$lreml)])
            
            # Print
            if ((i%%100)==0) {print(paste("Working on vertex ",i))}
            
            # Output
            return(stat)
        }
        
        if (numcore==1)
        {      
            for (i in c(1:nv))
            {
                out[[i]]<-build_the_model(i,Y,X,Xcols,Zcols,ni,e)
            }
        }
        else
        {
            if (numcore>1)
            {
            out<-parallel::mclapply(c(1:nv),build_the_model,Y,X,Xcols,Zcols,ni,e,mc.cores=numcore)
            }
        }
    }
    
    #
    
    return(out)
    
}
