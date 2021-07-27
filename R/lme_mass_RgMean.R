lme_mass_RgMean<-function(Rgs,Data)
{
    
    npar<-nrow(Data)
    nv<-ncol(Data)
    
    RgMeans<-matrix(0,npar,nv)
    Rgnums<-unique(sort(Rgs))
    nRg<-length(Rgnums)
    
    for (i in c(1:nRg))
    {
        Rgvtxs<-which(Rgs==Rgnums[i])
        mRgData<-rowMeans(Data[,Rgvtxs,drop=F])
        RgMeans[,Rgvtxs] <- kronecker(matrix(1,1,length(Rgvtxs)),mRgData) 
    }

    # Output

    out<-NULL
    out$RgMeans<-RgMeans
    out$nRg<-nRg
    
    return(out)
        
}