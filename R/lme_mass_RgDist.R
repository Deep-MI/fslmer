lme_mass_RgDist<-function(Surf,maskvtx,vtxInd,Dtype="euc")
{

    nv<-length(vtxInd)
    
    if (nv== 1)
    {
        Dist<-0
    }
    else
    {
        if (Dtype=='euc')
        {
            Coord<-Surf$vertices[,maskvtx[vtxInd],drop=F]
            Dist<-matrix(0,nv,nv)
            aux<-matrix(1,1,nv)
            for (i in c(1:nv))
            {
                c<-Coord[,i,drop=F]
                Dist[i,] <- sqrt(colSums((Coord - kronecker(aux,c))^2)) 
            }
        }
        else
        {
            if (Dtype=='surf')
            {    
                stop("The geodesic algorithm is currently not implemented.")
            }
        }
            
    }
    
    return(Dist)
    
}