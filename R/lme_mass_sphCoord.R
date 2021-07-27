lme_mass_sphCoord<-function(Coord3D)
{
    
    X<-Coord3D[,1,drop=F]
    Y<-Coord3D[,2,drop=F]
    Z<-Coord3D[,3,drop=F]

    phi<-atan2(Y,X)
    theta<-atan2(sqrt(X^2+Y^2),Z)
    
    # Output
    
    out<-NULL
    out$phi<-phi
    out$theta<-theta
    
    return(out)
    
}
    
    
