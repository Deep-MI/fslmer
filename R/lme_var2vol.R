lme_var2vol<-function(data,v=1,type=3,dof=1)
{
    
    if (!is.array(data)) data<-as.array(data)
    
    if (length(dim(data))<4) ndim4<-1 else ndim4<-dim(data)[4]
    if (length(dim(data))<3) ndim3<-1 else ndim3<-dim(data)[3]
    if (length(dim(data))<2) ndim2<-1 else ndim2<-dim(data)[2]
    ndim1<-dim(data)[1]
    
    data<-array(data,dim=c(ndim1,ndim2,ndim3,ndim4))
    
    out<-list(x=data,v=v,ndim1=ndim1,ndim2=ndim2,ndim3=ndim3,nframes=ndim4,type=type,dof=dof)
              
    return(out)
}


