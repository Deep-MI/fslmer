lme_readlabel<-function(fname)
{
    # reads the label file 'fname' 
    # output will be nvertices-by-5, where each column means:
    # (1) vertex number, (2-4) xyz at each vertex, (5) stat
    #
    # IMPORTANT: the vertex number is 0-based.
    
    fid<-file(fname,'r')
    
    readLines(fid,n=1)
    
    nv<-scan(fid,what=integer(),n=1)
    
    d<-simplify2array(scan(fid,what=list(integer(),numeric(),numeric(),numeric(),numeric()),n=5*nv))
    
    close(fid)
    
    return(d)
}