#' Title
#'
#' @param fname
#'
#' @return
#' @export
#'
#' @examples

lme_readsurf<-function(fname)
{

    TRIANGLE_FILE_MAGIC_NUMBER<-16777214
    QUAD_FILE_MAGIC_NUMBER<-16777215

    fid<-file(fname,"rb")

    magic<-strtoi(paste("0x",paste(readBin(fid,raw(),n=3,endian="big"),collapse=""),sep=""))

    if(magic == QUAD_FILE_MAGIC_NUMBER)
    {
        stop("reading of QUAD files is currently not supported")

        # this is how the matlab code looks like
        # vnum = fread3(fid) ;
        # fnum = fread3(fid) ;
        # vertex_coords = fread(fid, vnum*3, 'int16') ./ 100 ;
        #
        #if (nargout > 1)
        #    for i=1:fnum
        #       for n=1:4
        #            faces(i,n) = fread3(fid) ;
        #        end
        #    end
        #end
    }
    else
    {
        if (magic == TRIANGLE_FILE_MAGIC_NUMBER)
        {
            readLines(fid,n=1)
            readLines(fid,n=1)

            vnum<-readBin(fid,integer(),size=4,endian="big")
            fnum<-readBin(fid,integer(),size=4,endian="big")
            vertex_coords<-readBin(fid,numeric(),size=4,n=vnum*3,endian="big")
            faces<-readBin(fid,integer(),size=4,endian="big",n=fnum*3)

        }

    }

    close(fid)

    vertex_coords<-matrix(vertex_coords,3)
    faces<-t(matrix(faces,3))+1

    # output

    out<-NULL
    out$vertices<-vertex_coords
    out$faces<-faces

    return(out)
}
