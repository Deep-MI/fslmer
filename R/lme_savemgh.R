lme_savemgh<-function(vol, fname)
{
    # R translation of save_mgh.m

    MRI.UCHAR<-0
    MRI.INT<-1
    MRI.LONG<-2
    MRI.FLOAT<-3
    MRI.SHORT<-4
    MRI.BITMAP<-5
    MRI.TENSOR<-6

    fid<-file(fname,open = "wb", blocking = TRUE)

    width<-vol$ndim1
    height<-vol$ndim2
    depth<-vol$ndim3
    nframes<-vol$nframes

    writeBin(as.integer(1),fid,size = 4, endian = "big")            # magic

    writeBin(as.integer(width),fid,size = 4, endian = "big")        # x
    writeBin(as.integer(height),fid,size = 4, endian = "big")       # y
    writeBin(as.integer(depth),fid,size = 4, endian = "big")        # z
    writeBin(as.integer(nframes),fid,size = 4, endian = "big")      # frames

    writeBin(as.integer(MRI.FLOAT),fid,size = 4, endian = "big")    # type (currently only MRI.FLOAT implemented)

    writeBin(as.integer(0),fid,size = 4, endian = "big")            # dof (default: 1)

    #

    UNUSED.SPACE.SIZE<-256
    USED.SPACE.SIZE<-1*2                                            # space for flag only
    #USED.SPACE.SIZE<-(3*4+4*3*4+1*2)                               # space for ras transform plus flag

    unused.space.size<-UNUSED.SPACE.SIZE-USED.SPACE.SIZE

    writeBin(as.integer(0),fid,size = 2, endian = "big")            # ras-good-flag; change to 1 if RAS transform is implemented

    # fill up unused sapce

    writeBin(as.integer(rep(0,unused.space.size)),fid,size = 1)

    # now we'll write the data

    for (k in c(1:nframes))
    {
        for (j in c(1:depth))
        {
            for (i in c(1:height))
            {
                writeBin(vol$x[,i,j,k],fid, size = 4, endian = "big")
            }
        }
    }

    # that's all

    close(fid)

}
