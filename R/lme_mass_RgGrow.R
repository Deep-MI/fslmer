#' Region-growing algorithm to identify spatially homogeneous regions
#'
#' @param SphSurf Spherical surface
#' @param Re Output from \code{lme_mass_fit_init}
#' @param Theta Output from \code{lme_mass_fit_init}
#' @param maskvtx Indices (one-based) for vertices to include in the analysis (default: all)
#' @param nst Internal parameter (default: 2)
#' @param prc Internal parameter (default: 95)
#'
#' @return
#' A list with spatially homogeneous regions and their means.
#' 
#' @export
#'
#' @examples
#' RgGrow <- lme_mass_RgGrow(SphSurf, Re, Theta, maskvtx=NA, nst=2, prc=95)

lme_mass_RgGrow<-function(SphSurf,Re,Theta,maskvtx=NA,nst=2,prc=95)
{
    # -------------------------------------------------------------------------
    # Auxiliary functions

    # vtxwGrowing

    vtxwGrowing<-function(Re,Params,coord,Adj,maskvtx,nst,prc)
    {
        # seeds

        print('Computing seeds ...',quote=F)
        outSeeds<-seeds(Re[,maskvtx,drop=F],Params[,maskvtx,drop=F],coord[,maskvtx,drop=F],nst,prc)
        Rgs1<-outSeeds$Rgs
        Rgseed<-outSeeds$Rgseed
        Rgstd<-outSeeds$Rgstd
        Rgseed<-maskvtx[Rgseed]

        ns<-length(Rgseed)
        print(paste(ns,'seeds were computed after splitting the surface'))

        np<-nrow(Params)
        nv<-ncol(Params)

        Rgs<-matrix(0,1,nv)
        Rgs[maskvtx]<-Rgs1
        notIncl<-matrix(TRUE,1,nv)
        notIncl[Rgs==0]<-FALSE
        notIncl[Rgseed]<-FALSE
        maxRgsz<-90 # Do not form regions with size greater than maxRgsz
        Rgvtxs<-matrix(0,ns,maxRgsz)
        Rgvtxs[,1]<-Rgseed
        Rgvtxs_prvtxind<-matrix(0,1,ns)
        Rgvtxs_lastind<-matrix(1,1,ns)
        Rgind<-c(1:ns)
        thr<-nst*Rgstd
        growing_ns<-ns

        while (growing_ns>0)
        {

            print(paste(growing_ns,'seeds for growing'),quote=F)
            print(paste('Current maximum region size',max(Rgvtxs_lastind),'vertices'))

            i<-1
            while (i<=growing_ns)
            {

                Rgvtxs_prvtxind[Rgind[i]]<-Rgvtxs_prvtxind[Rgind[i]]+1
                prvtx<-Rgvtxs[Rgind[i],Rgvtxs_prvtxind[Rgind[i]]]
                RgParams<-Params[,Rgvtxs[Rgind[i],c(1:Rgvtxs_lastind[Rgind[i]]),drop=F],drop=F]
                mRgParams<-rowMeans(RgParams)
                Rgsz<-Rgvtxs_lastind[Rgind[i]]

                # processing the current vertex prvtx
                adjvtxs<-Adj[prvtx,Adj[prvtx,,drop=F]>0]

                # try to include only the adjacent vertices not included in other regions or in this region itself
                adjvtxs<-adjvtxs[notIncl[adjvtxs]]

                # try to include first the adjacent vertices with the param closest to the region mean param
                loc<-which.min(colSums(abs(Params[,adjvtxs] - kronecker(matrix(1,1,length(adjvtxs)),mRgParams))))
                NewRgParams<-cbind(RgParams,Params[,adjvtxs[loc]])
                mNewRgParams<-rowMeans(NewRgParams)
                nRgnv<-ncol(NewRgParams)
                DistNewRgParams<-abs(NewRgParams-kronecker(matrix(1,1,nRgnv),mNewRgParams))

                # an adjacent vertex is included in the region if the difference
                # between its params and the region mean params is bellow the
                # threshold thr or if the region still contains a percent
                # prc of vertices bellow the threshold.

                nadjv <- length(adjvtxs)
                j <- 1

                while ((Rgsz<maxRgsz) && (j<=nadjv) && sum(colSums(DistNewRgParams<=kronecker(matrix(1,1,nRgnv),thr[,Rgind[i],drop=F])) == np) >= prc*(nRgnv)/100)
                {
                    # In addition the vertex's residuals must be correlated above 0.5
                    # with any other vertex's residual in the region.

                    if (min(min(stats::cor(cbind(Re[,Rgvtxs[Rgind[i],1:Rgvtxs_lastind[Rgind[i]]]],Re[,adjvtxs[loc],drop=F]))))>= 0.5)
                    {
                        Rgvtxs_lastind[Rgind[i]]<-Rgvtxs_lastind[Rgind[i]]+1
                        Rgsz<-Rgsz+1
                        Rgvtxs[Rgind[i],Rgvtxs_lastind[Rgind[i]]]<-adjvtxs[loc]
                        notIncl[adjvtxs[loc]]<-FALSE
                        RgParams<-Params[,Rgvtxs[Rgind[i],c(1:Rgvtxs_lastind[Rgind[i]])],drop=F]
                        mRgParams<-rowMeans(RgParams)
                    }

                    # adjvtxs<-c(adjvtxs[c(1:(loc-1))],adjvtxs[c((loc+1):length(adjvtxs))]) # KD
                    # changed since indexing in R is different from Matlab, KD
                    if ((((loc-1)>0)) & ((loc+1)<=length(adjvtxs)))
                    {
                        adjvtxs<-c(adjvtxs[c(1:(loc-1))],adjvtxs[c((loc+1):length(adjvtxs))])
                    }
                    else
                    {
                        if ((((loc-1)>0)) & (!((loc+1)<=length(adjvtxs))))
                        {
                            adjvtxs<-adjvtxs[c(1:(loc-1))]
                        }
                        else
                        {
                            if ((!(((loc-1)>0))) & ((loc+1)<=length(adjvtxs)))
                            {
                                adjvtxs<-adjvtxs[c((loc+1):length(adjvtxs))]
                            }
                            else
                            {
                                adjvtxs<-NULL
                            }
                        }
                    }
                    loc<-which.min(colSums(abs(Params[,adjvtxs] - kronecker(matrix(1,1,length(adjvtxs)),mRgParams))))
                    NewRgParams<-cbind(RgParams,Params[,adjvtxs[loc],drop=F])
                    mNewRgParams<-rowMeans(NewRgParams)
                    nRgnv<-ncol(NewRgParams)
                    DistNewRgParams<-abs(NewRgParams-kronecker(matrix(1,1,nRgnv),mNewRgParams))
                    j<-j+1

                }

                #

                if ((Rgvtxs_prvtxind[Rgind[i]]==Rgvtxs_lastind[Rgind[i]]) || (Rgvtxs_lastind[Rgind[i]]==maxRgsz))
                {
                    # Rgind<-c(Rgind[1:(i-1)],Rgind[(i+1):length(Rgind)]) # changed because indexing in R is different from Matlab, KD
                    if ((((i-1)>0)) & ((i+1)<=length(Rgind)))
                    {
                        Rgind<-c(Rgind[c(1:(i-1))],Rgind[c((i+1):length(Rgind))])
                    }
                    else
                    {
                        if ((((i-1)>0)) & (!((i+1)<=length(Rgind))))
                        {
                            Rgind<-Rgind[c(1:(i-1))]
                        }
                        else
                        {
                            if ((!(((i-1)>0))) & ((i+1)<=length(Rgind)))
                            {
                                Rgind<-Rgind[c((i+1):length(Rgind))]
                            }
                            else
                            {
                                Rgind<-NULL
                            }
                        }
                    }
                    growing_ns<-growing_ns-1
                }
                else
                {
                    i<-i+1
                }

            }

        }

        for (i in c(1:ns))
        {
            Rgs[Rgvtxs[i,1:Rgvtxs_lastind[i]]]<-i
        }

        nni<-sum(notIncl)
        notInclv<-which(notIncl==1)
        Rgs[notInclv]<-c((ns+1):(ns+nni))

        # try to add holes to their most similar adjacent homogeneous region
        print(paste(nni,'holes (unassigned vertices). Trying to form new regions among them or add them to their most similar adjacent region ...'),quote=F)

        for (i in c(1:nni))
        {

            hadjvtxs<-Adj[notInclv[i],Adj[notInclv[i],,drop=F]>0]

            if (!(length(hadjvtxs)==0))
            {
                closestdistRParams<-sum(abs(Params[,notInclv[i],drop=F]-rowMeans(Params[,Rgs==Rgs[hadjvtxs[1]],drop=F])))
                nadjvtxs<-length(hadjvtxs)
                m<-1; j<-2
                while (j<=nadjvtxs)
                {
                    distRParams<-sum(abs(Params[,notInclv[i],drop=F]-rowMeans(Params[,Rgs==Rgs[hadjvtxs[j]],drop=F])))
                    if (distRParams<closestdistRParams)
                    {
                        closestdistRParams<-distRParams
                        m<-j
                    }
                    j<-j+1
                }

                outMerge<-merge(Rgs,Re,Params,Rgs[notInclv[i]],Rgs[hadjvtxs[m]],nst,prc)
                Rgs<-outMerge$Rgs
                mrg<-outMerge$tf

                if (mrg) notIncl[notInclv[i]]<-0

            }

        }

        nni<-sum(notIncl)
        if (nni>0)
        {
            print(paste(nni,'holes (unassigned vertices) remained. Will be treated as independent regions of size 1.'),quote=F)
        }
        else
        {
            print('No holes (unassigned vertices) remained.',quote=F)
        }

        # Output

        return(Rgs)

    }

    #

    # seeds

    seeds<-function(Re,Params,coord,nst,prc)
    {

        np<-nrow(Params)
        nv<-ncol(Params)

        Rgs<-matrix(1,1,nv)
        spl<-TRUE

        while (spl)
        {

            # splRgs<-unique(Rgs)
            splRgs<-unique(sort(Rgs)) # KD
            ntospl<-length(splRgs)
            spltf<-matrix(0,1,ntospl)
            nRgs<-ntospl

            for (i in c(1:ntospl))
            {
                outSplit<-split(Rgs,Re,Params,coord,splRgs[i],nst,prc)

                Rgs<-outSplit$Rgs
                spltf[i]<-outSplit$tf
                nnRgs<-outSplit$nRnv

                nRgs<-nRgs+nnRgs
            }
            if (sum(spltf)==0) spl<-FALSE

        }

        # splRgs<-unique(Rgs)
        splRgs<-unique(sort(Rgs)) # KD
        nRgs<-length(splRgs)
        Rgseed<-matrix(0,1,nRgs)
        Rgstd<-matrix(0,np,nRgs)

        for (i in c(1:nRgs))
        {
            Rgv<-which(Rgs==splRgs[i])
            RgParams<-Params[,Rgv,drop=F]
            Rgstd[,i]<-sqrt(diag(stats::var(t(RgParams))))
            Rgstd[is.na(Rgstd)]<-0 # since var(0)==NA in R, unlike Matlab, where var(0)==0, KD
            nv<-length(Rgv)
            spp<-rowSums(coord[,Rgv,drop=F])/nv
            spp<-round(spp,6)
            loc<-which.min(colSums(abs(coord[,Rgv] - kronecker(matrix(1,1,nv),spp))))
            Rgseed[i]<-Rgv[loc]
        }

        # Output

        out<-NULL
        out$Rgs<-Rgs
        out$Rgseed<-Rgseed
        out$Rgstd<-Rgstd

        return(out)

    }

    #

    # split

    split<-function(Rgs,Re,Params,coord,Rnum,nst,prc)
    {

        if (homogenity(Re,Params,Rgs,Rnum,nst,prc))
        {
            nRnv<-0
            tf<-FALSE
        }
        else
        {
            # split the non-homogeneous region
            Rvts<-which(Rgs==Rnum)
            nv<-length(Rvts)
            spp<-rowSums(coord[,Rvts,drop=F])/nv
            spp<-round(spp,6)
            nR<-max(Rgs)
            for (i in c(1:nv))
            {
                if ((coord[1,Rvts[i]]<=spp[1]) && (coord[2,Rvts[i]]<=spp[2]))
                {
                    Rgs[Rvts[i]]<-1
                }
                else
                {
                    if ((coord[1,Rvts[i]]<=spp[1]) && (coord[2,Rvts[i]]>spp[2]))
                    {
                        Rgs[Rvts[i]]<-2
                    }
                    else
                    {
                        if ((coord[1,Rvts[i]]>spp[1]) && (coord[2,Rvts[i]]<=spp[2]))
                        {
                            Rgs[Rvts[i]]<-3
                        }
                        else
                        {
                            if ((coord[1,Rvts[i]]>spp[1]) && (coord[2,Rvts[i]]>spp[2]))
                            {
                                Rgs[Rvts[i]]<-4
                            }
                        }
                    }
                }
            }

            # Rnv<-unique(Rgs[Rvts])
            Rnv<-unique(sort(Rgs[Rvts])) # KD
            nRnv<-length(Rnv)-1

            for (i in c(1:nRnv)) Rgs[Rvts[Rgs[Rvts]==Rnv[i+1]]]<-nR+i

            Rgs[Rvts[Rgs[Rvts]==Rnv[1]]]<-Rnum
            mrgRgs<-c(Rnum,(nR+1):(nR+nRnv))

            # merge homogeneous subregions
            nmtomrg<-length(mrgRgs)
            i<-1
            while (i<=nmtomrg)
            {
                j<-i+1
                mrg<-FALSE
                while ((!mrg) && (j<=nmtomrg))
                {
                    outMerge<-merge(Rgs,Re,Params,mrgRgs[i],mrgRgs[j],nst,prc)
                    Rgs<-outMerge$Rgs
                    mrg<-outMerge$tf
                    j<-j+1
                }
                if (mrg)
                {
                    mrgRgs<-c(mrgRgs[1:(j-2)],mrgRgs[j:length(mrgRgs)])
                    nmtomrg<-nmtomrg-1
                }
                else
                {
                    i<-i+1
                }

            }

            tf<-TRUE

        }

        # Output

        out<-NULL
        out$Rgs<-Rgs
        out$tf<-tf
        out$nRnv<-nRnv

        return(out)

    }

    #

    # merge

    merge<-function(Rgs,Re,Params,Rn1,Rn2,nst,prc)
    {

        if (Rn1<=Rn2)
        {
            Rkn<-Rn1
            Ren<-Rn2
        }
        else
        {
            Rkn<-Rn2
            Ren<-Rn1
        }

        auxRgs<-Rgs
        auxRgs[auxRgs==Ren]<-Rkn

        if (homogenity(Re,Params,auxRgs,Rkn,nst,prc))
        {
            Rgs<-auxRgs
            tf<-TRUE
        }
        else
        {
            tf<-FALSE
        }

        # Output

        out<-NULL
        out$Rgs<-Rgs
        out$tf<-tf

        return(out)

    }

    #

    # homogeneity

    homogenity<-function(Re,Params,Rgs,Rn,nst,prc)
    {

        tf<-FALSE
        maxRgsz<-90

        # Regions with size greater than maxRgsz are not considered homogeneous

        Rgv<-which(Rgs==Rn)

        # if (length(Rgv)<=maxRgsz)
        if ((length(Rgv)>0) && (length(Rgv)<=maxRgsz)) # KD
        {
            RgParams<-Params[,Rgv,drop=F]
            mRgParams<-rowMeans(RgParams)
            stRgParams<-sqrt(diag(stats::var(t(RgParams))))
            stRgParams[is.na(stRgParams)]<-0 # since var(0)==NA in R, unlike Matlab, where var(0)==0, KD
            thr<-nst*stRgParams
            np<-nrow(RgParams)
            nv<-ncol(RgParams)
            DistRgParams = abs(RgParams-kronecker(matrix(1,1,nv),mRgParams))

            tmpcor<-stats::cor(Re[,Rgv,drop=F]); if (length(tmpcor)>1) {diag(tmpcor)<-NA}; if (all(is.na(tmpcor))) {tmpcor<-FALSE} else {tmpcor<-(min(tmpcor,na.rm=T) >= 0.5)} # KD
            # if ((sum(colSums(DistRgParams <= kronecker(matrix(1,1,nv),thr)) == np) >= prc*nv/100) && min(stats::cor(Re[,Rgv,drop=F]),na.rm=T) >= 0.5) # KD
            if ((sum(colSums(DistRgParams <= kronecker(matrix(1,1,nv),thr)) == np) >= prc*nv/100) && tmpcor) # KD
            {
                tf<-TRUE
            }
            rm(tmpcor)

        }

        # Output

        return(tf)

    }

    #

    # -------------------------------------------------------------------------
    # Main function

    nv<-ncol(Theta)

    if (all(is.na(maskvtx))) maskvtx<-c(1:nv)

    outAdjM<-lme_mass_AdjMtx(SphSurf,maskvtx)
    AdjM<-outAdjM$AdjM

    outSphC<-lme_mass_sphCoord(t(SphSurf$vertices))
    sphcoord<-t(cbind(outSphC$phi,outSphC$theta))

    Regions<-vtxwGrowing(Re,Theta,sphcoord,AdjM,maskvtx,nst,prc)

    outRgMean<-lme_mass_RgMean(Regions,Theta)

    # output

    out<-NULL
    out$Regions<-Regions
    out$RgMeans<-outRgMean$RgMeans

    return(out)

}

#
