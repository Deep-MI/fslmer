#' 	Mass-univariate linear mixed-effects model fitting using spatially homogeneous regions
#'
#' @param X Design matrix
#' @param Zcols Vector of random effect indices in design matrix
#' @param Y Outcome variable
#' @param ni Vector indicating the repeated observations of each subject
#' @param Th0 Output from \code{lme_mass_fit_init}
#' @param Rgs Output from \code{lme_mass_RgGrow}
#' @param Surf Spherical surface from \code{lme_readsurf}
#' @param fname (currently ignored)
#' @param Dtype Distance (default: euclidean)
#' @param sptm Spatial measure (default: exp)
#' @param prs Number of cores for parallel computing (default: 1)
#' @param e Tolerance (default: 10^-1)
#'
#' @return
#' This function returns a list of lists, with entries st (=status) and stats,
#' with the following entries: Bhat, CovBhat, phisqhat, Dhat, Zcols, invEI,
#' Pth, Qthth, lreml
#'
#' @export
#'
#' @examples
#' \dontrun{Surf <- lme_readsurf(...)}
#' \dontrun{fitInit <- lme_mass_fit_init(...)}
#' \dontrun{RgGrow <- lme_mass_RgGrow(...)}
#' \dontrun{fitRgw <- lme_mass_fit_Rgw(X, Zcols, Y, ni, fitInit$Theta0, RgGrow$Rgs, Surf)}

lme_mass_fit_Rgw<-function(X,Zcols,Y,ni,Th0,Rgs,Surf,fname=NA,Dtype="euc",sptm="exp",prs=1,e=0.1)
{

    # --------------------------------------------------------------------------
    # check if parallel computing is feasible

    if (prs==1) print("No parallel computing enabled (not recommended)",quote=F)

    # --------------------------------------------------------------------------
    # Auxiliary functions

    # convergence

    convergence<-function(Rgst,lreml,i,i_N=NA)
    {
        tf<-FALSE

        if (Rgst  && (lreml[length(lreml)] >= lreml[1]))
        {
            print(paste('Region ', i, '/', i_N, ': Convergence at iteration ', length(lreml), '. Initial and final likelihoods: ', lreml[1], ', ', lreml[length(lreml)], '.',sep=""),quote=F)
            tf<-TRUE
        }
        else
        {
            if (Rgst && (lreml[length(lreml)] < lreml[1]))
            {
                print(paste('Region ', i, '/', i_N, ': Convergence to a saddle point at iteration ', length(lreml), '. Initial and final likelihoods: ', lreml[1], ', ', lreml[length(lreml)], '.',sep=""),quote=F)
            }
            else
            {
                print(paste('Region ', i, '/', i_N, ': Algorithm did not converge. Initial and final likelihoods: ', lreml[1], ', ', lreml[length(lreml)], '.', sep=""),quote=F)
            }
        }

        return(tf)
    }

    # balanceRgind

    balanceRgind<-function(prs,nRg)
    {
        bRgind <- matrix(0,prs,ceiling(nRg/prs))
        prsnRg <- matrix(0,prs,1)
        reverse <- FALSE

        i <- 1
        while (i <= nRg)
        {
            if (reverse)
            {
                j <- prs
                while ((i <= nRg) && (j >=1))
                {
                    prsnRg[j] = prsnRg[j] + 1
                    bRgind[j,prsnRg[j]] <- i
                    j <- j - 1
                    i <- i + 1
                }
                reverse <- FALSE
            }
            else
            {
                j <- 1
                while ((i <= nRg) && (j <= prs))
                {
                    prsnRg[j] <- prsnRg[j] + 1
                    bRgind[j,prsnRg[j]] <- i
                    j <- j + 1
                    i <- i + 1
                }
                reverse <- TRUE
            }
        }

        #

        out<-NULL
        out$bRgind <- bRgind
        out$prsnRg <- prsnRg
        return(out)

    }

    # parallel processing

    do_estimate<-function(j,prsnRg,Rgs,bRgind,RgVtxInd,nRgvtx,Surf,maskvtx,Dtype,X,Zcols,Yj,Th0j,Dist,sptm,ni,e)
    {

        progress<-0
        posi<-0

        for (i in c(1:prsnRg[j]))
        {
            RgVtxInd <- which(Rgs == bRgind[j,i])
            nRgvtx <- length(RgVtxInd)
            posf <- posi+nRgvtx

            tryCatch(
            {
                if (nRgvtx > 1)
                {
                    Dist <- lme_mass_RgDist(Surf,maskvtx,RgVtxInd,Dtype)
                    outRgFSfit<-lme_RgFSfit(X,Zcols,Yj[[j]][,(posi+1):posf,drop=F],ni,Th0j[[j]][,(posi+1):posf],Dist,sptm)

                    Rgstats[(posi+1):posf]<-outRgFSfit$stats
                    Rgsts[(posi+1):posf]<-outRgFSfit$st
                }
                else
                {
                    outFSfit <- lme_FSfit(X,Zcols,Yj[[j]][,posf,drop=F],ni,e)
                    Rgstats[[posf]]<-outFSfit[c(1:4,6:10)]
                    Rgsts[posf]<-outFSfit$st
                }

                if (convergence(Rgsts[[posf]],Rgstats[[posf]]$lreml,bRgind[j,i],length(unique(sort(Rgs))))) sRgst <- sRgst + 1 # moved this up
            },
                 error=function(e){print(paste("Region ", bRgind[j,i], ': did not run',e,sep=""),quote=F)}
            )

            posi <- posf

        }

        # output

        out<-NULL

        out$Rgstats<-Rgstats
        out$Rgsts<-Rgsts

        return(out)
    }

    # --------------------------------------------------------------------------
    # Main function

    n<-sum(ni)
    nv0<-ncol(Y)

    maskvtx<-which(Rgs!=0)
    nv<-length(maskvtx)
    Rgs<-Rgs[maskvtx]

    Y <- Y[,maskvtx,drop=F]
    Th0<-Th0[,maskvtx,drop=F]

    Rgnums<-unique(sort(Rgs))
    nRg<-length(Rgnums)
    Rglth <- matrix(0,1,nRg)

    for (i in c(1:nRg)) Rglth[i] <- sum(Rgs == Rgnums[i])

    # Sort regions by length in descend order
    ix<-order(Rglth,decreasing=T)
    aux<-Rgs
    for (i in c(1:nRg)) Rgs[aux == Rgnums[ix[i]]]<-i

    # Balance region size across workers
    if (prs>nRg) {prs<-nRg; print(paste0("Changed number of cores to" ,prs," since there are only ",nRg," regions to work on."))}

    outBal<-balanceRgind(prs,nRg)
    bRgind<-outBal$bRgind
    prsnRg<-outBal$prsnRg

    # Initialization

    statstruct <- list(Bhat=NA,CovBhat=NA,phisqhat=NA,Dhat=NA,Zcols=NA,invEI=NA,Pth=NA,Qthth=NA,lreml=-10^10)

    prsnumloc<-matrix(0,prs,1)

    for (j in c(1:prs))
    {
        for (i in c(1:prsnRg[j]))
        {
            nRgvtx <- sum(Rgs == bRgind[j,i])
            prsnumloc[j] <- prsnumloc[j] + nRgvtx
        }
    }

    nrf<-nrow(Th0)

    Rgstats<-NULL
    Rgsts<-NULL
    Yj<-NULL
    Th0j<-NULL

    for (j in c(1:prs))
    {
        Rgstats[[j]]<-vector("list",prsnumloc[j])
        for (k in c(1:prsnumloc[j])) Rgstats[[j]][[k]] <- statstruct
        Rgsts[j] <- list(matrix(FALSE,prsnumloc[j],1))
        Yj[j] <- list(matrix(0,n,prsnumloc[j]))
        Th0j[j] <- list(matrix(0,nrf,prsnumloc[j]))
    }

    for (j in c(1:prs))
    {
        posi <- 0
        for (i in c(1:prsnRg[j]))
        {
            RgVtxInd <- which(Rgs == bRgind[j,i])
            nRgvtx <- length(RgVtxInd)
            posf <- posi+nRgvtx
            if (nRgvtx)
            {
                Yj[[j]][,(posi+1):(posf)]<-Y[,RgVtxInd]
                Th0j[[j]][,(posi+1):(posf)]<-Th0[,RgVtxInd]
            }
            else
            {
                print(paste("Found region of size 0:",i,bRgind[j,i],posi+1,posf,length(RgVtxInd))) # KD
            }
            posi <- posf
        }
    }

    rm(Y,Th0)

    # Estimation

    print('Starting model fitting at each region ...',quote=F)

    sRgst <- 0

    if (prs>1)
    {
        outEst<-parallel::mclapply(c(1:prs),function(x) do_estimate(x,prsnRg,Rgs,bRgind,RgVtxInd,nRgvtx,Surf,maskvtx,Dtype,X,Zcols,Yj,Th0j,Dist,sptm,ni,e),mc.cores=prs)

        for (j in c(1:prs))
        {
            Rgstats[[j]]<-outEst[[j]]$Rgstats
            Rgsts[[j]]<-outEst[[j]]$Rgsts
        }

        rm(outEst)
    }
    else
    {
        for (j in c(1:prs))
        {
            outEst<-do_estimate(j,prsnRg,Rgs,bRgind,RgVtxInd,nRgvtx,Surf,maskvtx,Dtype,X,Zcols,Yj,Th0j,Dist,sptm,ni,e)

            Rgstats[[j]]<-outEst$Rgstats
            Rgsts[[j]]<-outEst$Rgsts

            rm(outEst)
        }
    }

    #

    stats<-NULL
    for (i in c(1:nv0)) stats[[i]] <- statstruct

    st <- matrix(NA,nv0,1)

    for (j in c(1:prs))
    {
        posi <- 0
        for (i in c(1:prsnRg[j]))
        {
            RgVtxInd <- which(Rgs == bRgind[j,i])
            nRgvtx <- length(RgVtxInd)
            posf <- posi+nRgvtx
            Rgstatsji <- Rgstats[[j]][(posi+1):posf]
            Rgstsji <- Rgsts[[j]][(posi+1):posf]

            for (k in c(1:nRgvtx))
            {
                stats[maskvtx[RgVtxInd[k]]] <- Rgstatsji[k]
                stats[[maskvtx[RgVtxInd[k]]]]$lreml <- stats[[maskvtx[RgVtxInd[k]]]]$lreml[length(stats[[maskvtx[RgVtxInd[k]]]]$lreml)]
                st[maskvtx[RgVtxInd[k]]] <- Rgstsji[k]
            }
            posi <- posf

        }
    }

    # Output

    out<-NULL
    out$stats<-stats
    out$st<-st

    return(out)

}
