#' Estimation of subject-specific random effects estimates at each vertex
#'
#' @param stats Structure array containing statistics for every voxel/vertex (generated with either lme_mass_fit_Rgw or lme_mass_fit_vw)
#' @param X Ordered design matrix (according to time for each subject)
#' @param Zcols Vector with the indices of the colums of X that are considered as random effects
#' @param Y Ordered data matrix (n x nv, n=total number of scans and nv=number of vertices/voxels)
#' @param ni Vector whose m entries are the number of repeated measures for each subject (ordered according to X)
#' @param maskvtx Mask's vertices (1-based). Default: NA (all vertices included)
#' @param prs Number of cores for parallel computing (default: 1)
#'
#' @return
#' This function returns the subject-specific random effects estimates at each vertex. The output is a list of lists, with the following entries:
#' Rfx: Estimated subject-specific random effects matrix (m x nrfx*nv). The columns of this matrix are grouped by vertex. For example if
#' there are two random effects in the model then the first two columns contain the subject-specific random effect coefficients for the first vertex,
#' then the next two columns contain the subject-specific random effect coefficients for the second vertex and so on ...
#' nrfx: Number of random effects (length(Zcols)).
#' Bhat: Population-level regression coefficients in stats stacked in one matrix.
#'
#' @export
#'
#' @examples
#' \dontrun{fitRgw <- lme_mass_fit_Rgw(X, Zcols, Y, ni, fitInit$Theta0, RgGrow$Rgs, Surf)}
#' \dontrun{rfx <- lme_mass_rfx(fitRgw$stats, X, Zcols, Y, ni, maskvtx)}

lme_mass_rfx <- function(stats, X, Zcols, Y, ni, maskvtx=NA, prs=1){

    # parallel processing

    do_estimate <- function(j, stats, X, Zcols, maskvtx, nrfx, m){

        #
        print(paste0("Extracting coefficients at vertex ", j), quote=F)

        #
        tryCatch({
            if ((!is.null(stats[[maskvtx[j]]])) && (!is.na(stats[[maskvtx[j]]]$Bhat)) && (!is.na(stats[[maskvtx[j]]]$Dhat)) && (!is.na(stats[[maskvtx[j]]]$phisqhat))) {
                D <- stats[[maskvtx[j]]]$Dhat
                phisq <- stats[[maskvtx[j]]]$phisqhat
                r <- Y[, maskvtx[j]] - X %*% stats[[maskvtx[j]]]$Bhat
                posi <- 1
                scInvD <- solve(D) * phisq
                ijcol <- (maskvtx[j]-1) * nrfx + 1
                fjcol <- ijcol + nrfx - 1
                Rfxj <- matrix(0, nrow=m, ncol=nrfx)
                for (i in c(1:m)) {
                    posf <- posi + ni[i] - 1
                    Zi <- X[posi:posf, Zcols, drop=F]
                    Wi <- (diag(ni[i]) - (lme_mrdivide(Zi, ((t(Zi)%*%Zi) + scInvD), ginv.do=F)%*%t(Zi))) / phisq
                    ri <- r[posi:posf]
                    Rfxj[i, ] <- t(D%*%t(Zi)%*%Wi%*%ri)
                    posi <- posf + 1
                }
                Bhatj <- stats[[maskvtx[j]]]$Bhat[Zcols]
            }
        },
        error=function(e){print(paste("Vertex ", j, ': did not run: ', e, sep=""), quote=F)}
        )

        #
        return(list(Bhatj=Bhatj, Rfxj=Rfxj, ijcol=ijcol, fjcol=fjcol, maskvtxj=maskvtx[j]))
    }

    # Initialization

    print('Initialization ...', quote = F)

    nv0 <-  dim(Y)[2]
    if (any(is.na(maskvtx))) {maskvtx <- seq(1, nv0)}
    m <- length(ni)
    nrfx <- length(Zcols)
    Rfx <- matrix(0, nrow=m, ncol=nv0*nrfx)
    Bhat <- matrix(0, nrow=nrfx, ncol=nv0)

    # Estimation

    print('Extracting coefficients at each vertex ...', quote = F)

    if (prs > 1)
    {
        outEst <- bettermc::mclapply(c(1:length(maskvtx)), function(x) {
            do_estimate(x, stats, X, Zcols, maskvtx, nrfx, m)
        }, mc.cores = prs, mc.progress=FALSE, mc.stdout="output")
        for (j in c(1:length(maskvtx))) {
            Rfx[, outEst[[j]]$ijcol:outEst[[j]]$fjcol] <- outEst[[j]]$Rfxj
            Bhat[, outEst[[j]]$maskvtxj] <- outEst[[j]]$Bhatj
        }
        rm(outEst)
    }
    else
    {
        for (j in c(1:length(maskvtx))) {
            outEst <- do_estimate(j, stats, X, Zcols, maskvtx, nrfx, m)
            Rfx[, outEst$ijcol:outEst$fjcol] <- outEst$Rfxj
            Bhat[, outEst$maskvtxj] <- outEst$Bhatj
            rm(outEst)
        }
    }

    # Output

    print('Storing results ...', quote = F)

    out <- NULL
    out$Rfx <- Rfx
    out$nrfx <- nrfx
    out$Bhat <- Bhat
    return(out)
}
