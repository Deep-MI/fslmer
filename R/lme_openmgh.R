#' Read FreeSurfer image (mgh) files
#'
#' @param fname Filename of a FreeSurfer mgh file.
#'
#' @return
#' The function returns a list with entries
#' x, v, ndim1, ndim2, ndim3, nframes, type, and dof.
#'
#' @export
#'
#' @examples
#' \dontrun{mgh <- lme_openmgh(filename)}

lme_openmgh <- function(fname) {

    to.read <- file(fname, "rb")

    v <- readBin(to.read, integer(), endian = "big")
    ndim1 <- readBin(to.read, integer(), endian = "big")
    ndim2 <- readBin(to.read, integer(), endian = "big")
    ndim3 <- readBin(to.read, integer(), endian = "big")
    nframes <- readBin(to.read, integer(), endian = "big")
    type <- readBin(to.read, integer(), endian = "big")
    dof <- readBin(to.read, integer(), endian = "big")

    close(to.read)

    to.read <- file(fname, "rb")
    dump <-
        readBin(to.read,
                double(),
                size = 4,
                n = 71,
                endian = "big")

    x <- array(NA, c(ndim1, ndim2, ndim3, nframes))

    for (k in c(1:nframes)) {
        for (j in c(1:ndim3)) {
            for (i in c(1:ndim2)) {
                x[, i, j, k] <-
                    readBin(
                        to.read,
                        double(),
                        size = 4,
                        n = ndim1,
                        endian = "big"
                    )
            }
        }
    }

    close(to.read)

    return(list(
        x = x,
        v = v,
        ndim1 = ndim1,
        ndim2 = ndim2,
        ndim3 = ndim3,
        nframes = nframes,
        type = type,
        dof = dof
    ))

}
