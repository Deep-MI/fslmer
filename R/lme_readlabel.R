#' Read FreeSurfer label files
#'
#' @param fname Filename of a FreeSurfer label file.
#'
#' @return
#' The function returns a nvertices-by-5 matrix, where each column means:
#' (1) vertex number, (2-4) xyz at each vertex, (5) label value. The
#' vertex number is 0-based.
#'
#' @export
#'
#' @examples
#' \dontrun{labels <- lme_readlabel(filename)}

lme_readlabel <- function(fname) {

    fid <- file(fname, 'r')

    readLines(fid, n=1)

    nv <- scan(fid, what = integer(), n = 1)

    d <- simplify2array(scan(fid, what=list(integer(), numeric(), numeric(), numeric(), numeric()), n=5*nv))

    close(fid)

    return(d)

}
