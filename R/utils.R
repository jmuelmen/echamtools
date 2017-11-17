#' A function to label vertically contiguous features (internal use only in cfodd.R and cosp.R)
label.vertical.features <- function(vfm) {
    x <- vfm
    if (length(x) == 0)
        return(x)
    diff.x <- diff(c(-1, x)) ## guarantee that the first group of 1's is preceded by a transition
    labels <- cumsum(diff.x != 0 & x != 0) * (x != 0) ## count up the edges
    labels
}

