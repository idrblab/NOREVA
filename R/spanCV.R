### Modified from code by ucfagls http://www.r-bloggers.com/whats-wrong-with-loess-for-palaeo-data/
#' @importFrom stats resid
loessGCV <- function(x) {
    if (!(inherits(x, "loess")))
        stop("Error: argument must be a loess object")
    span <- x$pars$span
    n <- x$n
    traceL <- x$trace.hat
    sigma2 <- sum(resid(x)^2)/(n - 1)
    gcv <- n * sigma2/(n - traceL)^2
    result <- list(span = span, gcv = gcv)
    result
}
