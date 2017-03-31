#' Generate RNA-seq data
#'
#' From a log-normal distribution
#'
#' @param len number of observations to be generated
#' @param meanlog ln(mean)
#' @param sdlog ln(sd)
#'
#' @return random observations
#' @export
#'
#' @examples
#' gen.rna.seq(1)
#' gen.rna.seq(100)
#' gen.rna.seq(10, log(1.2), log(1 + 0.5))
gen.rna.seq <- function(len,
                        meanlog = log(runif(1, min = 0, max= 3)),
                        sdlog   = log(1 + runif(1, min = 0.08, max= 1.5))) {
  res <- rlnorm(len, meanlog = meanlog, sdlog = sdlog)
  if (meanlog <= log(1e-1)) {
    # higher change if low expressed genes
    res[runif(length(res)) >= 0.9] <- 0
    sdlog <- log(1) + sdlog / 10
  }
  return(res)
}
