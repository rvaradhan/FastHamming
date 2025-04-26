#' Pairwise Hamming distances
#'
#' Computes the pairwise Hamming distances between rows of a binary matrix.
#'
#' @param X A binary (0/1) numeric matrix.
#' @param nthreads Integer; number of OpenMP threads to use. If `NULL` (the default) use all available cores,
#' @return An integer matrix of pairwise Hamming distances.
#' @useDynLib FastHamming, .registration = TRUE
#' @export
#' @examples
#' \donttest{
#' n <- 10000
#' m <- 1000
#' set.seed(2468)
#' X <- matrix(sample(0:1, n * m, replace = TRUE), nrow = n)
#' # Use all available threads
#' system.time(result <- hamming_distance(X))
#' # limit to 2 threads
#' system.time(hamming_distance(X, nthreads = 2))
#' }

hamming_distance <- function(X, nthreads=NULL) {
  stopifnot(is.matrix(X), is.numeric(X))
  if (any(X != 0L & X != 1L)) stop("Matrix X must contain only 0s and 1s.")
  storage.mode(X) <- "integer"
  if (!is.null(nthreads)) {
    stopifnot(length(nthreads) == 1,
              is.numeric(nthreads),
              nthreads >= 1)
    nthreads <- as.integer(nthreads)
  }
  .Call("c_hamming_distance", X, as.integer(nrow(X)), as.integer(ncol(X)),
        nthreads,
        PACKAGE = "FastHamming")
}


