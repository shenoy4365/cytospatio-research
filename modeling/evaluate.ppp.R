evaluate.ppp = function (X, U, EqualPairs, pairpot, potpars, correction, splitInf = FALSE, 
                           ..., Reach = NULL, precomputed = NULL, savecomputed = FALSE, 
                           pot.only = FALSE) 
{
  pt <- PairPotentialType(pairpot)
  
  if (length(correction) > 1) 
    stop("Only one edge correction allowed at a time!")
  if (!any(correction == c("periodic", "border", "translate", 
                           "translation", "isotropic", "Ripley", "none"))) 
    stop(paste("Unrecognised edge correction", sQuote(correction)))
  no.correction <- use.closepairs <- (correction %in% c("none", 
                                                        "border", "translate", "translation")) && !is.null(Reach) && 
    is.finite(Reach) && is.null(precomputed) && !savecomputed
  if (!is.null(precomputed)) {
    X <- precomputed$X
    U <- precomputed$U
    EqualPairs <- precomputed$E
    M <- precomputed$M
  }
  else {
    U <- as.ppp(U, X$window)
    if (!use.closepairs) 
      M <- crossdist(X, U, periodic = (correction == "periodic"))
  }
  nX <- npoints(X)
  nU <- npoints(U)
  dimM <- c(nX, nU)
  if (use.closepairs) {
    V <- evalPairPotential.ppp(X, U, EqualPairs, pairpot,
                             potpars, Reach)
  }
  else {
    V <- do.call.matched(pairpot, list(d = M, tx = marks(X),
                                         tu = marks(U), par = potpars))
  }
  # if (use.closepairs) {
  #   POT <- evalPairPotential.pp3(X, U, EqualPairs, pairpot, 
  #                            potpars, Reach)
  # }
  # else {
  #   POT <- do.call.matched(pairpot, list(d = M, tx = marks(X), 
  #                                        tu = marks(U), par = potpars))
  # }
  # IsOffset <- attr(POT, "IsOffset")
  # if (!is.matrix(POT) && !is.array(POT)) {
  #   if (length(POT) == 0 && X$n == 0) 
  #     POT <- array(POT, dim = c(dimM, 1))
  #   else stop("Pair potential did not return a matrix or array")
  # }
  # if (length(dim(POT)) == 1 || any(dim(POT)[1:2] != dimM)) {
  #   whinge <- paste0("The pair potential function ", short.deparse(substitute(pairpot)), 
  #                    " must produce a matrix or array with its first two dimensions\n", 
  #                    "the same as the dimensions of its input.\n")
  #   stop(whinge)
  # }
  # if (length(dim(POT)) == 2) 
  #   POT <- array(POT, dim = c(dim(POT), 1), dimnames = NULL)
  # if (splitInf) {
  #   IsNegInf <- (POT == -Inf)
  #   POT[IsNegInf] <- 0
  # }
  # if (correction == "translate" || correction == "translation") {
  #   edgewt <- edge.Trans(X, U)
  #   if (!is.matrix(edgewt)) 
  #     stop("internal error: edge.Trans() did not yield a matrix")
  #   if (nrow(edgewt) != X$n || ncol(edgewt) != length(U$x)) 
  #     stop("internal error: edge weights matrix returned by edge.Trans() has wrong dimensions")
  #   POT <- c(edgewt) * POT
  # }
  # else if (correction == "isotropic" || correction == "Ripley") {
  #   edgewt <- t(edge.Ripley.pp3(U, t(M), X$window))
  #   if (!is.matrix(edgewt)) 
  #     stop("internal error: edge.Ripley() did not return a matrix")
  #   if (nrow(edgewt) != X$n || ncol(edgewt) != length(U$x)) 
  #     stop("internal error: edge weights matrix returned by edge.Ripley() has wrong dimensions")
  #   POT <- c(edgewt) * POT
  # }
  # if (length(EqualPairs) > 0) {
  #   nplanes <- dim(POT)[3]
  #   for (k in 1:nplanes) {
  #     POT[cbind(EqualPairs, k)] <- 0
  #     if (splitInf) 
  #       IsNegInf[cbind(EqualPairs, k)] <- FALSE
  #   }
  # }
  # if (splitInf) 
  #   attr(POT, "IsNegInf") <- IsNegInf
  # if (pot.only) 
  #   return(POT)
  # V <- apply(POT, c(2, 3), sum)
  # if (splitInf) 
  #   attr(V, "-Inf") <- apply(IsNegInf, 2, any)
  # # attr(V, "POT") <- POT
  # rm(POT)
  # gc()
  # attr(V, "IsOffset") <- IsOffset
  # if (savecomputed) 
  #   attr(V, "computed") <- list(E = EqualPairs, M = M)
  return(V)
}

