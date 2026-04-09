evalInteraction.ppp = function (X, P, E = equalpairs(P, X), interaction, correction, 
          splitInf = FALSE, ..., precomputed = NULL, savecomputed = FALSE) 
{
  verifyclass(interaction, "interact")
  if (is.poisson(interaction)) {
    out <- matrix(numeric(0), nrow = npoints(P), ncol = 0)
    attr(out, "IsOffset") <- logical(0)
    if (splitInf) 
      attr(out, "-Inf") <- logical(nrow(out))
    return(out)
  }
  dofast <- (spatstat.options("fasteval") %in% c("on", "test")) && 
    !is.null(cando <- interaction$can.do.fast) && cando(X, correction, interaction$par) && !is.null(interaction$fasteval)
  if (dofast) {
    dosplit <- FALSE
  }
  else {
    needsplit <- oversize.quad(nU = npoints(P), nX = npoints(X))
    dosplit <- needsplit && !savecomputed
    if (needsplit && savecomputed) 
      warning(paste("Oversize quadscheme cannot be split into blocks", 
                    "because savecomputed=TRUE;", "memory allocation error may occur"))
  }
  if (!dosplit) {
    V <- evalInterEngine(X = X, P = P, E = E, interaction = interaction, 
                         correction = correction, splitInf = splitInf, ..., 
                         precomputed = precomputed, savecomputed = savecomputed)
  }
  else {
    nX <- npoints(X)
    nP <- npoints(P)
    Pdata <- E[, 2]
    Pall <- seq_len(nP)
    Pdummy <- if (length(Pdata) > 0) 
      Pall[-Pdata]
    else Pall
    nD <- length(Pdummy)
    bls <- quadBlockSizes.ppp(nX, nD, announce = TRUE)
    nblocks <- bls$nblocks
    nperblock <- bls$nperblock
    seqX <- seq_len(nX)
    # EX <- cbind(seqX, seqX)
    EX = NULL
    V_list = list()
    for (iblock in 0:nblocks) {
      print(iblock)
      if (iblock == 0){
        Pi <- X
      } else{
        first <- min(nD, (iblock - 1) * nperblock + 1)
        last <- min(nD, iblock * nperblock)
        Di <- P[Pdummy[first:last]]
        Pi = Di
      }
      Vi <- evalInterEngine.ppp(X = X, P = Pi, E = EX, interaction = interaction, 
                            correction = correction, splitInf = splitInf, 
                            ..., savecomputed = FALSE)
      Mi <- attr(Vi, "-Inf")
      V_list[[iblock+1]] = as.data.frame(Vi)
      if (iblock == 0) {
        # V <- Vi
        M <- Mi
      }
      else {
        # V <- rbind(V, Vi[-seqX, , drop = FALSE])
        # print('test')
        # V <- rbind(V, Vi)
        # V <- as.matrix(rbindlist(list(as.data.frame(V), as.data.frame(Vi[-seqX, , drop = FALSE]))))
        if (splitInf && !is.null(M)) 
          M <- c(M, Mi[-seqX])
      }
      rm(Vi)
      rm(Mi)
      gc()
    }
    V = rbindlist(V_list)
    V = as.matrix(V)
    if (length(Pdata) == 0) {
      V <- V[-seqX, , drop = FALSE]
      if (splitInf && !is.null(M)) 
        M <- M[-seqX]
    }
    else {
      ii <- integer(nP)
      ii[Pdata] <- seqX
      ii[Pdummy] <- (nX + 1):nrow(V)
      V <- V[ii, , drop = FALSE]
      if (splitInf && !is.null(M)) 
        M <- M[ii]
    }
    attr(V, "-Inf") <- M
  }
  return(V)
}
