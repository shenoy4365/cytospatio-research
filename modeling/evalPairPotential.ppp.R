evalPairPotential.ppp = function (X, P, E, pairpot, potpars, R) 
{
  nX <- npoints(X)
  nP <- npoints(P)
  pt <- PairPotentialType(pairpot)
  fakePOT <- do.call.matched(pairpot, list(d = matrix(, 0, 0), tx = marks(X)[integer(0)], tu = marks(P)[integer(0)], par = potpars))
  IsOffset <- attr(fakePOT, "IsOffset")
  fakePOT <- ensure3Darray(fakePOT)
  Vnames <- dimnames(fakePOT)[[3]]
  p <- dim(fakePOT)[3]
  cl <- crosspairs(X, P, R, what = "ijd")
  I <- cl$i
  J <- cl$j
  D <- cl$d

  # if (identical(data.frame(X), data.frame(P))){
  #   same_value_indices <- which(I == J)
  #   I <- I[-same_value_indices]
  #   J <- J[-same_value_indices]
  #   D <- D[-same_value_indices]
  # } 
  same_value_indices <- which(D == 0)
  I <- I[-same_value_indices]
  J <- J[-same_value_indices]
  D <- D[-same_value_indices]
  

  D <- matrix(D, ncol = 1)
  
  if (nX == 0 || nP == 0 || length(I) == 0) {
    di <- c(nP, p)
    dn <- list(NULL, Vnames)
    result <- as.matrix(array(0, dim = di, dimnames = dn))
    # attr(result, "IsOffset") <- IsOffset
    return(result)
  }
  if (!pt$marked) {
    POT <- do.call.matched(pairpot, list(d = D, par = potpars))
    IsOffset <- attr(POT, "IsOffset")
  }
  else {
    marX <- marks(X)
    marP <- marks(P)
    if (!identical(levels(marX), levels(marP))) 
      stop("Internal error: marks of X and P have different levels")
    types <- levels(marX)
    mI <- marX[I]
    mJ <- marP[J]
    POT <- NULL
    for (k in types) {
      relevant <- which(mJ == k)
      if (length(relevant) > 0) {
        fk <- factor(k, levels = types)
        POTk <- do.call.matched(pairpot, list(d = D[relevant, 
                                                    , drop = FALSE], tx = mI[relevant], tu = fk, 
                                              par = potpars))

        if (is.null(POT)) {
          POT <- array(, dim = c(length(I), 1, dim(POTk)[3]))
          IsOffset <- attr(POTk, "IsOffset")
          Vnames <- dimnames(POTk)[[3]]
        }
        POT[relevant, , ] <- POTk
      }
    }
  }
  POT <- ensure3Darray(POT)
  p <- dim(POT)[3]
  
  
  # result_list = list()
  # for (p_index in 1:P$n){
  #   result_list[[p_index]] = as.list(get_interactions(p_index, POT, J, p))
  # }
  # POT[!is.finite(POT)] = 0
  
  # lookup <- lapply(1:P$n, function(x) which(J == x))
  # lookup <- split(seq_along(J), J)
  # lookup <- lookup[1:P$n]
  lookup <- vector("list", P$n)
  indices <- split(seq_along(J), J)
  
  for(i in seq_along(indices)) {
    lookup[[as.numeric(names(indices[i]))]] <- indices[[i]]
  }
  # print(length(unique(I)))
  # print(length(unique(J)))
  # print(X)
  # print(P)
  # print(R)
  # print(1:P$n)
  
  result <- lapply(1:P$n, function(p_index) as.list(get_interactions(p_index, POT, p, lookup)))
  
  result = rbindlist(result)

  result = as.matrix(result)
  colnames(result) = Vnames
  rownames(result) = NULL
  # result <- array(0, dim = c(npoints(X), npoints(P), p), dimnames = list(NULL, 
  #                                                                        NULL, Vnames))
  # II <- rep(I, p)
  # JJ <- rep(J, p)
  # KK <- rep(1:p, each = length(I))
  # IJK <- cbind(II, JJ, KK)
  # result[IJK] <- POT
  # if (length(E) > 0) {
  #   E.rep <- apply(E, 2, rep, times = p)
  #   p.rep <- rep(1:p, each = nrow(E))
  #   result[cbind(E.rep, p.rep)] <- 0
  # }
  # attr(result, "IsOffset") <- IsOffset
  return(result)
}

get_interactions <- function(index, POT, p, lookup) {
  index_list <- lookup[[index]]
  POT_subset <- POT[index_list, 1, 1:p]
  
  if (!is.vector(POT_subset)) {
    count_list <- colSums(POT_subset)
  } else{
    count_list <- POT_subset
  }
  return(count_list)
}
