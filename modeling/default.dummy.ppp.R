default.dummy.ppp = function (X, nd = NULL, random = FALSE, ntile = NULL, npix = NULL, 
  quasi = FALSE, ..., eps = NULL, verbose = FALSE) 
{
  orig <- list(nd = nd, eps = eps, ntile = ntile, npix = npix)
  orig <- orig[!sapply(orig, is.null)]
  X <- as.ppp(X)
  win <- X$window
  a <- default.n.tiling.ppp(X, nd = nd, ntile = ntile, npix = npix, 
    eps = eps, random = random, quasi = quasi, verbose = verbose)
  nd <- a$nd
  ntile <- a$ntile
  npix <- a$npix
  periodsample <- !quasi && !random && is.mask(win) && all(nd%%win$dim == 
    0)
  #dummy <- if (quasi) 
  #  rQuasi(prod(nd), as.rectangle(win))
  #else if (random) 
  #  stratrand(win, nd[1L], nd[2L], 1)
  #else cellmiddles(win, nd[1L], nd[2L], npix)
  #dummy <- as.ppp(dummy, win, check = FALSE)
  #if (!is.rectangle(win) && !periodsample) 
  #  dummy <- dummy[win]
  corn <- as.ppp(corners(win), win, check = FALSE)
  corn <- corn[win]
  #dummy <- superimpose(dummy, corn, W = win, check = FALSE)
  dummy <- corn
  if (dummy$n == 0) 
    stop("None of the dummy points lies inside the window")
  attr(dummy, "weight.parameters") <- append(list(...), list(ntile = ntile, 
    verbose = verbose, npix = npix))
  attr(dummy, "dummy.parameters") <- list(nd = nd, random = random, 
    quasi = quasi, verbose = verbose, orig = orig)
  return(dummy)
}

