default.n.tiling.ppp = function (X, nd = NULL, ntile = NULL, npix = NULL, eps = NULL, 
  random = FALSE, quasi = FALSE, verbose = TRUE) 
{
  verifyclass(X, "ppp")
  win <- X$window
  pixels <- (win$type != "rectangle")
  if (nd.given <- !is.null(nd)) 
    nd <- ensure2print(nd, verbose)
  if (ntile.given <- !is.null(ntile)) 
    ntile <- ensure2print(ntile, verbose)
  if (npix.given <- !is.null(npix)) 
    npix <- ensure2print(npix, verbose)
  if (pixels) 
    sonpixel <- rev(ensure2print(spatstat.options("npixel"), 
      verbose, ""))
  ndummy.min <- ensure2print(spatstat.options("ndummy.min"), 
    verbose, "")
  #print(ndummy.min)
  #print(5 * ceiling(2 * sqrt(X$n)/10))
  #ndminX <- pmax(ndummy.min, 5 * ceiling(2 * sqrt(X$n)/10))
  #ndminX <- ensure2vector(ndminX)
  ndminX <- 0 
  if (!is.null(eps)) {
    eps <- ensure2print(eps, verbose)
    Xbox <- as.rectangle(as.owin(X))
    sides <- with(Xbox, c(diff(xrange), diff(yrange)))
    ndminX <- pmax(ndminX, ceiling(sides/eps))
  }
  if (npix.given) 
    Nmin <- Nmax <- npix
  else switch(win$type, rectangle = {
    Nmin <- ensure2vector(X$n)
    Nmax <- Inf
  }, polygonal = {
    Nmin <- sonpixel
    Nmax <- 4 * sonpixel
  }, mask = {
    nmask <- rev(win$dim)
    Nmin <- nmask
    Nmax <- pmax(2 * nmask, 4 * sonpixel)
  })
  if (nd.given && !ntile.given) {
    if (any(nd > Nmax)) 
      warning("number of dummy points nd exceeds maximum pixel dimensions")
    ntile <- min2div(nd, ndminX, nd)
  }
  else if (!nd.given && ntile.given) {
    nd <- min2mul(ntile, ndminX, Nmin)
    if (any(nd >= Nmin)) 
      nd <- ntile
  }
  else if (!nd.given && !ntile.given) {
    if (!pixels) {
      nd <- ntile <- ensure2vector(ndminX)
      if (verbose) 
        cat(paste("nd and ntile default to", nd[1L], 
          "*", nd[2L], "\n"))
    }
    else {
      nd <- ntile <- min2div(Nmin, ndminX, Nmax)
      if (any(nd >= Nmin)) {
        if (verbose) 
          cat("No suitable divisor of pixel dimensions\n")
        nd <- ntile <- ndminX
      }
    }
  }
  else {
    if (any(ntile > nd)) 
      warning("the number of tiles (ntile) exceeds the number of dummy points (nd)")
  }
  if (!ntile.given && quasi) {
    if (verbose) 
      cat("Adjusting ntile because quasi=TRUE\n")
    ntile <- maxdiv(ntile, if (pixels) 
      2L
    else 1L)
  }
  if (!npix.given && pixels) 
    npix <- min2mul(nd, Nmin, Nmax)
  if (verbose) {
    if (!quasi) 
      cat(paste("dummy points:", paste0(if (random) 
        "stratified random in"
      else NULL, "grid"), nd[1L], "x", nd[2L], "\n"))
    else cat(paste("dummy points:", nd[1L], "x", nd[2L], 
      "=", prod(nd), "quasirandom points\n"))
    cat(paste("weighting tiles", ntile[1L], "x", ntile[2L], 
      "\n"))
    if (pixels) 
      cat(paste("pixel grid", npix[1L], "x", npix[2L], 
        "\n"))
  }
  if (pixels) 
    return(list(nd = nd, ntile = ntile, npix = npix))
  else return(list(nd = nd, ntile = ntile, npix = npix))
}
ensure2print <- function(x, verbose=TRUE, blah="user specified") {
  xname <- short.deparse(substitute(x))
  x <- ensure2vector(x)
  if(verbose)
    cat(paste(blah, xname, "=", x[1L], "*", x[2L], "\n"))
  x
}

