quadBlockSizes.ppp = function (nX, nD, p = 1, nMAX = spatstat.options("maxmatrix")/p, 
          announce = TRUE) 
{
  if (is.quad(nX) && missing(nD)) {
    nD <- npoints(nX$dummy)
    nX <- npoints(nX$data)
  }
  # nperblock <- max(1, abs(floor(nMAX/nX - nX)))
  nperblock <- 5000
  nblocks <- ceiling(nD/nperblock)
  nperblock <- min(nperblock, ceiling(nD/nblocks))
  if (announce && nblocks > 1) {
    msg <- paste("Large quadrature scheme", "split into blocks to avoid memory size limits;", 
                 nD, "dummy points split into", nblocks, "blocks,")
    nfull <- nblocks - 1
    nlastblock <- nD - nperblock * nfull
    if (nlastblock == nperblock) {
      msg <- paste(msg, "each containing", nperblock, 
                   "dummy points")
    }
    else {
      msg <- paste(msg, "the first", ngettext(nfull, "block", 
                                              paste(nfull, "blocks")), "containing", nperblock, 
                   ngettext(nperblock, "dummy point", "dummy points"), 
                   "and the last block containing", nlastblock, 
                   ngettext(nlastblock, "dummy point", "dummy points"))
    }
    message(msg)
  }
  else nlastblock <- nperblock
  return(list(nblocks = nblocks, nperblock = nperblock, nlastblock = nlastblock))
}
