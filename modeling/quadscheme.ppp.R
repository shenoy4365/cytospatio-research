quadscheme.ppp = function (data, dummy, method = "d", ...) 
{
  data <- as.ppp(data)
  print(data)
  if (missing(dummy)) {
    dummy <- default.dummy.ppp(data, method = method,  ...)
    dp <- attr(dummy, "dummy.parameters")
    wp <- attr(dummy, "weight.parameters")
  }
  else {
    if (!is.ppp(dummy)) {
      dummy <- as.ppp(dummy, data$window, check = FALSE)
      dummy <- dummy[data$window]
      wp <- dp <- list()
    }
    else {
      dp <- attr(dummy, "dummy.parameters")
      wp <- attr(dummy, "weight.parameters")
    }
  }
  wp <- resolve.defaults(list(method = method), list(...), 
    wp)
  mX <- is.marked(data)
  mD <- is.marked(dummy)
  if (!mX && !mD) 
    Q <- do.call(quadscheme.spatial, append(list(data, dummy, 
      check = FALSE), wp))
  else if (mX && !mD) 
    Q <- do.call(quadscheme.replicated, append(list(data, 
      dummy, check = FALSE), wp))
  else if (!mX && mD) 
    stop("dummy points are marked but data are unmarked")
  else stop("marked data and marked dummy points -- sorry, this case is not implemented")
  Q$param$dummy <- dp
  return(Q)
}

