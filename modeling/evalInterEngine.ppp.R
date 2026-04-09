evalInterEngine.ppp = function (X, P, E, interaction, correction, splitInf = FALSE, 
          ..., Reach = NULL, precomputed = NULL, savecomputed = FALSE) 
{
  fasteval <- interaction$fasteval
  cando <- interaction$can.do.fast
  par <- interaction$par
  feopt <- spatstat.options("fasteval")
  dofast <- !is.null(fasteval) && (is.null(cando) || cando(X, 
                                                           correction, par)) && (feopt %in% c("on", "test")) && 
    (!splitInf || ("splitInf" %in% names(formals(fasteval))))
  V <- NULL
  if (dofast) {
    if (feopt == "test") 
      message("Calling fasteval")
    V <- fasteval(X, P, E, interaction$pot, interaction$par, 
                  correction, splitInf = splitInf, ...)
  }
  if (is.null(V)) {
    evaluate <- interaction$family$eval
    evalargs <- names(formals(evaluate))
    if (splitInf && !("splitInf" %in% evalargs)) 
      stop("Sorry, the", interaction$family$name, "interaction family", 
           "does not support calculation of the positive part", 
           call. = FALSE)
    if (is.null(Reach)) 
      Reach <- reach(interaction)
    if ("precomputed" %in% evalargs) {
      V <- evaluate.ppp(X, P, E, interaction$pot, interaction$par, 
                    correction = correction, splitInf = splitInf, 
                    ..., Reach = Reach, precomputed = precomputed, 
                    savecomputed = savecomputed)
    }
    else {
      V <- evaluate.ppp(X, P, E, interaction$pot, interaction$par, 
                    correction = correction, splitInf = splitInf, 
                    ..., Reach = Reach)
    }
  }
  return(V)
}
