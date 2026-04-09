mpl.prepare.ppp <- function(Q, X, P, trend, interaction, covariates, 
                        want.trend, want.inter, correction, rbord,
                        Pname="quadrature points", callstring="",
                        ...,
                        subsetexpr=NULL,
                        covfunargs=list(),
                        allcovar=FALSE,
                        precomputed=NULL, savecomputed=FALSE,
                        vnamebase=c("Interaction", "Interact."),
                        vnameprefix=NULL,
                        warn.illegal=TRUE,
                        warn.unidentifiable=TRUE,
                        weightfactor=NULL,
                        skip.border=FALSE,
                        clip.interaction=TRUE,
                        splitInf=FALSE) {
  ## Q: quadrature scheme
  ## X = data.quad(Q)
  ## P = union.quad(Q)
  
    if (missing(want.trend)) 
      want.trend <- !is.null(trend) && !identical.formulae(trend, 
                                                           ~1)
    if (missing(want.inter)) 
      want.inter <- !is.null(interaction) && !is.null(interaction$family)
    want.subset <- !is.null(subsetexpr)
    computed <- list()
    problems <- list()
    names.precomputed <- names(precomputed)
    likelihood.is.zero <- FALSE
    is.identifiable <- TRUE
    if (!missing(vnamebase)) {
      if (length(vnamebase) == 1) 
        vnamebase <- rep.int(vnamebase, 2)
      if (!is.character(vnamebase) || length(vnamebase) != 
          2) 
        stop("Internal error: illegal format of vnamebase")
    }
    if (!is.null(vnameprefix)) {
      if (!is.character(vnameprefix) || length(vnameprefix) != 
          1) 
        stop("Internal error: illegal format of vnameprefix")
    }
    updatecovariates <- FALSE
    covariates.df <- NULL
    if (allcovar || want.trend || want.subset) {
      if ("covariates.df" %in% names.precomputed) {
        covariates.df <- precomputed$covariates.df
      }
      else {
        if (!is.data.frame(covariates)) {
          covnames <- variablesinformula(trend)
          if (want.subset) 
            covnames <- union(covnames, all.vars(subsetexpr))
          if (allcovar) 
            covnames <- union(covnames, names(covariates))
          covnames <- setdiff(covnames, c("x", "y", "marks"))
          tenv <- environment(trend)
          covariates <- getdataobjects(covnames, tenv, 
                                       covariates, fatal = TRUE)
          updatecovariates <- any(attr(covariates, "external"))
        }
        covariates.df <- mpl.get.covariates(covariates, 
                                            P, Pname, covfunargs)
      }
      if (savecomputed) 
        computed$covariates.df <- covariates.df
    }
    if ("dotmplbase" %in% names.precomputed) 
      .mpl <- precomputed$dotmplbase
    else {
      nQ <- n.quad(Q)
      wQ <- w.quad(Q)
      mQ <- marks.quad(Q)
      zQ <- is.data(Q)
      yQ <- numeric(nQ)
      yQ[zQ] <- 1/wQ[zQ]
      zeroes <- attr(wQ, "zeroes")
      sQ <- if (is.null(zeroes)) 
        rep.int(TRUE, nQ)
      else !zeroes
      if (!is.null(weightfactor)) 
        wQ <- wQ * weightfactor
      .mpl <- list(W = wQ, Z = zQ, Y = yQ, MARKS = mQ, SUBSET = sQ)
    }
    if (savecomputed) 
      computed$dotmplbase <- .mpl
    glmdata <- data.frame(.mpl.W = .mpl$W, .mpl.Y = .mpl$Y)
    izdat <- .mpl$Z[.mpl$SUBSET]
    ndata <- sum(izdat)
    if (correction == "border") {
      bdP <- if ("bdP" %in% names.precomputed) 
        precomputed$bdP
      else bdist.points(P)
      if (savecomputed) 
        computed$bdP <- bdP
      .mpl$DOMAIN <- (bdP >= rbord)
    }
    skip.border <- skip.border && (correction == "border")
    internal.names <- c(".mpl.W", ".mpl.Y", ".mpl.Z", ".mpl.SUBSET", 
                        "SUBSET", ".mpl")
    reserved.names <- c("x", "y", "marks", internal.names)
    if (allcovar || want.trend || want.subset) {
      trendvariables <- variablesinformula(trend)
      cc <- check.clashes(internal.names, trendvariables, 
                          "the model formula")
      if (cc != "") 
        stop(cc)
      if (want.subset) {
        subsetvariables <- all.vars(subsetexpr)
        cc <- check.clashes(internal.names, trendvariables, 
                            "the subset expression")
        if (cc != "") 
          stop(cc)
        trendvariables <- union(trendvariables, subsetvariables)
      }
      if (allcovar || "x" %in% trendvariables) 
        glmdata <- data.frame(glmdata, x = P$x)
      if (allcovar || "y" %in% trendvariables) 
        glmdata <- data.frame(glmdata, y = P$y)
      if (("marks" %in% trendvariables) || !is.null(.mpl$MARKS)) {
        if (is.null(.mpl$MARKS)) 
          stop("Model formula depends on marks, but data do not have marks", 
               call. = FALSE)
        glmdata <- data.frame(glmdata, marks = .mpl$MARKS)
      }
      if (!is.null(covariates.df)) {
        cc <- check.clashes(reserved.names, names(covariates), 
                            sQuote("covariates"))
        if (cc != "") 
          stop(cc)
        if (!allcovar) 
          needed <- names(covariates.df) %in% trendvariables
        else needed <- rep.int(TRUE, ncol(covariates.df))
        if (any(needed)) {
          covariates.needed <- covariates.df[, needed, 
                                             drop = FALSE]
          glmdata <- data.frame(glmdata, covariates.needed)
          nbg <- is.na(covariates.needed)
          if (any(nbg)) {
            offending <- matcolany(nbg)
            covnames.na <- names(covariates.needed)[offending]
            quadpoints.na <- matrowany(nbg)
            n.na <- sum(quadpoints.na)
            n.tot <- length(quadpoints.na)
            errate <- n.na/n.tot
            pcerror <- round(signif(100 * errate, 2), 
                             2)
            complaint <- paste("Values of the", ngettext(length(covnames.na), 
                                                         "covariate", "covariates"), paste(sQuote(covnames.na), 
                                                                                           collapse = ", "), "were NA or undefined at", 
                               paste(pcerror, "%", " (", n.na, " out of ", 
                                     n.tot, ")", sep = ""), "of the", Pname)
            warning(paste(complaint, ". Occurred while executing: ", 
                          callstring, sep = ""), call. = FALSE)
            .mpl$SUBSET <- .mpl$SUBSET & !quadpoints.na
            details <- list(covnames.na = covnames.na, 
                            quadpoints.na = quadpoints.na, print = complaint)
            problems <- append(problems, list(na.covariates = details))
          }
        }
      }
    }
  
  ###################### I n t e r a c t i o n ####################
  
  Vnames <- NULL
  IsOffset <- NULL
  forbid <- NULL
  
  if(want.inter) {
    ## Form the matrix of "regression variables" V.
    ## The rows of V correspond to the rows of P (quadrature points)
    ## while the column(s) of V are the regression variables (log-potentials)
    
    E <- precomputed$E %orifnull% equalpairs.quad(Q)
    
    if(!skip.border) {
      ## usual case
      V <- evalInteraction.ppp(X, P, E, interaction, correction,
                           ...,
                           splitInf=splitInf,
                           precomputed=precomputed,
                           savecomputed=savecomputed)
    } else {
      ## evaluate only in eroded domain
      if(all(c("Esub", "Usub", "Retain") %in% names.precomputed)) {
        ## use precomputed data
        Psub <- precomputed$Usub
        Esub <- precomputed$Esub
        Retain <- precomputed$Retain
      } else {
        ## extract subset of quadrature points
        Retain <- .mpl$DOMAIN | is.data(Q)
        Psub <- P[Retain]
        ## map serial numbers in P to serial numbers in Psub
        Pmap <- cumsum(Retain)
        ## extract subset of equal-pairs matrix
        keepE <- Retain[ E[,2] ]
        Esub <- E[ keepE, , drop=FALSE]
        ## adjust indices in equal pairs matrix
        Esub[,2] <- Pmap[Esub[,2]]
      }
      ## call evaluator on reduced data
      if(all(c("X", "Q", "U") %in% names.precomputed)) {
        subcomputed <- resolve.defaults(list(E=Esub, U=Psub, Q=Q[Retain]),
                                        precomputed)
      } else subcomputed <- NULL
      if(clip.interaction) {
        ## normal
        V <- evalInteraction(X, Psub, Esub, interaction, correction,
                             ...,
                             splitInf=splitInf,
                             precomputed=subcomputed,
                             savecomputed=savecomputed)
      } else {
        ## ignore window when calculating interaction
        ## by setting 'W=NULL' (currently detected only by AreaInter)
        V <- evalInteraction(X, Psub, Esub, interaction, correction,
                             ...,
                             W=NULL,
                             splitInf=splitInf,
                             precomputed=subcomputed,
                             savecomputed=savecomputed)
      }
      if(savecomputed) {
        computed$Usub <- Psub
        computed$Esub <- Esub
        computed$Retain <- Retain
      }
    }
    
    if(!is.matrix(V))
      stop("interaction evaluator did not return a matrix")
    
    ## extract information about offsets
    IsOffset <- attr(V, "IsOffset")
    if(is.null(IsOffset)) IsOffset <- FALSE
    
    if(splitInf) {
      ## extract information about hard core terms
      forbid <- attr(V, "-Inf") %orifnull% logical(nrow(V))
    } 
    
    if(skip.border) {
      ## fill in the values in the border region with zeroes.
      Vnew <- matrix(0, nrow=npoints(P), ncol=ncol(V))
      colnames(Vnew) <- colnames(V)
      Vnew[Retain, ] <- V
      ## retain attributes
      attr(Vnew, "IsOffset") <- IsOffset
      attr(Vnew, "computed") <- attr(V, "computed")
      attr(Vnew, "POT") <- attr(V, "POT")
      V <- Vnew
      if(splitInf) {
        fnew <- logical(nrow(Vnew))
        fnew[Retain] <- forbid
        forbid <- fnew
      }
    }
    
    ## extract intermediate computation results 
    if(savecomputed)
      computed <- resolve.defaults(attr(V, "computed"), computed)
    
    ## Augment data frame by appending the regression variables
    ## for interactions.
    ##
    ## First determine the names of the variables
    ##
    Vnames <- dimnames(V)[[2]]
    if(is.null(Vnames)) {
      ## No names were provided for the columns of V.
      ## Give them default names.
      ## In ppm the names will be "Interaction"
      ##   or "Interact.1", "Interact.2", ...
      ## In mppm an alternative tag will be specified by vnamebase.
      nc <- ncol(V)
      Vnames <- if(nc == 1) vnamebase[1] else paste0(vnamebase[2], 1:nc)
      dimnames(V) <- list(dimnames(V)[[1]], Vnames)
    } else if(!is.null(vnameprefix)) {
      ## Variable names were provided by the evaluator (e.g. MultiStrauss).
      ## Prefix the variable names by a string
      ## (typically required by mppm)
      Vnames <- paste(vnameprefix, Vnames, sep="")
      dimnames(V) <- list(dimnames(V)[[1]], Vnames)
    }
    
    ## Check the names are valid as column names in a dataframe
    okVnames <- make.names(Vnames, unique=TRUE)
    if(any(Vnames != okVnames)) {
      warning(paste("Names of interaction terms",
                    "contained illegal characters;",
                    "names have been repaired."))
      Vnames <- okVnames
    }
    
    ##   Check for name clashes between the interaction variables
    ##   and the formula
    cc <- check.clashes(Vnames, termsinformula(trend), "model formula")
    if(cc != "") stop(cc)
    ##   and with the variables in 'covariates'
    if(!is.null(covariates)) {
      cc <- check.clashes(Vnames, names(covariates), sQuote("covariates"))
      if(cc != "") stop(cc)
    }
    
    ## OK. append variables.
    glmdata <- data.frame(glmdata, V)   
    
    ## check IsOffset matches Vnames
    if(length(IsOffset) != length(Vnames)) {
      if(length(IsOffset) == 1)
        IsOffset <- rep.int(IsOffset, length(Vnames))
      else
        stop("Internal error: IsOffset has wrong length", call.=FALSE)
    }
    
    ## Keep only those quadrature points for which the
    ## conditional intensity is nonzero. 
    
    ##KEEP  <- apply(V != -Inf, 1, all)
    .mpl$KEEP  <- matrowall(V != -Inf)
    
    .mpl$SUBSET <- .mpl$SUBSET & .mpl$KEEP
    
    ## Check that there are at least some data and dummy points remaining
    datremain <- .mpl$Z[.mpl$SUBSET]
    somedat <- any(datremain)
    somedum <- !all(datremain)
    if(warn.unidentifiable && !(somedat && somedum)) {
      ## Model would be unidentifiable if it were fitted.
      ## Register problem
      is.identifiable <- FALSE
      if(ndata == 0) {
        complaint <- "model is unidentifiable: data pattern is empty"
      } else {
        offending <- !c(somedat, somedum)
        offending <- c("all data points", "all dummy points")[offending]
        offending <- paste(offending, collapse=" and ")
        complaint <- paste("model is unidentifiable:",
                           offending, "have zero conditional intensity")
      }
      details <- list(data=!somedat,
                      dummy=!somedum,
                      print=complaint)
      problems <- append(problems, list(unidentifiable=details))
    }
    
    ## check whether the model has zero likelihood:
    ## check whether ANY data points have zero conditional intensity
    if(any(.mpl$Z & !.mpl$KEEP)) {
      howmany <- sum(.mpl$Z & !.mpl$KEEP)
      complaint <- paste(howmany,
                         "data point(s) are illegal",
                         "(zero conditional intensity under the model)")
      details <- list(illegal=howmany,
                      print=complaint)
      problems <- append(problems, list(zerolikelihood=details))
      if(warn.illegal && is.identifiable)
        warning(paste(complaint,
                      ". Occurred while executing: ",
                      callstring, sep=""),
                call. = FALSE)
      likelihood.is.zero <- TRUE
    }
  }
  
  ##################     S u b s e t   ###################
  
  if(correction == "border") 
    .mpl$SUBSET <- .mpl$SUBSET & .mpl$DOMAIN
  
  if(!is.null(subsetexpr)) {
    ## user-defined subset expression
    USER.SUBSET <- eval(subsetexpr, glmdata, environment(trend))
    if(is.owin(USER.SUBSET)) {
      USER.SUBSET <- inside.owin(P$x, P$y, USER.SUBSET)
    } else if(is.im(USER.SUBSET)) {
      USER.SUBSET <- as.logical(USER.SUBSET[P, drop=FALSE])
      if(anyNA(USER.SUBSET)) 
        USER.SUBSET[is.na(USER.SUBSET)] <- FALSE
    }
    if(!(is.logical(USER.SUBSET) || is.numeric(USER.SUBSET)))
      stop("Argument 'subset' should yield logical values", call.=FALSE)
    if(anyNA(USER.SUBSET)) {
      USER.SUBSET[is.na(USER.SUBSET)] <- FALSE
      warning("NA values in argument 'subset' were changed to FALSE",
              call.=FALSE)
    }
    .mpl$SUBSET <- .mpl$SUBSET & USER.SUBSET
  }
  
  glmdata <- cbind(glmdata,
                   data.frame(.mpl.SUBSET=.mpl$SUBSET,
                              stringsAsFactors=FALSE))
  
  #################  F o r m u l a   ##################################
  
  if(!want.trend) trend <- ~1 
  trendpart <- paste(as.character(trend), collapse=" ")
  if(!want.inter)
    rhs <- trendpart
  else {
    VN <- Vnames
    ## enclose offset potentials in 'offset(.)'
    if(any(IsOffset))
      VN[IsOffset] <- paste("offset(", VN[IsOffset], ")", sep="")
    rhs <- paste(c(trendpart, VN), collapse= "+")
  }
  fmla <- paste(".mpl.Y ", rhs)
  fmla <- as.formula(fmla)
  
  ##  character string of trend formula (without Vnames)
  trendfmla <- paste(".mpl.Y ", trendpart)
  
  ####
  result <- list(fmla=fmla, trendfmla=trendfmla,
                 covariates=if(updatecovariates) covariates else NULL,
                 glmdata=glmdata, Vnames=Vnames, IsOffset=IsOffset,
                 subsetexpr=subsetexpr,
                 problems=problems,
                 likelihood.is.zero=likelihood.is.zero,
                 is.identifiable=is.identifiable,
                 computed=computed,
                 vnamebase=vnamebase, vnameprefix=vnameprefix,
                 forbid=forbid)
  return(result)
}

check.clashes <- function(forbidden, offered, where) {
  name.match <- outer(forbidden, offered, "==")
  if(any(name.match)) {
    is.matched <- apply(name.match, 2, any)
    matched.names <- (offered)[is.matched]
    if(sum(is.matched) == 1) {
      return(paste("The variable",sQuote(matched.names),
                   "in", where, "is a reserved name"))
    } else {
      return(paste("The variables",
                   paste(sQuote(matched.names), collapse=", "),
                   "in", where, "are reserved names"))
    }
  }
  return("")
}