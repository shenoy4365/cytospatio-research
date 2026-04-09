mppm.ppp = function (formula, data, interaction = Poisson(), ..., iformula = NULL, 
          random = NULL, weights = NULL, use.gam = FALSE, reltol.pql = 0.001, 
          gcontrol = list()) 
{
  cl <- match.call()
  callstring <- paste(short.deparse(sys.call()), collapse = "")
  if (!inherits(formula, "formula")) 
    stop(paste("Argument", dQuote("formula"), "should be a formula"))
  stopifnot(is.hyperframe(data))
  data.sumry <- summary(data, brief = TRUE)
  npat <- data.sumry$ncases
  if (npat == 0) 
    stop(paste("Hyperframe", sQuote("data"), "has zero rows"))
  if (!is.null(iformula) && !inherits(iformula, "formula")) 
    stop(paste("Argument", sQuote("iformula"), "should be a formula or NULL"))
  if (has.random <- !is.null(random)) {
    if (!inherits(random, "formula")) 
      stop(paste(sQuote("random"), "should be a formula or NULL"))
    if (use.gam) 
      stop("Sorry, random effects are not available in GAMs")
  }
  if (!(is.interact(interaction) || is.hyperframe(interaction))) 
    stop(paste("The argument", sQuote("interaction"), "should be a point process interaction object (class", 
               dQuote("interact"), "), or a hyperframe containing such objects", 
               sep = ""))
  if (is.null(weights)) {
    weights <- rep(1, npat)
  }
  else {
    check.nvector(weights, npat, things = "rows of data", 
                  oneok = TRUE)
    if (length(weights) == 1L) 
      weights <- rep(weights, npat)
  }
  backdoor <- list(...)$backdoor
  if (is.null(backdoor) || !is.logical(backdoor)) 
    backdoor <- FALSE
  checkvars(formula, data.sumry$col.names, extra = c("x", 
                                                     "y", "id", "marks"), bname = "data")
  if (length(formula) < 3) 
    stop(paste("Argument", sQuote("formula"), "must have a left hand side"))
  Yname <- formula[[2]]
  trend <- formula[c(1, 3)]
  if (!is.name(Yname)) 
    stop("Left hand side of formula should be a single name")
  Yname <- paste(Yname)
  if (!inherits(trend, "formula")) 
    stop("Internal error: failed to extract RHS of formula")
  allvars <- variablesinformula(trend)
  itags <- if (is.hyperframe(interaction)) 
    names(interaction)
  else "Interaction"
  ninteract <- length(itags)
  if (is.null(iformula)) {
    if (ninteract > 1) 
      stop(paste("interaction hyperframe has more than 1 column;", 
                 "you must specify the choice of interaction", 
                 "using argument", sQuote("iformula")))
    iused <- TRUE
    iformula <- as.formula(paste("~", itags))
  }
  else {
    if (length(iformula) > 2) 
      stop(paste("The interaction formula", sQuote("iformula"), 
                 "should not have a left hand side"))
    permitted <- paste(sQuote("interaction"), "or permitted name in", 
                       sQuote("data"))
    checkvars(iformula, itags, extra = c(data.sumry$dfnames, 
                                         "id"), bname = permitted)
    ivars <- variablesinformula(iformula)
    iused <- itags %in% ivars
    if (sum(iused) == 0) 
      stop("No interaction specified in iformula")
    allvars <- c(allvars, ivars)
  }
  if (!is.null(random)) {
    if (length(random) > 2) 
      stop(paste("The random effects formula", sQuote("random"), 
                 "should not have a left hand side"))
    checkvars(random, itags, extra = c(data.sumry$col.names, 
                                       "x", "y", "id"), bname = "either data or interaction")
    allvars <- c(allvars, variablesinformula(random))
  }
  allvars <- unique(allvars)
  data <- cbind.hyperframe(data, id = factor(1:npat))
  data.sumry <- summary(data, brief = TRUE)
  allvars <- unique(c(allvars, "id"))
  Y <- data[, Yname, drop = TRUE]
  if (npat == 1) 
    Y <- solist(Y)
  Yclass <- data.sumry$classes[Yname]
  if (Yclass == "ppp") {
    Y <- solapply(Y, quadscheme, ...)
  }
  else if (Yclass == "quad") {
    Y <- as.solist(Y)
  }
  else {
    stop(paste("Column", dQuote(Yname), "of data", "does not consist of point patterns (class ppp)", 
               "nor of quadrature schemes (class quad)"), call. = FALSE)
  }
  datanames <- names(data)
  used.cov.names <- allvars[allvars %in% datanames]
  has.covar <- (length(used.cov.names) > 0)
  if (has.covar) {
    dfvar <- used.cov.names %in% data.sumry$dfnames
    imvar <- data.sumry$types[used.cov.names] == "im"
    if (any(nbg <- !(dfvar | imvar))) 
      stop(paste("Inappropriate format for", ngettext(sum(nbg), 
                                                      "covariate", "covariates"), paste(sQuote(used.cov.names[nbg]), 
                                                                                        collapse = ", "), ": should contain image objects or vector/factor"))
    covariates.hf <- data[, used.cov.names, drop = FALSE]
    has.design <- any(dfvar)
    dfvarnames <- used.cov.names[dfvar]
    datadf <- if (has.design) 
      as.data.frame(covariates.hf, discard = TRUE, warn = FALSE)
    else NULL
    if (has.design) {
      if (any(nbg <- matcolany(is.na(datadf)))) 
        stop(paste("There are NA's in the", ngettext(sum(nbg), 
                                                     "covariate", "covariates"), commasep(dQuote(names(datadf)[nbg]))))
    }
  }
  else {
    has.design <- FALSE
    datadf <- NULL
  }
  if (is.interact(interaction)) {
    ninteract <- 1
    processes <- list(Interaction = interaction$name)
    interaction <- hyperframe(Interaction = interaction, 
                              id = 1:npat)[, 1]
    constant <- c(Interaction = TRUE)
  }
  else if (is.hyperframe(interaction)) {
    inter.sumry <- summary(interaction)
    ninteract <- inter.sumry$nvars
    nr <- inter.sumry$ncases
    if (nr == 1 && npat > 1) {
      interaction <- cbind.hyperframe(id = 1:npat, interaction)[, 
                                                                -1]
      inter.sumry <- summary(interaction)
    }
    else if (nr != npat) 
      stop(paste("Number of rows in", sQuote("interaction"), 
                 "=", nr, "!=", npat, "=", "number of rows in", 
                 sQuote("data")))
    ok <- (inter.sumry$classes == "interact")
    if (!all(ok)) {
      nbg <- names(interaction)[!ok]
      nn <- sum(!ok)
      stop(paste(ngettext(nn, "Column", "Columns"), paste(sQuote(nbg), 
                                                          collapse = ", "), ngettext(nn, "does", "do"), 
                 "not consist of interaction objects"))
    }
    ok <- unlist(lapply(as.list(interaction), consistentname))
    if (!all(ok)) {
      nbg <- names(interaction)[!ok]
      stop(paste("Different interactions may not appear in a single column.", 
                 "Violated by", paste(sQuote(nbg), collapse = ", ")))
    }
    processes <- lapply(as.list(interaction), firstname)
    constant <- (inter.sumry$storage == "hyperatom")
    if (any(!constant)) {
      others <- interaction[, !constant]
      constant[!constant] <- sapply(lapply(as.list(others), 
                                           unique), length) == 1
    }
  }
  trivial <- unlist(lapply(as.list(interaction), allpoisson))
  nondfnames <- datanames[!(datanames %in% data.sumry$dfnames)]
  ip <- impliedpresence(itags, iformula, datadf, nondfnames)
  if (any(rowSums(ip) > 1)) 
    stop("iformula invokes more than one interaction on a single row")
  Vnamelist <- rep(list(character(0)), ninteract)
  names(Vnamelist) <- itags
  Isoffsetlist <- rep(list(logical(0)), ninteract)
  names(Isoffsetlist) <- itags
  for (i in 1:npat) {
    Yi <- Y[[i]]
    covariates <- if (has.covar) 
      covariates.hf[i, , drop = TRUE, strip = FALSE]
    else NULL
    # if (has.design) {
    #   covariates[dfvarname] <- lapply(as.list(as.data.frame(covariates[dfvarnames])), 
    #                                    as.im, W = Yi$data$window)
    # }
    prep0 <- bt.frame(Yi, trend, Poisson(), ..., covariates = NULL, 
                      allcovar = TRUE, use.gam = use.gam)
    glmdat <- prep0$glmdata
    for (j in (1:ninteract)[iused & !trivial]) {
      inter <- interaction[i, j, drop = TRUE]
      if (!is.null(ss <- inter$selfstart)) 
        interaction[i, j] <- inter <- ss(Yi$data, inter)
      prepj <- bt.frame(Yi, ~1, inter, ..., covariates = NULL, 
                        allcovar = TRUE, use.gam = use.gam, vnamebase = itags[j], 
                        vnameprefix = itags[j])
      vnameij <- prepj$Vnames
      if (i == 1) 
        Vnamelist[[j]] <- vnameij
      else if (!identical(vnameij, Vnamelist[[j]])) 
        stop("Internal error: Unexpected conflict in glm variable names")
      isoffset.ij <- prepj$IsOffset
      if (i == 1) 
        Isoffsetlist[[j]] <- isoffset.ij
      else if (!identical(isoffset.ij, Isoffsetlist[[j]])) 
        stop("Internal error: Unexpected conflict in offset indicators")
      glmdatj <- prepj$glmdata
      if (nrow(glmdatj) != nrow(glmdat)) 
        stop("Internal error: differing numbers of rows in glm data frame")
      iterms.ij <- glmdatj[vnameij]
      subset.ij <- glmdatj$.mpl.SUBSET
      glmdat <- cbind(glmdat, iterms.ij)
      glmdat$.mpl.SUBSET <- glmdat$.mpl.SUBSET & subset.ij
    }
    if (i == 1) {
      moadf <- glmdat
    }
    else {
      recognised <- names(glmdat) %in% names(moadf)
      if (any(!recognised)) {
        newnames <- names(glmdat)[!recognised]
        zeroes <- as.data.frame(matrix(0, nrow(moadf), 
                                       length(newnames)))
        names(zeroes) <- newnames
        moadf <- cbind(moadf, zeroes)
      }
      provided <- names(moadf) %in% names(glmdat)
      if (any(!provided)) {
        absentnames <- names(moadf)[!provided]
        zeroes <- as.data.frame(matrix(0, nrow(glmdat), 
                                       length(absentnames)))
        names(zeroes) <- absentnames
        glmdat <- cbind(glmdat, zeroes)
      }
      m.isfac <- sapply(as.list(glmdat), is.factor)
      g.isfac <- sapply(as.list(glmdat), is.factor)
      if (any(uhoh <- (m.isfac != g.isfac))) 
        errorInconsistentRows("values (factor and non-factor)", 
                              colnames(moadf)[uhoh])
      if (any(m.isfac)) {
        m.levels <- lapply(as.list(moadf)[m.isfac], 
                           levels)
        g.levels <- lapply(as.list(glmdat)[g.isfac], 
                           levels)
        clash <- !mapply(identical, x = m.levels, y = g.levels)
        if (any(clash)) 
          errorInconsistentRows("factor levels", (colnames(moadf)[m.isfac])[clash])
      }
      moadf <- rbind(moadf, glmdat)
    }
  }
  # moadf$caseweight <- weights[moadf$id]
  moadf$caseweight <- weights
  if (backdoor) 
    return(moadf)
  fmla <- prep0$trendfmla
  if (!all(trivial)) 
    fmla <- paste(fmla, "+", as.character(iformula)[[2]])
  fmla <- as.formula(fmla)
  for (j in (1:ninteract)[iused]) {
    vnames <- Vnamelist[[j]]
    tag <- itags[j]
    isoffset <- Isoffsetlist[[j]]
    if (any(isoffset)) {
      vnames[isoffset] <- paste("offset(", vnames[isoffset], 
                                ")", sep = "")
    }
    if (trivial[j]) 
      moadf[[tag]] <- 0
    else if (!identical(vnames, tag)) {
      if (length(vnames) == 1) 
        vn <- paste("~", vnames[1])
      else vn <- paste("~(", paste(vnames, collapse = " + "), 
                       ")")
      vnr <- as.formula(vn)[[2]]
      vnsub <- list(vnr)
      names(vnsub) <- tag
      fmla <- eval(substitute(substitute(fom, vnsub), 
                              list(fom = fmla)))
      if (has.random && tag %in% variablesinformula(random)) 
        random <- eval(substitute(substitute(fom, vnsub), 
                                  list(fom = random)))
    }
  }
  fmla <- as.formula(fmla)
  assign("glmmsubset", moadf$.mpl.SUBSET, envir = environment(fmla))
  for (nama in colnames(moadf)) assign(nama, moadf[[nama]], 
                                       envir = environment(fmla))
  glmmsubset <- .mpl.SUBSET <- moadf$.mpl.SUBSET
  .mpl.W <- moadf$.mpl.W
  caseweight <- moadf$caseweight
  want.trend <- prep0$info$want.trend
  if (want.trend && use.gam) {
    fitter <- "gam"
    ctrl <- do.call(gam.control, resolve.defaults(gcontrol, 
                                                  list(maxit = 50)))
    FIT <- gam(fmla, family = quasi(link = log, variance = mu), 
               weights = .mpl.W * caseweight, data = moadf, subset = (.mpl.SUBSET == 
                                                                        "TRUE"), control = ctrl)
    deviants <- deviance(FIT)
  }
  else if (!is.null(random)) {
    fitter <- "glmmPQL"
    ctrl <- do.call(lmeControl, resolve.defaults(gcontrol, 
                                                 list(maxIter = 50)))
    attr(fmla, "ctrl") <- ctrl
    fixed <- 42
    FIT <- hackglmmPQL(fmla, random = random, family = quasi(link = log, 
                                                             variance = mu), weights = .mpl.W * caseweight, data = moadf, 
                       subset = glmmsubset, control = attr(fixed, "ctrl"), 
                       reltol = reltol.pql)
    deviants <- -2 * logLik(FIT)
  }
  else {
    fitter <- "glm"
    ctrl <- do.call(glm.control, resolve.defaults(gcontrol, 
                                                  list(maxit = 50)))
    FIT <- glm(fmla, family = quasi(link = "log", variance = "mu"), 
               weights = .mpl.W * caseweight, data = moadf, subset = (.mpl.SUBSET == 
                                                                        "TRUE"), control = ctrl)
    deviants <- deviance(FIT)
  }
  env <- list2env(moadf, parent = sys.frame(sys.nframe()))
  environment(FIT$terms) <- env
  W <- with(moadf, .mpl.W * caseweight)
  SUBSET <- moadf$.mpl.SUBSET
  Z <- (moadf$.mpl.Y != 0)
  maxlogpl <- -(deviants/2 + sum(log(W[Z & SUBSET])) + sum(Z & 
                                                             SUBSET))
  result <- list(Call = list(callstring = callstring, cl = cl), 
                 Info = list(has.random = has.random, has.covar = has.covar, 
                             has.design = has.design, Yname = Yname, used.cov.names = used.cov.names, 
                             allvars = allvars, names.data = names(data), is.df.column = (data.sumry$storage == 
                                                                                            "dfcolumn"), rownames = row.names(data), correction = prep0$info$correction, 
                             rbord = prep0$info$rbord), Fit = list(fitter = fitter, 
                                                                   use.gam = use.gam, fmla = fmla, FIT = FIT, moadf = moadf, 
                                                                   Vnamelist = Vnamelist, Isoffsetlist = Isoffsetlist), 
                 Inter = list(ninteract = ninteract, interaction = interaction, 
                              iformula = iformula, iused = iused, itags = itags, 
                              processes = processes, trivial = trivial, constant = constant), 
                 formula = formula, trend = trend, iformula = iformula, 
                 random = random, npat = npat, data = data, Y = Y, maxlogpl = maxlogpl, 
                 datadf = datadf)
  class(result) <- c("mppm", class(result))
  return(result)
}

checkvars <- function(f, b, extra=NULL, bname=short.deparse(substitute(b))){
  fname <- short.deparse(substitute(f))
  fvars <- variablesinformula(f)
  bvars <- if(is.character(b)) b else names(b)
  bvars <- c(bvars, extra)
  nbg <- !(fvars %in% bvars)
  if(any(nbg)) {
    nn <- sum(nbg)
    stop(paste(ngettext(nn, "Variable", "Variables"),
               commasep(dQuote(fvars[nbg])),
               "in", fname,
               ngettext(nn, "is not one of the", "are not"),
               "names in", bname))
  }
  return(NULL)
}

consistentname <- function(x) {
  xnames <- unlist(lapply(x, getElement, name="name"))
  return(length(unique(xnames)) == 1)
}

firstname <- function(z) { z[[1]]$name }

allpoisson <- function(x) all(sapply(x, is.poisson.interact))

marklevels <- function(x) { levels(marks(x)) }

errorInconsistentRows <- function(what, offending) {
  stop(paste("There are inconsistent",
             what,
             "for the",
             ngettext(length(offending), "variable", "variables"),
             commasep(sQuote(offending)),
             "between different rows of the hyperframe 'data'"),
       call.=FALSE)
}