#' @title Spatial Variance Inflation Factor
#'
#' @description Calculate the VIF
#'
#' @param base_model model to use as comparison.
#' @param model model to compare.
#'
#' @return res data.frame with VIF for fixed parameters
#'
#' @export

SVIF <- function(base_model, model) {
  fixed_names <- rownames(base_model$summary_fixed)

  vif <- (model$summary_fixed$sd^2)/(base_model$summary_fixed$sd^2)
  vif <- data.frame(vif)

  res <- data.frame(fixed_names, vif, stringsAsFactors = FALSE)
  names(res) <- c("parameter", "VIF")

  return(res)
}

#' @title Spatial Variance Retraction Factor
#'
#' @description Calculate the VRF
#'
#' @param base_model model to use as comparison
#' @param model model to compare
#'
#' @return res data.frame with VRF for fixed parameters
#'
#' @export

SVRF <- function(base_model, model) {
  res <- SVIF(base_model = base_model, model = model)
  res$VIF <- 1- 1/res$VIF

  names(res) <- c(c("parameter", "VRF"))

  return(res)
}

#' @title Projection matrix
#'
#' @description Calculate the projection matrix under several approaches
#'
#' @param X covariate matrix
#' @param groups ids for the subjects
#' @param method projection's method
#'
#' @return Px Projection matrix
#' Px_ort (I - Px)

proj_mat <- function(X, groups = NULL, method = "rhz") {
  N <- nrow(X)
  if(is.null(groups)) groups <- 1:N

  n_groups <- length(unique(groups))

  ##-- Projection matrices ----
  if(method == "rhz"){
    if(n_groups != N){
      n_group <- diag(1/tapply(groups, groups, length))

      XX_inv <- solve(t(X)%*%X)
      Paux <- XX_inv%*%(t(X) %r% groups)
      Px <- (X %r% groups)%*%Paux
      Px_ort <- diag(1, nrow = n_groups, ncol = n_groups) - n_group%*%Px
    } else{
      XX_inv <- solve(t(X)%*%X)
      Paux <- XX_inv%*%t(X)
      Px <- X%*%Paux
      Px_ort <- diag(nrow(X)) - Px
    }
  }

  return(list(Paux = Paux, Px = Px, Px_ort = Px_ort))
}

#' @title Generate data from ICAR model
#'
#' @description Generate data from ICAR model
#'
#' @param W adjcency matrix
#' @param sig standard deviation
#'
#' @importFrom stats rnorm
#'
#' @return x

ricar <- function(W, sig = 1) {
  n <- ncol(W)
  num <- rowSums(W)

  Q <- -W
  diag(Q) <- num

  Q_aux <- eigen(Q)$vectors[, order(eigen(Q)$values)]

  D_aux <- sort(eigen(Q)$values)

  rnd <- rnorm(n-1, 0, sqrt(sig*(1/D_aux[-1])))
  rnd <- Q_aux%*%c(0, rnd)

  return(as.vector(rnd))
}

#' #' @title Generate data from CAR model
#' #'
#' #' @description Generate data from CAR model
#' #'
#' #' @param W adjcency matrix
#' #' @param sig standard deviation
#' #' @param rho dependence parameter
#' #'
#' #' @return x
#'
#' rcar <- function(W, sig, rho = 0.9999){
#'   D <- diag(colSums(W))
#'   Q <- sig*(D - rho*W)
#'   sigma <- solve(Q)
#'
#'   samp <- as.numeric(rmvnorm(n = 1, mean = rep(0, ncol(W)), sigma = sigma))
#'   return(samp)
#' }

#' @title select_marginal
#'
#' @description Select the desired marginals on a INLA models
#'
#' @param samp a sample from ?inla.posterior.sample
#' @param ids ids to restore

select_marginal <- function(samp, ids) {
  row_names <- row.names(samp$latent)
  samp <- c(samp$latent[row_names %in% ids])

  return(samp)
}

#' @title Append two lists
#'
#' @description Get commom parameters in two list and generate one append list
#'
#' @param x list base
#' @param y second list
#'
#' @return x

append_list <- function (x, y) {
  xnames <- names(x)
  for (v in names(y)) {
    if(v %in% xnames && is.list(x[[v]]) && is.list(y[[v]])){
      x[[v]] <- append_list(x[[v]], y[[v]])
    } else{
      if(!is.null(y[[v]])){
        x[[v]] <- y[[v]]
      }
    }
  }
  return(x)
}

#' @title meang
#'
#' @description Mean by groups
#'
#' @param x a numeric vector
#' @param g group indexes
#' @param weighted TRUE for weighted mean

meang <- function(x, g, weighted = FALSE) {
  if(weighted){
    res <- tapply(X = x, INDEX = g, FUN = function(x) mean(x)*length(x))
  } else{
    res <- tapply(X = x, INDEX = g, FUN = mean)
  }

  return(res)
}

#' @title Reduction operator
#'
#' @description Reduction operator
#'
#' @param x a numeric vector or a numeric matrix
#' @param g_index group indexes
#'
#' @export

`%r%` <- function(x, g_index) {

  if(!(class(x) %in% c("numeric", "matrix"))) stop("x must be a vector or a matrix")

  if(is.matrix(x)){
    n <- nrow(x)
    p <- ncol(x)

    if(!(length(g_index) %in% c(n, p))) stop("g_index should be equal to nrow(x) or ncol(x)")
    dim_reduction <- ifelse(length(g_index) == p, 1, 2)

    reduction <- apply(X = x, MARGIN = dim_reduction, FUN = function(y) tapply(X = y, INDEX = g_index, FUN = sum))

    if(dim_reduction == 1) reduction <- t(reduction)
  } else{
    reduction <- tapply(X = x, INDEX = g_index, FUN = sum)
  }

  return(reduction)
}

#' @title Enlargement operator
#'
#' @description Enlargement operator
#'
#' @param x a numeric vector or a numeric matrix
#' @param g_index group indexes
#'
#' @export

`%e%` <- function(x, g_index) {

  if(!(class(x) %in% c("numeric", "matrix"))) stop("x must be a vector or a matrix")

  if(is.matrix(x)) if(is.null(colnames(x)) & is.null(row.names(x))) stop("x must be a named matrix")
  if(!is.matrix(x)) if(is.null(names(x))) stop("x must be a named vector.")

  if(!all(names(x) %in% unique(g_index))) stop("names(x) and g_index does not match.")

  if(is.matrix(x)){
    n <- nrow(x)
    p <- ncol(x)
    ng <- length(unique(g_index))

    if(!(ng %in% c(n, p))) stop("Number of g_index groups should be equal to nrow(x) or ncol(x)")

    if(ng == p){
      dim_enlargement <- 1
      names_x <- colnames(x)
    } else{
      dim_enlargement <- 2
      names_x <- row.names(x)
    }

    n_group <- table(g_index)[names_x]

    enlargement <- apply(X = x, MARGIN = dim_enlargement, FUN = function(x) (x/n_group)[g_index])
    if(dim_enlargement == 1) enlargement <- t(enlargement)
  } else{
    names_x <- names(x)
    n_group <- table(g_index)[names_x]
    enlargement <- (x/table(g_index)[names(x)])[g_index]

    enlargement <- as.numeric(enlargement)
    names(enlargement) <- g_index
  }

  return(enlargement)
}

#' @title Updating INLA formula
#'
#' @description Updating INLA formula
#'
#' @param formula a formula to be updated to INLA format
#'
#' @importFrom stats terms.formula

update_inla_formula <- function(formula) {
  ##-- Checking formula
  terms_formula <- terms.formula(formula, specials = c("f"), data = NULL)
  terms_labels <- paste(attr(terms_formula, "variables"))
  terms_f <- attr(terms_formula, "specials")$f + 1 ##-- + 1 for the list parameter

  pos_restricted <- grep(x = terms_labels, pattern = "restricted|r_")
  pos_unrestricted <- grep(x = terms_labels,
                           pattern = "\"(iid|besag|besag2|besagproper|besagproper2)")

  ##-- Updating formula
  if(length(pos_restricted) > 0){
    formula_char <- format(formula)
    formula_char <- gsub(pattern = "restricted_besag|r_besag", replacement = "besag", x = formula_char)
    formula_new <- as.formula(formula_char)

    terms_formula <- terms.formula(formula_new, specials = c("f"), data = NULL)
    terms_labels <- paste(attr(terms_formula, "variables"))
  } else{
    formula_new <- formula
  }

  if(length(terms_f) > 0){
    var_f <- list()
    for(i in seq_along(terms_f)){
      var_f[[i]] = eval(expr = parse(text = gsub(pattern = "^f\\(",
                                                 replacement = "INLA::f(",
                                                 x = terms_labels[terms_f[i]])),
                        envir = parent.frame(n = 2))
    }

    ##-- Restricted and unrestricted components
    list_models <- lapply(var_f, "[", c("label", "model", "n"))

    list_restricted <- list_models[terms_f %in% pos_restricted]
    var_restricted <- unlist(lapply(list_restricted, FUN = "[[", "label"))
    size_restricted <- unlist(lapply(list_restricted, FUN = "[[", "n"))

    list_unrestricted <- list_models[terms_f %in% pos_unrestricted]
    var_unrestricted <- unlist(lapply(list_unrestricted, FUN = "[[", "label"))
    size_unrestricted <- unlist(lapply(list_unrestricted, FUN = "[[", "n"))
  } else{
    var_restricted <- NULL
    size_restricted <- NULL
    var_unrestricted <- NULL
    size_unrestricted <- NULL
  }

  return(list(formula = formula_new,
              var_restricted = var_restricted,
              vars_unrestricted = var_unrestricted,
              size_restricted = size_restricted,
              size_unrestricted = size_unrestricted))
}

#' @title Deviance Information Criterion
#'
#' @description Get the Deviance Information Criterion (DIC) from a model
#'
#' @param object a object from ?rsglmm, ?rscm or ?rsfm
#'
#' @return DIC
#'
#' @export

DIC <- function(object) {
  out <- object$out

  if(class(out) == "inla") {
    return(out$dic$dic)
  }

  if(class(out) == "sparse.sglmm") {
    return(out$dic)
  }

  stop(sprintf("Don't know how to deal with an object of class %s. Did you fit a model using rsglmm, rscm or rsfm?", class(out)))
}

#' @title Watanabe–Akaike information criterion
#'
#' @description Get the Watanabe–Akaike information criterion (WAIC) from a model
#'
#' @param object a object from ?rsglmm, ?rscm or ?rsfm
#'
#' @return WAIC
#'
#' @export

WAIC <- function(object) {
  out <- object$out

  if(class(out) == "inla") {
    return(out$waic$waic)
  }

  if(class(out) == "sparse.sglmm") {
    return(NA_real_)
  }

  stop(sprintf("Don't know how to deal with an object of class %s. Did you fit a model using rsglmm, rscm or rsfm?", class(out)))
}

## Copied from R2OpenBugs
findOpenBUGS <- function(){
  dir <- Sys.getenv("OpenBUGS_PATH")
  if(nchar(dir))
    return(dir)

  dir <- getOption("R2OpenBUGS.pgm")
  if(!is.null(dir))
    return(dir)

  if(.Platform$OS.type != "windows")
    return(Sys.which('OpenBUGS'))

  deps <- utils::packageDescription("R2OpenBUGS", fields="SystemRequirements")
  version.req <- gsub(".*OpenBUGS ?\\(>= ?(.+)\\).*", "\\1", deps)
  ob.reg <- try(utils::readRegistry("Software\\OpenBUGS", "HLM", view = "32-bit"), silent = TRUE)
  if (inherits(ob.reg, "try-error"))
    return(NA)

  rnames <- names(ob.reg)
  ver <- gsub("OpenBUGS ", "", rnames)
  version.inst <- gsub("(.+)e$","\\1", ver)

  if(length(version.inst > 1)){
    id <- which(apply(outer(version.inst, version.inst, Vectorize(compareVersion, c("a", "b"))), 1, function(x) all(x >= 0)))
    version.inst <- version.inst[id]
    rnames <- rnames[id]
  }

  if (compareVersion(max(version.inst), version.req) < 0)
    warning("Found OpenBUGS version ", version.inst, ".\n Requires ", version.req, " or greater.")

  utils::readRegistry(paste("Software", "OpenBUGS", rnames, sep="\\"), "HLM", view = "32-bit")[["InstallPath"]]
}
