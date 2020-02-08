#' @title Spatial Variance Inflation Factor
#'
#' @description Calculate the VIF
#'
#' @param base_model Model to use as comparison
#' @param model Model to compare
#'
#' @return res data.frame with VIF for fixed parameters
#'
#' @export

SVIF <- function(base_model, model){
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
#' @param base_model Model to use as comparison
#' @param model Model to compare
#'
#' @return res data.frame with VRF for fixed parameters
#'
#' @export

SVRF <- function(base_model, model){
  res <- SVIF(base_model = base_model, model = model)
  res$VIF <- 1- 1/res$VIF

  names(res) <- c(c("parameter", "VRF"))

  return(res)
}

#' @title Projection matrix
#'
#' @description Calculate the projection matrix under several approaches
#'
#' @param X Covariate matrix
#' @param groups Ids for the subjects
#' @param method Projection's method
#'
#' @return Px Projection matrix
#' Px_ort (I - Px)

proj_mat <- function(X, groups = NULL, method = "rhz"){
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
#' @param W Adjcency matrix
#' @param sig Standard deviation
#'
#' @return x

ricar <- function(W, sig = 1){
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

#' @title Generate data from CAR model
#'
#' @description Generate data from CAR model
#'
#' @param W Adjcency matrix
#' @param sig Standard deviation
#' @param rho Dependence parameter
#'
#' @return x

rcar <- function(W, sig, rho = 0.9999){
  D <- diag(colSums(W))
  Q <- sig*(D - rho*W)
  sigma <- solve(Q)

  samp <- as.numeric(rmvnorm(n = 1, mean = rep(0, ncol(W)), sigma = sigma))
  return(samp)
}

#' @title select_marginal
#'
#' @description Select the desired marginals on a INLA models

select_marginal <- function(samp, ids){
  row_names <- row.names(samp$latent)
  samp <- c(samp$latent[row_names %in% ids])

  return(samp)
}

#' @title inla_stats
#'
#' @description Get inla statistics from a marginal

inla_stats <- function(marginal){
  z_stats <- inla.zmarginal(marginal = marginal, silent = TRUE)
  e_stats <- inla.mmarginal(marginal = marginal)
  kld_stat <- NA

  summary_df <- data.frame(z_stats$"mean", z_stats$"sd", z_stats$"quant0.025", z_stats$"quant0.5", z_stats$"quant0.75", e_stats, kld_stat)
  colnames(summary_df) <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant", "mode")

  return(summary_df)
}

#' @title Append two lists
#'
#' @description Get commom parameters in two list and generate one append list
#'
#' @param x List base
#' @param y Second list
#'
#' @return x

append_list <- function (x, y){
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

##-- Filters
meang <- function(x, g, weighted = FALSE){
  if(weighted){
    res <- tapply(X = x, INDEX = g, FUN = function(x) mean(x)*length(x))
  } else{
    res <- tapply(X = x, INDEX = g, FUN = mean)
  }

  return(res)
}

##-- Reduction ----
`%r%` <- function(x, g_index){

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

##-- Enlargement ----
`%e%` <- function(x, g_index){

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

##-- Updating INLA formula ----
update_inla_formula <- function(formula){
  ##-- Checking formula
  terms_formula <- terms.formula(formula, specials = c("f"), data = NULL)
  terms_labels <- attr(terms_formula, "term.labels")
  terms_f <- attr(terms_formula, "specials")$f - 1

  pos_restricted <- grep(x = terms_labels, pattern = "restricted")

  ##-- Updating formula
  if(length(pos_restricted) > 0){
    formula_char <- INLA:::inla.formula2character(formula)
    formula_char <- gsub(pattern = "restricted_besag", replacement = "besag", x = formula_char)
    formula_new <- as.formula(formula_char)

    terms_formula <- terms.formula(formula_new, specials = c("f"), data = NULL)
    terms_labels <- attr(terms_formula, "term.labels")
  } else{
    formula_new <- formula
  }

  if(length(terms_f) > 0){
    var_f <- list()
    for(i in seq_along(terms_f)){
      var_f[[i]] = eval(expr = parse(text = gsub(pattern = "^f\\(",
                                                 replacement = "INLA::f(",
                                                 x = terms_labels[terms_f[i]])))
    }

    ##-- Restricted components
    list_models <- lapply(var_f, "[", c("label", "model", "n"))
    list_restricted <- list_models[terms_f %in% pos_restricted]
    var_restricted <- unlist(lapply(list_restricted, FUN = "[[", "label"))
    size_restricted <- unlist(lapply(list_restricted, FUN = "[[", "n"))
  } else{
    var_restricted <- NULL
    size_restricted <- NULL
  }

  return(list(formula = formula_new, var_restricted = var_restricted, size_restricted = size_restricted))
}

#' @title Restricted Spatial Frailty Model in INLA
#'
#' @description Fit a Restricted Spatial Frailty model using INLA
#'
#' @param f INLA formula ?inla.surv
#' @param data data.frame
#' @param family 'exponential', 'weibull', 'weibullcure', 'loglogistic', 'gamma', 'lognormal' or 'pwe'
#' @param cure Desired cure fraction
#' @param proj 'none', 'rhz', 'hh' or 'spock'
#' @param fast To use the reduction operator
#' @param nsamp Sample size to use the projection approach
#' @param ... Other parameters used in ?inla
#'
#' @return INLA object with corrected parameters
#'
#' @import INLA
#' @importFrom stats density
#'
#' @export

rsfm_inla_old <- function(f, data, family, cure = FALSE, proj = "none", fast = TRUE, nsamp = 1000, ...){
  ##-- Time
  time_start <- Sys.time()

  ##-- Model
  if(cure & family == "weibull") family <- "weibullcure"

  family <- switch(family,
                   "exponential" = "exponential.surv",
                   "weibull" = "weibull.surv",
                   "weibullcure" = "weibullcure",
                   "loglogistic" = "loglogistic.surv",
                   "gamma" = "gamma.surv",
                   "lognormal" = "lognormal.surv",
                   "pwe" = "coxph",
                   stop(family, "family is not implemented."))

  ##-- Updating formula
  inla_formula <- update_inla_formula(formula = f)
  f <- inla_formula$formula
  reg_name <- inla_formula$var_restricted
  reg_size <- inla_formula$size_restricted

  time_start_inla <- Sys.time()
  mod <- inla(formula = f, data = data, family = family,
              control.compute = list(config = TRUE), ...)
  time_end_inla <- Sys.time()

  mod$summary.fixed$VRF <- 0

  out <- list()
  out$inla <- mod

  ##-- Correcting the model
  if(proj != "none" & length(reg_name) > 0){
    X <- as.matrix(mod$model.matrix)
    fixed_vars <- mod$names.fixed

    reg_name <- reg_name
    reg_size <- reg_size
    reg_pos <- data[[reg_name]]

    ##-- Getting a posterior sample from the model
    time_start_sampling <- Sys.time()
    model_sample <- inla.posterior.sample(result = mod, n = nsamp)
    time_end_sampling <- Sys.time()

    id_sample <- c(paste0(reg_name, ":", 1:reg_size), paste0(fixed_vars, ":", 1))

    model_sample <- do.call(args = lapply(model_sample, select_marginal, ids = id_sample), what = "rbind")

    W_samp <- model_sample[, 1:reg_size]
    colnames(W_samp) <- 1:reg_size
    beta_samp <- model_sample[, -c(1:(reg_size))]

    n_dens <- 100

    ##-- Projection matrices
    time_start_correction <- Sys.time()
    if(fast){
      proj_aux <- proj_mat(X = X, groups = reg_pos, method = proj)
    } else{
      proj_aux <- proj_mat(X = X, groups = NULL, method = proj)
    }

    ##-- Fixed effects
    if(fast){
      W_aux <- t(proj_aux$Paux%*%t(W_samp))
    } else{
      W_aux <- t(proj_aux$Paux%*%t(W_samp[, reg_pos]))
    }

    beta_ast <- matrix(beta_samp + W_aux, ncol = ncol(X))

    for(i in 1:length(fixed_vars)){
      fixed_var <- fixed_vars[i]

      dens <- density(beta_ast[, i], n = n_dens)
      dens_df <- data.frame(x = dens$x, y = dens$y)

      mod$marginals.fixed[[fixed_var]] <- dens_df
      summary_df <- suppressWarnings(inla_stats(marginal = dens_df))

      summary_df$VRF <- (out$inla$summary.fixed$sd[i]^2 - summary_df$sd^2)/(out$inla$summary.fixed$sd[i]^2)
      mod$summary.fixed[i, ] <- summary_df
    }

    ##-- Random effects
    if(fast){
      W_ast <- proj_aux$Px_ort%*%t(W_samp)
      Z_ast <- W_samp[, reg_pos] - t(W_ast[reg_pos, ])
    } else{
      W_ast <- proj_aux$Px_ort%*%t(W_samp[, reg_pos])
      Z_ast <- W_samp[, reg_pos] - t(W_ast)

      W_ast <- apply(X = W_ast, MARGIN = 2, meang, g = reg_pos)
    }

    ##-- ICAR
    for(i in 1:reg_size){
      dens <- density(W_ast[i, ], n = n_dens)
      dens_df <- data.frame(x = dens$x, y = dens$y)

      mod$marginals.random[[reg_name]][[sprintf("index.%s", i)]] <- dens_df

      summary_df <- suppressWarnings(inla_stats(marginal = dens_df))
      mod$summary.random[[reg_name]][i, -1] <- summary_df[1, ]
    }

    out[["restricted"]] <- mod

    vif_tab <- VIF(base_model = out[["inla"]], model = out[["restricted"]])
    out$VIF <- vif_tab

    ##-- Time
    time_end_correction <- Sys.time()
  } else{
    Z_ast <- NULL
    W_ast <- NULL
    time_start_correction <- time_start_sampling <- Sys.time()
    time_end_correction <- time_end_sampling <- time_start_correction
  }

  time_tab <- data.frame(INLA = as.numeric(difftime(time1 = time_end_inla, time2 = time_start_inla, units = "secs")),
                         sampling = as.numeric(difftime(time1 = time_end_sampling, time2 = time_start_sampling, units = "secs")),
                         correction = as.numeric(difftime(time1 = time_end_correction, time2 = time_start_correction, units = "secs")),
                         total = as.numeric(difftime(time1 = Sys.time(), time2 = time_start, units = "secs")))

  out$time <- time_tab

  return(out)
}

