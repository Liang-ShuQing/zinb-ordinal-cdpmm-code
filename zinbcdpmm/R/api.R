
# Public API wrappers and S3 methods

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Simulate longitudinal ZINB + ordinal data
#'
#' Thin wrapper around \code{generate_data}. Returns an object of class
#' \code{zinbcdpmm_data}.
#'
#' @param n Number of subjects.
#' @param nis Visits per subject (scalar or length-n vector).
#' @param seed RNG seed.
#' @param random_nis If TRUE, sample visit counts from \code{nis_range}.
#' @param nis_range Integer range for random visit counts.
#' @param scenario Simulation scenario (1 or 2).
#' @param C Number of ordinal categories.
#' @param re_dist Random-effects truth: "normal", "mixture", or "mvt".
#' @param mvt_df Degrees of freedom when \code{re_dist = "mvt"}.
#' @param rho Common RE correlation (default 0.5).
#' @param zi_intercept Optional override of zero-inflation intercept.
#' @param ... Passed to \code{generate_data}.
#' @return A list with class \code{zinbcdpmm_data}.
#' @export
simulate_zinb_ordinal <- function(n = 100, nis = 10, seed = 2026,
                                  random_nis = FALSE, nis_range = 1:20,
                                  scenario = 1, C = 5,
                                  re_dist = "normal",
                                  mvt_df = 4,
                                  rho = NULL,
                                  zi_intercept = NULL,
                                  ...) {
  dat <- generate_data(
    n = n, nis = nis, seed = seed,
    random_nis = random_nis, nis_range = nis_range,
    scenario = scenario, C = C,
    re_dist = re_dist, mvt_df = mvt_df,
    rho = rho, zi_intercept = zi_intercept, ...
  )
  class(dat) <- c("zinbcdpmm_data", class(dat))
  dat
}

#' Fit joint ZINB-ordinal model with CDPMM random effects
#'
#' Thin wrapper around \code{fit_joint_model}.
#'
#' @param data A list / \code{zinbcdpmm_data} from \code{simulate_zinb_ordinal}
#'   or with the same fields as \code{generate_data} output.
#' @param chain MCMC iterations.
#' @param burn Burn-in.
#' @param thin Thinning interval.
#' @param G Truncation level of the stick-breaking CDPMM (default 8).
#' @param ... Additional arguments passed to \code{fit_joint_model}
#'   (e.g. \code{iw_nu0}, \code{kappa0}, \code{r_update}).
#' @return An object of class \code{zinbcdpmm_fit}.
#' @export
fit_zinb_ordinal <- function(data, chain = 2000, burn = 1000, thin = 5,
                             G = 8, ...) {
  cl <- match.call()
  fit <- fit_joint_model(
    dat = data, chain = chain, burn = burn, thin = thin, G = G, ...
  )
  fit$call <- cl
  fit$n <- data$n
  fit$N <- data$N
  fit$chain <- chain
  fit$burn <- burn
  fit$thin <- thin
  fit$G <- G
  class(fit) <- c("zinbcdpmm_fit", class(fit))
  fit
}

#' @export
print.zinbcdpmm_fit <- function(x, ...) {
  cat("zinbcdpmm joint CDPMM fit\n")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
  }
  n_keep <- if (!is.null(x$alpha_samples)) nrow(x$alpha_samples) else NA_integer_
  cat(sprintf(
    "Subjects n=%s, observations N=%s, kept draws=%s\n",
    as.character(x$n %||% NA), as.character(x$N %||% NA), as.character(n_keep)
  ))
  if (!is.null(x$dic)) cat(sprintf("DIC=%.2f  WAIC=%.2f  LOOIC=%.2f\n",
                                   x$dic, x$waic, x$looic))
  invisible(x)
}

#' @export
summary.zinbcdpmm_fit <- function(object, ...) {
  sm <- list(
    call = object$call,
    alpha = cbind(
      mean = object$alpha_est,
      `2.5%` = object$alpha_ci_lower,
      `97.5%` = object$alpha_ci_upper
    ),
    beta = cbind(
      mean = object$beta_est,
      `2.5%` = object$beta_ci_lower,
      `97.5%` = object$beta_ci_upper
    ),
    gamma = cbind(
      mean = object$gamma_est,
      `2.5%` = object$gamma_ci_lower,
      `97.5%` = object$gamma_ci_upper
    ),
    delta = cbind(
      mean = object$delta_est,
      `2.5%` = object$delta_ci_lower,
      `97.5%` = object$delta_ci_upper
    ),
    r = c(
      mean = object$r_est,
      `2.5%` = unname(object$r_ci_lower),
      `97.5%` = unname(object$r_ci_upper)
    ),
    Sigma = object$Sigma_est,
    Rho = object$Rho_est,
    dic = object$dic,
    waic = object$waic,
    looic = object$looic,
    nclust_est = object$nclust_est
  )
  class(sm) <- "summary.zinbcdpmm_fit"
  sm
}

#' @export
print.summary.zinbcdpmm_fit <- function(x, digits = 3, ...) {
  cat("Summary of zinbcdpmm joint CDPMM fit\n\n")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  cat("Fixed effects (zero / alpha):\n")
  print(round(x$alpha, digits))
  cat("\nFixed effects (count / beta):\n")
  print(round(x$beta, digits))
  cat("\nFixed effects (ordinal / gamma):\n")
  print(round(x$gamma, digits))
  cat("\nThresholds (delta):\n")
  print(round(x$delta, digits))
  cat("\nNB size r:\n")
  print(round(x$r, digits))
  cat("\nImplied Sigma:\n")
  print(round(x$Sigma, digits))
  cat("\nRho (zero-count, zero-ord, count-ord):\n")
  print(round(x$Rho, digits))
  cat(sprintf("\nDIC=%.2f  WAIC=%.2f  LOOIC=%.2f  mean nclust=%.2f\n",
              x$dic, x$waic, x$looic, x$nclust_est))
  invisible(x)
}

