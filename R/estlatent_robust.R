#' Estimate latent-outcome effects (two-step robust GMM with Godambe sandwich)
#'
#' @param mod Optional regression specification for latent regression.
#'   Recommended form: "~Z+x1+x2". Also accepts "Y~..." or "eta~...".
#'   If NULL, defaults to "~Z".
#' @param Z Treatment input: vector/1-col matrix-data.frame, or a single column name when data is provided.
#' @param Y Measured outcomes: matrix-data.frame with >=2 columns, or character vector of names when data is provided.
#' @param data Optional data frame to resolve character inputs and covariates in mod.
#' @param method "sem" or "gmm".
#' @param IV_Y Optional loading-IV names for stage-1 loading GMM. Default uses treatment Z only.
#' @param opt For each GMM stage: TRUE=two-step efficient GMM, FALSE=one-step GMM.
#'@examples
#' \dontrun{
#' data("test_data", package = "LatentOutcomes")  # input data
#' m1 <- estlatent_ave(mod="~Z", Z="Z", Y=c("Y1","Y2","Y3"), data=test_dat, method="sem")
#' summary(m1) # sem model
#' m2 <- estlatent_ave(mod="~Z", Z="Z", Y=c("Y1","Y2","Y3"), data=test_dat, method="gmm", IV_Y=c("Z","Y2"), opt=TRUE)
#' summary(m2) # gmm model, using Z and Y2 as IVs for the factor loadings
#' m3 <- estlatent_ave(mod="~Z+x1+x2", Z="Z", Y=c("Y1","Y2","Y3"), data=test_dat, method="gmm", IV_Y=c("Z","Y1","Y2","Y3"), opt=TRUE) # control covariates x1 and x2, using Z, Y1, Y2, and Y3 as IVs for the factor loadings
#' summary(m3)
#'}
#'@references Fu, Jiawei, and Donald P. Green. "Causal Inference for Experiments with Latent Outcomes: Key Results and Their Implications for Design and Analysis." (2025).
#'@import lavaan
#'@import momentfit
#'@importFrom stats pnorm printCoefmat
#'@export

estlatent_robust <- function(mod = NULL, Z, Y, data = NULL, method = "sem", IV_Y = NULL, opt = TRUE) {
  .pinv <- function(M, tol = 1e-8) {
    s <- svd(M)
    if (length(s$d) == 0) {
      return(matrix(0, ncol = nrow(M), nrow = ncol(M)))
    }
    dmax <- max(s$d)
    keep <- s$d > tol * dmax
    if (!any(keep)) {
      return(matrix(0, ncol = nrow(M), nrow = ncol(M)))
    }
    s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
  }

  .jacobian_fd <- function(fun, par, rel_step = 1e-6) {
    par <- as.numeric(par)
    f0 <- as.numeric(fun(par))
    m <- length(f0)
    p <- length(par)
    J <- matrix(0, nrow = m, ncol = p)
    for (k in seq_len(p)) {
      h <- rel_step * (abs(par[k]) + 1)
      p_plus <- par
      p_minus <- par
      p_plus[k] <- p_plus[k] + h
      p_minus[k] <- p_minus[k] - h
      f_plus <- as.numeric(fun(p_plus))
      f_minus <- as.numeric(fun(p_minus))
      J[, k] <- (f_plus - f_minus) / (2 * h)
    }
    J
  }

  .fit_stage_gmm <- function(moment_fun, theta0, opt) {
    gbar_fun <- function(theta) colMeans(moment_fun(theta))
    Q_fun <- function(theta, W) {
      g <- gbar_fun(theta)
      as.numeric(crossprod(g, W %*% g))
    }

    q <- length(gbar_fun(theta0))
    W1 <- diag(q)
    fit1 <- optim(theta0, fn = function(th) Q_fun(th, W1), method = "BFGS",
                  control = list(maxit = 5000, reltol = 1e-12))
    theta_pre <- as.numeric(fit1$par)

    if (!opt) {
      theta_hat <- theta_pre
      W_hat <- W1
    } else {
      m_pre <- moment_fun(theta_pre)
      S_pre <- crossprod(m_pre) / nrow(m_pre)
      W2 <- .pinv(S_pre)
      fit2 <- optim(theta_pre, fn = function(th) Q_fun(th, W2), method = "BFGS",
                    control = list(maxit = 5000, reltol = 1e-12))
      theta_hat <- as.numeric(fit2$par)
      W_hat <- W2
    }

    list(
      par = theta_hat,
      W = W_hat
    )
  }

  if (!is.null(data)) {
    if (!is.data.frame(data)) stop("data must be a data.frame")
    if (is.character(Z) && length(Z) == 1) {
      if (!(Z %in% colnames(data))) stop("Z not found in data")
      Z <- data[[Z]]
    }
    if (is.character(Y)) {
      if (any(!Y %in% colnames(data))) stop("Some Y columns are not found in data")
      Y <- data[, Y, drop = FALSE]
    }
  }

  if (!is.data.frame(Z)) {
    if (is.matrix(Z)) {
      Z <- as.data.frame(Z, check.names = FALSE)
    } else {
      Z <- data.frame(Z = Z)
    }
  }
  if (!is.data.frame(Y)) {
    Y <- as.data.frame(Y, check.names = FALSE)
  }

  if (is.null(colnames(Z))) colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
  if (is.null(colnames(Y))) colnames(Y) <- paste0("Y", seq_len(ncol(Y)))

  if (any(is.na(Z))) stop("Z has NAs; please remove or fill in NAs first")
  if (nrow(Z) != nrow(Y)) stop("Y and Z have different numbers")
  if (!is.logical(opt) || length(opt) != 1 || is.na(opt)) stop("opt must be TRUE or FALSE")
  if (!(method %in% c("sem", "gmm"))) stop("method must be one of 'sem' or 'gmm'")

  n_y <- ncol(Y)
  if (n_y <= 1) stop("The function needs more than 1 outcome measures")
  for (j in 1:n_y) {
    if (any(is.na(Y[, j]))) stop("Y has NAs; please remove or fill in NAs first")
  }

  cov_dat <- if (!is.null(data)) data else Z

  if (is.null(mod)) {
    rhs_text <- colnames(Z)[1]
  } else {
    if (!is.character(mod) || length(mod) != 1) stop("mod must be a single character string")
    mod_txt <- trimws(mod)
    if (startsWith(mod_txt, "~")) {
      rhs_text <- trimws(sub("^~", "", mod_txt))
    } else if (grepl("~", mod_txt, fixed = TRUE)) {
      parts <- strsplit(mod_txt, "~", fixed = TRUE)[[1]]
      if (length(parts) != 2) stop("mod must be in the form 'lhs ~ rhs'")
      lhs <- trimws(parts[1])
      rhs_text <- trimws(parts[2])
      if (nchar(lhs) > 0 && !(lhs %in% c("Y", "eta"))) {
        stop("mod lhs must be 'Y' or 'eta' (or omit lhs and use '~rhs')")
      }
    } else {
      rhs_text <- mod_txt
    }
    if (nchar(rhs_text) == 0) stop("mod rhs is empty")
  }

  mm <- tryCatch(
    model.matrix(as.formula(paste0("~", rhs_text)), data = cov_dat),
    error = function(e) stop("Variables in mod must be available in `data` (or in `Z` when `data` is NULL).")
  )
  if ("(Intercept)" %in% colnames(mm)) mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  if (ncol(mm) == 0) stop("mod produced no predictors after removing intercept")
  colnames(mm) <- make.names(colnames(mm), unique = TRUE)
  cov_mm <- data.frame(mm, check.names = FALSE)

  z_name <- make.names(colnames(Z)[1])
  if (!(z_name %in% colnames(cov_mm))) stop("mod must include the treatment variable (first column of Z)")

  x_names <- setdiff(colnames(cov_mm), z_name)
  n_x <- length(x_names)
  w_names <- colnames(cov_mm)

  iv_pool <- data.frame(check.names = FALSE)
  if (!is.null(data)) {
    iv_pool <- data
  }
  for (nm in colnames(cov_mm)) {
    if (!(nm %in% colnames(iv_pool))) iv_pool[[nm]] <- cov_mm[[nm]]
  }
  for (nm in colnames(Y)) {
    if (!(nm %in% colnames(iv_pool))) iv_pool[[nm]] <- Y[[nm]]
  }
  for (nm in colnames(Z)) {
    if (!(nm %in% colnames(iv_pool))) iv_pool[[nm]] <- Z[[nm]]
  }
  if (is.null(IV_Y)) {
    iv_names_load <- z_name
  } else {
    if (!is.character(IV_Y)) stop("IV_Y must be a character vector of variable names")
    iv_names_load <- unique(IV_Y)
  }
  if (any(!iv_names_load %in% colnames(iv_pool))) {
    stop("Some IV names are not found in available variables (data, outcomes, treatment, and mod covariates)")
  }

  iv_load_mm <- iv_pool[, iv_names_load, drop = FALSE]
  iv_names_reg <- colnames(cov_mm)
  iv_reg_mm <- cov_mm[, iv_names_reg, drop = FALSE]

  y_names <- colnames(Y)
  n_iv_load <- ncol(iv_load_mm)
  n_iv_reg <- ncol(iv_reg_mm)
  n_w <- ncol(cov_mm)

  reg_text <- paste0("eta ~ ", paste(colnames(cov_mm), collapse = "+"))
  dat_sem <- cbind(cov_mm, Y)
  ld_text <- paste0("eta =~ 1*", paste0(colnames(Y), collapse = "+"))
  var_text <- paste0(z_name, "~~", z_name)
  mod_c <- paste(ld_text, reg_text, var_text, sep = "\n")

  sem_tmp <- sem(mod_c, data = dat_sem)
  pe <- summary(sem_tmp)$pe

  ld <- pe[pe$lhs == "eta" & pe$op == "=~", ]
  ld_est <- setNames(ld$est, ld$rhs)
  ld_se <- setNames(ld$se, ld$rhs)
  ld_z <- setNames(ld$z, ld$rhs)
  ld_p <- setNames(ld$pvalue, ld$rhs)
  sem_lambda_est <- as.numeric(ld_est[colnames(Y)])
  sem_lambda_se <- as.numeric(ld_se[colnames(Y)])
  sem_lambda_z <- as.numeric(ld_z[colnames(Y)])
  sem_lambda_p <- as.numeric(ld_p[colnames(Y)])

  rg <- pe[pe$lhs == "eta" & pe$op == "~", ]
  rg_est <- setNames(rg$est, rg$rhs)
  rg_se <- setNames(rg$se, rg$rhs)
  rg_z <- setNames(rg$z, rg$rhs)
  rg_p <- setNames(rg$pvalue, rg$rhs)

  if (!(z_name %in% names(rg_est))) stop("SEM result does not contain treatment coefficient")
  sem_beta_est <- as.numeric(rg_est[z_name])
  sem_beta_se <- as.numeric(rg_se[z_name])
  sem_beta_z <- as.numeric(rg_z[z_name])
  sem_beta_p <- as.numeric(rg_p[z_name])

  if (n_x > 0) {
    sem_x_est <- as.numeric(rg_est[x_names])
    sem_x_se <- as.numeric(rg_se[x_names])
    sem_x_z <- as.numeric(rg_z[x_names])
    sem_x_p <- as.numeric(rg_p[x_names])
  }

  if (method == "sem") {
    final_lambda_est <- sem_lambda_est
    final_lambda_se <- sem_lambda_se
    final_lambda_z <- sem_lambda_z
    final_lambda_p <- sem_lambda_p

    final_beta_est <- sem_beta_est
    final_beta_se <- sem_beta_se
    final_beta_z <- sem_beta_z
    final_beta_p <- sem_beta_p

    if (n_x > 0) {
      final_x_est <- sem_x_est
      final_x_se <- sem_x_se
      final_x_z <- sem_x_z
      final_x_p <- sem_x_p
    }
  }

  if (method == "gmm") {
    dat_gmm <- cbind(Y, cov_mm, iv_load_mm, iv_reg_mm)
    w_idx <- (n_y + 1):(n_y + n_w)
    iv_load_idx <- (n_y + n_w + 1):(n_y + n_w + n_iv_load)
    iv_reg_idx <- (n_y + n_w + n_iv_load + 1):(n_y + n_w + n_iv_load + n_iv_reg)

    iv_names_norm <- make.names(iv_names_load)
    y_names_norm <- make.names(y_names)

    n_load_mom <- 0
    for (j in 2:n_y) {
      bad_j <- which(iv_names_norm %in% c(y_names_norm[1], y_names_norm[j]))
      usable_j <- n_iv_load - length(bad_j)
      if (usable_j <= 0) {
        stop(paste0("No valid loading IV remains for ", y_names[j],
                    "; IV_Y cannot include only Y1 and this outcome."))
      }
      n_load_mom <- n_load_mom + usable_j
    }
    if (n_load_mom < (n_y - 1)) {
      stop("Stage-1 loading GMM under-identified: add more valid IVs in IV_Y")
    }

    m1_fun <- function(theta_load) {
      gmm_opt_robust(theta_load, dat_gmm,
                      n_y = n_y,
                      w_idx = w_idx,
                      iv_load_idx = iv_load_idx,
                      iv_reg_idx = iv_reg_idx,
                      iv_names_load = iv_names_load,
                      iv_names_reg = iv_names_reg,
                      y_names = y_names,
                      block = "load")
    }

    sem_alpha <- pe$est[pe$lhs == "eta" & pe$op == "~1"][1]
    if (length(sem_alpha) == 0 || is.na(sem_alpha)) sem_alpha <- 0

    theta1_start <- as.numeric(sem_lambda_est[2:n_y])
    if (any(!is.finite(theta1_start))) theta1_start <- rep(1, n_y - 1)
    fit1 <- .fit_stage_gmm(m1_fun, theta1_start, opt = opt)
    theta1_hat <- as.numeric(fit1$par)
    lambda_hat <- c(1, theta1_hat)

    m2_fun <- function(theta_reg, lambda_val = lambda_hat) {
      gmm_opt_robust(theta_reg, dat_gmm,
                      n_y = n_y,
                      w_idx = w_idx,
                      iv_load_idx = iv_load_idx,
                      iv_reg_idx = iv_reg_idx,
                      iv_names_load = iv_names_load,
                      iv_names_reg = iv_names_reg,
                      y_names = y_names,
                      block = "reg",
                      lambda_fix = lambda_val)
    }

    theta2_start <- c(sem_alpha, as.numeric(rg_est[w_names]))
    theta2_start[!is.finite(theta2_start)] <- 0
    fit2 <- .fit_stage_gmm(function(th) m2_fun(th, lambda_hat), theta2_start, opt = opt)
    theta2_hat <- as.numeric(fit2$par)

    # Stacked estimating equations with fixed stage weight matrices:
    # s1 = G1' W1 g1bar, s2 = G2' W2 g2bar.
    W1 <- fit1$W
    W2 <- fit2$W
    p1 <- length(theta1_hat)
    p2 <- length(theta2_hat)
    p_all <- p1 + p2
    n <- nrow(dat_gmm)

    g1_bar_fun <- function(th1) colMeans(m1_fun(th1))
    g2_bar_fun <- function(th2, th1) colMeans(m2_fun(th2, c(1, th1)))

    G1_hat <- .jacobian_fd(g1_bar_fun, theta1_hat)
    G2_hat <- .jacobian_fd(function(th2) g2_bar_fun(th2, theta1_hat), theta2_hat)
    M1_hat <- m1_fun(theta1_hat)
    M2_hat <- m2_fun(theta2_hat, lambda_hat)

    Psi1 <- M1_hat %*% W1 %*% G1_hat
    Psi2 <- M2_hat %*% W2 %*% G2_hat
    Psi <- cbind(Psi1, Psi2)
    B_hat <- crossprod(Psi) / n

    s_stack_fun <- function(theta_all) {
      th1 <- theta_all[1:p1]
      th2 <- theta_all[(p1 + 1):p_all]

      g1_bar <- g1_bar_fun(th1)
      G1 <- .jacobian_fd(g1_bar_fun, th1)
      s1 <- as.numeric(t(G1) %*% W1 %*% g1_bar)

      g2_bar_local_fun <- function(t2) g2_bar_fun(t2, th1)
      g2_bar <- g2_bar_local_fun(th2)
      G2 <- .jacobian_fd(g2_bar_local_fun, th2)
      s2 <- as.numeric(t(G2) %*% W2 %*% g2_bar)

      c(s1, s2)
    }

    theta_all_hat <- c(theta1_hat, theta2_hat)
    A_hat <- .jacobian_fd(s_stack_fun, theta_all_hat)
    A_inv <- .pinv(A_hat)
    V_all <- A_inv %*% B_hat %*% t(A_inv) / n
    se_all <- sqrt(pmax(diag(V_all), 0))

    final_lambda_est <- c(1, theta1_hat)
    final_lambda_se <- c(0, se_all[1:p1])
    final_lambda_z <- c(NA, final_lambda_est[-1] / final_lambda_se[-1])
    final_lambda_p <- c(NA, 2 * (1 - pnorm(abs(final_lambda_z[-1]))))

    b_hat <- theta2_hat[2:(1 + n_w)]
    b_se <- se_all[(p1 + 2):p_all]
    b_idx <- which(w_names == z_name)
    final_beta_est <- b_hat[b_idx]
    final_beta_se <- b_se[b_idx]
    final_beta_z <- final_beta_est / final_beta_se
    final_beta_p <- 2 * (1 - pnorm(abs(final_beta_z)))

    if (n_x > 0) {
      x_idx <- match(x_names, w_names)
      final_x_est <- b_hat[x_idx]
      final_x_se <- b_se[x_idx]
      final_x_z <- final_x_est / final_x_se
      final_x_p <- 2 * (1 - pnorm(abs(final_x_z)))
    }
  }

  lambda_name <- paste0("lambda_", seq_len(n_y))
  trt_name <- colnames(Z)[1]

  if (n_x > 0) {
    coef_name <- c(lambda_name, trt_name, x_names)
    coef_est <- c(final_lambda_est, final_beta_est, final_x_est)
    coef_se <- c(final_lambda_se, final_beta_se, final_x_se)
    coef_z <- c(final_lambda_z, final_beta_z, final_x_z)
    coef_p <- c(final_lambda_p, final_beta_p, final_x_p)
  } else {
    coef_name <- c(lambda_name, trt_name)
    coef_est <- c(final_lambda_est, final_beta_est)
    coef_se <- c(final_lambda_se, final_beta_se)
    coef_z <- c(final_lambda_z, final_beta_z)
    coef_p <- c(final_lambda_p, final_beta_p)
  }

  output <- cbind(Estimate = coef_est, SE = coef_se, Z = coef_z, `P-value` = coef_p)
  rownames(output) <- coef_name

  fit <- list(
    call = match.call(),
    method = if (method == "gmm") "gmm_robust" else method,
    mod = mod,
    IV = if (exists("iv_names_load")) iv_names_load else NULL,
    coefficients = output
  )
  class(fit) <- "estlatent"
  fit
}
