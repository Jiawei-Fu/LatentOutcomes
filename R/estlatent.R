#' Estimate latent-outcome effects (per-outcome GMM without inverse-variance aggregation)
#'
#' @param mod Required regression specification for latent regression.
#'   Recommended form: "~trt+x1+x2". Also accepts "Y~..." or "eta~...".
#' @param Y Measured outcomes: matrix-data.frame with >=2 columns, or character vector of names when data is provided.
#' @param data Optional data frame to resolve character inputs and covariates in mod.
#' @param method "sem" or "gmm".
#' @param IV_Y Optional loading-IV names for GMM. Default uses all regressors from `mod`.
#' @param opt For GMM: TRUE=two-step efficient GMM, FALSE=one-step GMM.
#' @export
estlatent <- function(mod, Y, data = NULL, method = "sem", IV_Y = NULL, opt = TRUE) {

  if (!is.null(data)) {
    if (!is.data.frame(data)) stop("data must be a data.frame")
    if (is.character(Y)) {
      if (any(!Y %in% colnames(data))) stop("Some Y columns are not found in data")
      Y <- data[, Y, drop = FALSE]
    }
  }
  if (missing(mod)) stop("mod must be provided")
  if (is.null(mod)) stop("mod cannot be NULL")
  if (!is.data.frame(Y)) {
    Y <- as.data.frame(Y, check.names = FALSE)
  }

  if (is.null(colnames(Y))) colnames(Y) <- paste0("Y", seq_len(ncol(Y)))

  if (!is.logical(opt) || length(opt) != 1 || is.na(opt)) stop("opt must be TRUE or FALSE")
  if (!(method %in% c("sem", "gmm"))) stop("method must be one of 'sem' or 'gmm'")

  n_y <- ncol(Y)
  if (n_y <= 1) stop("The function needs more than 1 outcome measures")
  for (j in 1:n_y) {
    if (any(is.na(Y[, j]))) stop("Y has NAs; please remove or fill in NAs first")
  }

  cov_dat <- if (!is.null(data)) data else parent.frame()

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

  mm <- tryCatch(
    model.matrix(as.formula(paste0("~", rhs_text)), data = cov_dat),
    error = function(e) stop("Variables in mod must be available in `data` (or current environment when `data` is NULL).")
  )
  if ("(Intercept)" %in% colnames(mm)) mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  if (ncol(mm) == 0) stop("mod produced no predictors after removing intercept")
  colnames(mm) <- make.names(colnames(mm), unique = TRUE)
  cov_mm <- data.frame(mm, check.names = FALSE)

  w_names <- colnames(cov_mm)

  # loading-IV candidate pool: allow any variable from data, plus outcomes and mod covariates.
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
  if (is.null(IV_Y)) {
    iv_names_load <- w_names
  } else {
    if (!is.character(IV_Y)) stop("IV_Y must be a character vector of variable names")
    iv_names_load <- unique(IV_Y)
  }
  if (any(!iv_names_load %in% colnames(iv_pool))) {
    stop("Some IV names are not found in available variables (data, outcomes, and mod covariates)")
  }

  iv_load_mm <- iv_pool[, iv_names_load, drop = FALSE]
  iv_names_reg <- colnames(cov_mm)
  iv_reg_mm <- cov_mm[, iv_names_reg, drop = FALSE]

  y_names <- colnames(Y)
  n_iv_load <- ncol(iv_load_mm)
  n_iv_reg <- ncol(iv_reg_mm)
  n_w <- ncol(cov_mm)

  # SEM benchmark (for method='sem' output and GMM starting values)
  reg_text <- paste0("eta ~ ", paste(colnames(cov_mm), collapse = "+"))
  dat_sem <- cbind(cov_mm, Y)
  ld_text <- paste0("eta =~ 1*", paste0(colnames(Y), collapse = "+"))
  mod_c <- paste(ld_text, reg_text, sep = "\n")

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

  if (any(!w_names %in% names(rg_est))) stop("SEM result does not contain all regression coefficients from mod")
  sem_b_est <- as.numeric(rg_est[w_names])
  sem_b_se <- as.numeric(rg_se[w_names])
  sem_b_z <- as.numeric(rg_z[w_names])
  sem_b_p <- as.numeric(rg_p[w_names])

  if (method == "sem") {
    final_lambda_est <- sem_lambda_est
    final_lambda_se <- sem_lambda_se
    final_lambda_z <- sem_lambda_z
    final_lambda_p <- sem_lambda_p

    final_b_est <- sem_b_est
    final_b_se <- sem_b_se
    final_b_z <- sem_b_z
    final_b_p <- sem_b_p
  }

  if (method == "gmm") {
    dat_gmm <- cbind(Y, cov_mm, iv_load_mm, iv_reg_mm)
    w_idx <- (n_y + 1):(n_y + n_w)
    iv_load_idx <- (n_y + n_w + 1):(n_y + n_w + n_iv_load)
    iv_reg_idx <- (n_y + n_w + n_iv_load + 1):(n_y + n_w + n_iv_load + n_iv_reg)

    iv_names_norm <- make.names(iv_names_load)
    iv_reg_names_norm <- make.names(iv_names_reg)
    y_names_norm <- make.names(y_names)
    overlap_norm <- intersect(iv_names_norm, iv_reg_names_norm)

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

    # Parameters: lambda_2..lambda_ny + alpha + beta(W)
    n_par <- (n_y - 1) + 1 + n_w
    # Moments: loading moments + per-outcome regression moments.
    # For j>=2, remove v_j * W_k moments if W_k name overlaps loading IV names.
    n_reg_mom <- 0
    for (j in 1:n_y) {
      if (j == 1) {
        reg_keep <- seq_len(n_iv_reg)
      } else {
        reg_keep <- setdiff(seq_len(n_iv_reg), which(iv_reg_names_norm %in% overlap_norm))
      }
      n_reg_mom <- n_reg_mom + 1 + length(reg_keep)
    }
    n_mom <- n_load_mom + n_reg_mom
    if (n_mom < n_par) stop("gmm_opt under-identified: add more valid IVs in IV_Y")

    g_opt <- function(theta, dat_gmm) {
      gmm_opt(theta, dat_gmm,
              n_y = n_y,
              w_idx = w_idx,
              iv_load_idx = iv_load_idx,
              iv_reg_idx = iv_reg_idx,
              iv_names_load = iv_names_load,
              iv_names_reg = iv_names_reg,
              y_names = y_names)
    }

    sem_alpha <- pe$est[pe$lhs == "eta" & pe$op == "~1"][1]
    if (length(sem_alpha) == 0 || is.na(sem_alpha)) sem_alpha <- 0

    theta_s <- c(sem_lambda_est[2:n_y], sem_alpha, as.numeric(rg_est[w_names]))
    mom_model <- momentModel(g_opt, dat_gmm, theta0 = theta_s, vcov = "iid")
    fit_type <- if (opt) "twostep" else "onestep"
    fit_gmm <- gmmFit(mom_model, type = fit_type)

    cf <- coef(fit_gmm)
    vv <- vcov(fit_gmm)
    se <- sqrt(diag(vv))

    final_lambda_est <- c(1, cf[1:(n_y - 1)])
    final_lambda_se <- c(0, se[1:(n_y - 1)])
    final_lambda_z <- c(NA, final_lambda_est[-1] / final_lambda_se[-1])
    final_lambda_p <- 2 * (1 - pnorm(abs(final_lambda_z)))

    b_idx <- (n_y + 1):(n_y + n_w)
    final_b_est <- cf[b_idx]
    final_b_se <- se[b_idx]
    final_b_z <- final_b_est / final_b_se
    final_b_p <- 2 * (1 - pnorm(abs(final_b_z)))
  }

  lambda_name <- paste0("lambda_", seq_len(n_y))
  reg_name <- w_names
  coef_name <- c(lambda_name, reg_name)
  coef_est <- c(final_lambda_est, final_b_est)
  coef_se <- c(final_lambda_se, final_b_se)
  coef_z <- c(final_lambda_z, final_b_z)
  coef_p <- c(final_lambda_p, final_b_p)

  output <- cbind(Estimate = coef_est, SE = coef_se, Z = coef_z, `P-value` = coef_p)
  rownames(output) <- coef_name

  fit <- list(
    call = match.call(),
    method = method,
    mod = mod,
    IV = if (exists("iv_names_load")) iv_names_load else NULL,
    coefficients = output
  )
  class(fit) <- "estlatent"
  fit
}
