#' @title GMM moments (per-outcome proxy, no inverse-variance aggregation)
#' @export
gmm_opt <- function(theta, dat, n_y, w_idx, iv_load_idx, iv_reg_idx,
                    iv_names_load, iv_names_reg, y_names) {
  Y <- as.matrix(dat[, 1:n_y, drop = FALSE])
  W <- as.matrix(dat[, w_idx, drop = FALSE])
  IV_load <- as.matrix(dat[, iv_load_idx, drop = FALSE])
  IV_reg <- as.matrix(dat[, iv_reg_idx, drop = FALSE])

  n_iv_load <- ncol(IV_load)
  n_iv_reg <- ncol(IV_reg)
  n_w <- ncol(W)

  iv_names_norm <- make.names(iv_names_load)
  iv_reg_names_norm <- make.names(iv_names_reg)
  y_names_norm <- make.names(y_names)
  overlap_norm <- intersect(iv_names_norm, iv_reg_names_norm)

  # Center outcomes to remove indicator means from loading moments.
  for (j in 1:n_y) {
    Y[, j] <- Y[, j] - mean(Y[, j])
  }

  # Parameter layout:
  # theta[1:(n_y-1)]      = lambda_2,...,lambda_ny
  # theta[n_y]            = alpha
  # theta[(n_y+1):(n_y+n_w)] = beta for W
  lambda <- c(1, theta[1:(n_y - 1)])
  alpha <- theta[n_y]
  b <- theta[(n_y + 1):(n_y + n_w)]

  # Count loading moments after excluding Y1 and Yj as IV for each lambda_j.
  n_load_mom <- 0
  for (j in 2:n_y) {
    bad_j <- which(iv_names_norm %in% c(y_names_norm[1], y_names_norm[j]))
    n_load_mom <- n_load_mom + (n_iv_load - length(bad_j))
  }

  # Moments:
  # 1) loading-IV moments
  # 2) regression moments for each rescaled proxy eta_j = Y_j / lambda_j:
  #    E[v_j]=0 and E[v_j * IV_reg_kept_j]=0.
  #    For j>=2, drop overlap moments where IV_reg variable name is in IV_load.
  n_reg_mom <- 0
  for (j in 1:n_y) {
    if (j == 1) {
      reg_keep <- seq_len(n_iv_reg)
    } else {
      reg_keep <- setdiff(seq_len(n_iv_reg), which(iv_reg_names_norm %in% overlap_norm))
    }
    n_reg_mom <- n_reg_mom + 1 + length(reg_keep)
  }
  tot_mom <- n_load_mom + n_reg_mom
  mom <- vector("list", tot_mom)
  idx <- 1

  # 1) Loading moments
  for (j in 2:n_y) {
    rj <- Y[, j] - lambda[j] * Y[, 1]
    bad_j <- which(iv_names_norm %in% c(y_names_norm[1], y_names_norm[j]))
    use_k <- setdiff(seq_len(n_iv_load), bad_j)
    for (k in use_k) {
      mom[[idx]] <- rj * IV_load[, k]
      idx <- idx + 1
    }
  }

  # 2) Regression moments
  wxb <- as.vector(W %*% b)
  for (j in 1:n_y) {
    eta_j <- Y[, j] / lambda[j]
    vj <- eta_j - alpha - wxb
    mom[[idx]] <- vj
    idx <- idx + 1
    if (j == 1) {
      reg_keep <- seq_len(n_iv_reg)
    } else {
      reg_keep <- setdiff(seq_len(n_iv_reg), which(iv_reg_names_norm %in% overlap_norm))
    }
    for (k in reg_keep) {
      mom[[idx]] <- vj * IV_reg[, k]
      idx <- idx + 1
    }
  }

  do.call(cbind, mom)
}
