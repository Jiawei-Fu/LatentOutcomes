#' @title GMM moments for robust two-step estimation
#' @export
gmm_opt_robust <- function(theta, dat, n_y, w_idx, iv_load_idx, iv_reg_idx,
                           iv_names_load, iv_names_reg, y_names,
                           block = c("joint", "load", "reg"),
                           lambda_fix = NULL) {
  block <- match.arg(block)

  Y <- as.matrix(dat[, 1:n_y, drop = FALSE])
  W <- as.matrix(dat[, w_idx, drop = FALSE])
  IV_load <- as.matrix(dat[, iv_load_idx, drop = FALSE])
  IV_reg <- as.matrix(dat[, iv_reg_idx, drop = FALSE])

  n_iv_load <- ncol(IV_load)
  n_iv_reg <- ncol(IV_reg)
  n_w <- ncol(W)

  iv_names_norm <- make.names(iv_names_load)
  y_names_norm <- make.names(y_names)

  for (j in 1:n_y) {
    Y[, j] <- Y[, j] - mean(Y[, j])
  }

  if (block %in% c("joint", "load")) {
    lambda <- c(1, theta[1:(n_y - 1)])
  } else {
    if (is.null(lambda_fix) || length(lambda_fix) != n_y) {
      stop("lambda_fix must be provided with length n_y when block='reg'")
    }
    lambda <- as.numeric(lambda_fix)
  }

  if (block %in% c("joint", "reg")) {
    if (block == "joint") {
      alpha <- theta[n_y]
      b <- theta[(n_y + 1):(n_y + n_w)]
    } else {
      alpha <- theta[1]
      b <- theta[2:(1 + n_w)]
    }
  }

  mom <- list()
  idx <- 1

  if (block %in% c("joint", "load")) {
    for (j in 2:n_y) {
      rj <- Y[, j] - lambda[j] * Y[, 1]
      bad_j <- which(iv_names_norm %in% c(y_names_norm[1], y_names_norm[j]))
      use_k <- setdiff(seq_len(n_iv_load), bad_j)
      for (k in use_k) {
        mom[[idx]] <- rj * IV_load[, k]
        idx <- idx + 1
      }
    }
  }

  if (block %in% c("joint", "reg")) {
    wxb <- as.vector(W %*% b)
    for (j in 1:n_y) {
      eta_j <- Y[, j] / lambda[j]
      vj <- eta_j - alpha - wxb
      mom[[idx]] <- vj
      idx <- idx + 1
      for (k in seq_len(n_iv_reg)) {
        mom[[idx]] <- vj * IV_reg[, k]
        idx <- idx + 1
      }
    }
  }

  do.call(cbind, mom)
}
