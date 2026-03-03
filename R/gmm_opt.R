#' @title GMM opt
#'export
gmm_opt <- function(theta,dat,mod=NULL,n_y,w_idx,iv_load_idx,iv_reg_idx,iv_names_load,y_names){
  if(!is.null(mod) && (!is.character(mod) || length(mod)!=1)){
    stop("mod must be NULL or a single model string")
  }

  ### theta 1, eta
  ### theta 2, n_y lambda
  ### theta n_y+1,..., 2(n_y)-1 var of epsilon

  Y <- as.matrix(dat[,1:n_y,drop=FALSE])
  W <- as.matrix(dat[,w_idx,drop=FALSE])
  IV_load <- as.matrix(dat[,iv_load_idx,drop=FALSE])
  IV_reg <- as.matrix(dat[,iv_reg_idx,drop=FALSE])
  n_iv_load <- ncol(IV_load)
  n_iv_reg <- ncol(IV_reg)
  n_w <- ncol(W)
  iv_names_norm <- make.names(iv_names_load)
  y_names_norm <- make.names(y_names)

  ## demean Y
  for (i in 1:n_y) {
    Y[,i] <- dat[,i]-mean(dat[,i])
  }

  n_load_mom <- 0
  for(j in 2:n_y){
    bad_j <- which(iv_names_norm %in% c(y_names_norm[1], y_names_norm[j]))
    n_load_mom <- n_load_mom + (n_iv_load - length(bad_j))
  }
  # loading-IV moments + cov(Y1,Yj) for j=2..n_y + variance moments + E[v] + E[v*IV_reg]
  tot_mom <- n_load_mom + (n_y - 1) + n_y + 1 + n_iv_reg

  mom <- vector("list", tot_mom )

  indx <- 1

  for (j in 2:n_y) {

    ### Moment conditions for lambda (IV estimation)
    rj <- Y[,j]-theta[j]*Y[,1]
    bad_j <- which(iv_names_norm %in% c(y_names_norm[1], y_names_norm[j]))
    use_k <- setdiff(seq_len(n_iv_load), bad_j)
    for(k in use_k){
      mom[[indx]] <- rj * IV_load[,k]
      indx <- indx + 1
    }
  }

  ### covariance moments: Cov(Y1, Yj) = lambda_j * Var(eta), for j = 2..n_y
  for (j in 2:n_y) {
    mom[[indx]] <- Y[,1]*Y[,j] - theta[j]*theta[1]
    indx <- indx +1
  }

  ### theta theta n_y+1,..., 2(n_y)-1 var of epsilon

  mom[[indx]] <- (Y[,1])^2-theta[1]- theta[n_y+1]
  indx <- indx +1

  for (h in 2:n_y) {
    mom[[indx]] <-(Y[,h])^2-(theta[h]^2)*theta[1] - theta[n_y+h]
    indx <- indx +1
  }


  ### optimal weight tau

  fake_lambda <- c(1,theta[2:n_y])
  w <- fake_lambda^2/theta[(n_y+1):(2*n_y)]
  omega <- w/sum(w)

  avgY <- rowSums( sweep( Y, 2, omega/fake_lambda, `*` ))

  # w_s <- sum(w)
  #Y[,1]*(w[1]/w_s) + (Y[,2]/2)*(w[2]/w_s) + (Y[,3]/3)*(w[3]/w_s)

  alpha <- theta[2*n_y+1]
  b_start <- 2*n_y+2
  b <- theta[b_start:(b_start+n_w-1)]
  v <- avgY - alpha - as.vector(W %*% b)
  mom[[indx]] <- v
  indx <- indx + 1
  for(k in 1:n_iv_reg){
    mom[[indx]] <- v * IV_reg[,k]
    indx <- indx + 1
  }

  do.call(cbind, mom)
}
