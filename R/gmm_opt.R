#' @title GMM opt
#'export
gmm_opt <- function(theta,dat,IV_Y,tau){

  ### theta 1, eta
  ### theta 2, n_y lambda
  ### theta n_y+1,..., 2(n_y)-1 var of epsilon

  ncol <- ncol(dat)
  Z <- dat[,ncol]
  n_y <- ncol-1 # how many outcome
  p <- mean(Z)
  Y <- dat[,1:(ncol-1)]

  ## demean Y
  for (i in 1:n_y) {
    Y[,i] <- dat[,i]-mean(dat[,i])
  }


  tot_mom <- NA

  if(IV_Y==T){
    tot_mom <-  (n_y-1)*(n_y-1)+1+n_y  # include var of eta, therefore +1; and var of epsilon
  }else{
    tot_mom <- n_y+n_y
  }

  if(tau==T){
    tot_mom <- tot_mom+1
  }else{
    tot_mom <- tot_mom+2
  }

  mom <- vector("list", tot_mom )

  indx <- 1

  for (j in 2:n_y) {

    ### Moment conditions for lambda (IV estimation)
    rj <- Y[,j]-theta[j]*Y[,1]
    mom[[indx]] <- rj*Z # IV equ 1: use Z as IV
    indx <- indx + 1

    if(IV_Y==T){ # use other Y as IV
      otherYs <- setdiff(2:n_y, j)
      for (k in otherYs) {
        mom[[indx]] <- rj * Y[, k]  # IV equ 2: use other Ys as IV
        indx <- indx + 1
      }
    }

  }

  ### theta 1 for variance of eta

  mom[[indx]] <- Y[,1]*Y[,2] - theta[2]*theta[1]
  indx <- indx +1

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

  if(tau==F){
    mom[[indx]] <- (avgY - theta[2*n_y+1]) * Z
    indx <- indx + 1
    mom[[indx]] <- (avgY - theta[2*n_y+2]) * (1-Z)
  }else{
    mom[[indx]]<- (Z/p-(1-Z)/(1-p))*avgY-theta[2*n_y+1]
  }

  do.call(cbind, mom)
}
