#' @title GMM Equal
#'export
gmm_equal <- function(theta,dat,IV_Y,tau){

  ### theta 1,2, n_y-1 intercept
  ### theta  n_y...2*n_y-2 lambda 2,3,...,n_y-1

  ncol <- ncol(dat)
  Z <- dat[,ncol]
  Y <- dat[,1:(ncol-1)]
  n_y <- ncol-1 # how many outcome

  p <- mean(Z)

  tot_mom <- NA

  if(IV_Y==T){
    tot_mom <-  (n_y-1)*(n_y)
  }else{
    tot_mom <- (n_y-1)*2
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
    rj <- Y[,j]-theta[j-1]-theta[n_y+(j-2)]*Y[,1] # IV equ 1: intercept
    mom[[indx]] <- rj
    indx <- indx + 1
    mom[[indx]] <- rj*Z # IV equ 2: use Z as IV
    indx <- indx + 1

    if(IV_Y==T){ # use other Y as IV
      otherYs <- setdiff(2:n_y, j)
      for (k in otherYs) {
        mom[[indx]] <- rj * Y[, k]  # IV equ 2: use other Ys as IV
        indx <- indx + 1
      }
    }

  }
  ### equal weight tau

  avgY <- rowMeans(sweep(Y, 2, c(1,theta[n_y:(2*n_y-2)]), "/"))

  if(tau==F){
    mom[[indx]]<- (avgY - theta[2*n_y-1]) * Z
    indx <- indx + 1
    mom[[indx]] <- (avgY - theta[2*n_y]) * (1 - Z)
  }else{
    mom[[indx]]<- (Z/p-(1-Z)/(1-p))*avgY-theta[2*n_y-1]
  }

  do.call(cbind, mom)
}

