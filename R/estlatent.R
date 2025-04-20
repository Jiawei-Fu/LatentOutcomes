#' Estimate the Causal Effects with Latent Outcomes
#'
#'@param Z a vector of Treatment variable.
#'@param Y a data frame that includes the measured outcomes. The outcome in the first column will be set as lambda=1.
#'@param X an optional data frame that includes other covariates. The defauly is NULL.
#'@param method estimators options.
#'@param eta tbd
#'@param tau For GMM only. If not, calculate sample average in two groups
#'@param IV_Y Whether use other outcome measures as IVs
#'
#'@examples
#' \dontrun{
#' data(test_data)  # input data
#' estlatent(test_dat$Z,test_dat[,1:3],X=NULL,eta = 1,method="sem",IV_Y=T,tau=T)
#'}
#'@import lavaan
#'@import momentfit
#'@importFrom stats pnorm printCoefmat
#'@export

estlatent <- function(Z,Y,X=NULL,eta = 1,method="sem",IV_Y=T,tau=T){

  n_z <- ncol(Z)
  n_y <- ncol(Y)

  if(n_y==1){stop("The function needs more than 1 outcome measures")}

  if(is.data.frame(Z)==FALSE){Z <- data.frame(Z=Z)}

  ### recall deal with NA

  ### sem prepare
  if(is.null(X)){
    reg_text <- paste0("eta ~ ",colnames(Z))
    dat <- cbind(Z,Y)
  }else{
    reg_text <- paste0("eta ~ ",colnames(Z),"+", paste0(colnames(X),collapse = "+"))
    dat <- cbind(Z,Y,X)
  }

  ld_text <-  paste0("eta =~ 1*", paste0(colnames(Y),collapse = "+"))
  var_text <- paste0(colnames(Z),"~~",  colnames(Z))

  mod_c <- paste(ld_text,reg_text,var_text,sep="\n")

  sem_tmp <- sem(mod_c,data=dat)

  ### sem results

  sem_lambda_est <- summary(sem_tmp)$pe[1:n_y,5]
  sem_lambda_se <- summary(sem_tmp)$pe[1:n_y,6]
  sem_lambda_p <- summary(sem_tmp)$pe[1:n_y,8]
  sem_lambda_z <- summary(sem_tmp)$pe[1:n_y,7]
  sem_beta_est <- summary(sem_tmp)$pe[n_y+1,5]
  sem_beta_se <- summary(sem_tmp)$pe[n_y+1,6]
  sem_beta_p <- summary(sem_tmp)$pe[n_y+1,8]
  sem_beta_z <- summary(sem_tmp)$pe[n_y+1,7]
  sem_eta_var <- summary(sem_tmp)$pe[2*n_y+3,5]

  sem_var_y <- summary(sem_tmp)$pe[(n_y+3):(2*n_y+2),5]

  if(method=="sem"){
    final_lambda_est <- sem_lambda_est
    final_lambda_se <- sem_lambda_se
    final_lambda_p <- sem_lambda_p
    final_lambda_z <- sem_lambda_z
    final_beta_est <- sem_beta_est
    final_beta_se <- sem_beta_se
    final_beta_z <- sem_beta_z
    final_beta_p <- sem_beta_p
  }

  ### equal
  if(method=="gmm_equal"){
    IV_Y <- IV_Y
    tau <- tau
    dat <- cbind(Y,Z)
    g_equal <- function(theta, dat) {
      gmm_equal(theta, dat, IV_Y = IV_Y, tau = tau)
    }
    # starting value
    if(tau==T){
    theta_s <- c(rep(1,n_y-1),sem_lambda_est[2:n_y],sem_beta_est)
    res_equal <- momentModel(g_equal,dat,theta0=theta_s,vcov="iid")
    rec_gmm_equal <- gmmFit(res_equal)
    final_lambda_est <- c(1,coef(rec_gmm_equal)[n_y:(2*n_y-2)])
    final_lambda_se <- c(0,sqrt(diag(vcov(rec_gmm_equal))[n_y:(2*n_y-2)]))
    final_lambda_z <- c(NA,final_lambda_est[-1]/final_lambda_se[-1])
    final_lambda_p <- 2 * (1 - pnorm(abs(final_lambda_z)))   # two‑sided
    final_beta_est <- coef(rec_gmm_equal)[2*n_y-1]
    final_beta_se <- sqrt(diag(vcov(rec_gmm_equal))[2*n_y-1])
    final_beta_z <- final_beta_est/final_beta_se
    final_beta_p <- 2 * (1 - pnorm(abs(final_beta_z)))   # two‑sided

    }else{
      theta_s <- c(rep(1,n_y-1),sem_lambda_est[2:n_y],mean(dat[,1][Z==1]),mean(dat[,1][Z==0]))
      res_equal <- momentModel(g_equal,dat,theta0=theta_s,vcov="iid")
      rec_gmm_equal <- gmmFit(res_equal)
      final_lambda_est <- c(1,coef(rec_gmm_equal)[n_y:(2*n_y-2)])
      final_lambda_se <- c(0,sqrt(diag(vcov(rec_gmm_equal))[n_y:(2*n_y-2)]))
      final_lambda_z <- c(NA,final_lambda_est[-1]/final_lambda_se[-1])
      final_lambda_p <- 2 * (1 - pnorm(abs(final_lambda_z)))   # two‑sided
      final_beta_est <- coef(rec_gmm_equal)[2*n_y-1]-coef(rec_gmm_equal)[2*n_y]
      final_beta_se <- sqrt(diag(vcov(rec_gmm_equal))[2*n_y-1] + diag(vcov(rec_gmm_equal))[2*n_y] - 2*vcov(rec_gmm_equal)[2*n_y-1,2*n_y])
      final_beta_z <- final_beta_est/final_beta_se
      final_beta_p <- 2 * (1 - pnorm(abs(final_beta_z)))   # two‑sided
    }
  }


  ### optimal

  if(method=="gmm_opt"){
    IV_Y <- IV_Y
    tau <- tau
    dat <- cbind(Y,Z)
    g_opt <- function(theta, dat) {
      gmm_equal(theta, dat, IV_Y = IV_Y, tau = tau)
    }
    # starting value
    if(tau==T){
    theta_s <- c(sem_eta_var,sem_lambda_est[2:n_y],sem_var_y,sem_beta_est)
    res_opt <-momentModel(g_opt,dat,theta0=theta_s,vcov="iid")
    rec_gmm_opt <- gmmFit(res_opt)
    final_lambda_est <- c(1,coef(rec_gmm_opt)[2:n_y])
    final_lambda_se <- c(0,sqrt(diag(vcov(rec_gmm_opt))[2:n_y]))
    final_lambda_z <- c(NA,final_lambda_est[-1]/final_lambda_se[-1])
    final_lambda_p <- 2 * (1 - pnorm(abs(final_lambda_z)))   # two‑sided
    final_beta_est <- coef(rec_gmm_opt)[2*n_y+1]
    final_beta_se <- sqrt(diag(vcov(rec_gmm_opt))[2*n_y+1])
    final_beta_z <- final_beta_est/final_beta_se
    final_beta_p <- 2 * (1 - pnorm(abs(final_beta_z)))   # two‑sided

    }else{
      theta_s <- c(sem_eta_var,sem_lambda_est[2:n_y],sem_var_y,mean(dat[,1][Z==1]),mean(dat[,1][Z==0]))
      res_opt <-momentModel(g_opt,dat,theta0=theta_s,vcov="iid")
      rec_gmm_opt <- gmmFit(res_opt)
      final_lambda_est <- c(1,coef(rec_gmm_opt)[2:n_y])
      final_lambda_se <- c(0,sqrt(diag(vcov(rec_gmm_opt))[2:n_y]))
      final_lambda_z <- c(NA,final_lambda_est[-1]/final_lambda_se[-1])
      final_lambda_p <- 2 * (1 - pnorm(abs(final_lambda_z)))   # two‑sided
      final_beta_est <- coef(rec_gmm_opt)[2*n_y+1]-coef(rec_gmm_opt)[2*n_y+2]
      final_beta_se <- sqrt(diag(vcov(rec_gmm_opt))[2*n_y+1] + diag(vcov(rec_gmm_opt))[2*n_y+2] - 2*vcov(rec_gmm_opt)[2*n_y+1,2*n_y+2])
      final_beta_z <- final_beta_est/final_beta_se
      final_beta_p <- 2 * (1 - pnorm(abs(final_beta_z)))   # two‑sided
    }

  }


  ### ICW


  ### output

  lamdba_name <- rep(NA,n_y)

  for (j in 1:n_y) {
    lamdba_name[j] <- paste0("lambda_",j)
  }

  coef_name <- c(lamdba_name,"beta")
  coef_est <- c(final_lambda_est,final_beta_est)
  coef_se <- c(final_lambda_se,final_beta_se)
  coef_z <- c(final_lambda_z,final_beta_z)
  coef_p <- c(final_lambda_p,final_beta_p)

  output <- matrix(NA, nrow = n_y+1, ncol = 4)
  output[,1] <- coef_est
  output[,2] <- coef_se
  output[,3] <- coef_z
  output[,4] <- coef_p

  rownames(output) <- coef_name
  colnames(output) <- c("Estimate", "SE","Z","P-value")

  printCoefmat(output, P.values = TRUE, has.Pvalue = TRUE)

  # for extract
  output2 <- list()
  output2$estimate <- coef_est
  output2$se <- coef_se
  output2$p <- coef_p
  invisible(output2)

}
