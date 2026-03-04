#' Estimate the Causal Effects with Latent Outcomes
#'
#' @param mod Optional regression specification for the latent regression.
#'   Recommended form is RHS-only, e.g. `"~Z+x1+x2"`. Also accepts `"Y~..."` or `"eta~..."`.
#'   If `NULL`, the default model uses all columns of `Z` on the RHS.
#' @param Z Treatment input. Can be:
#'   a numeric/logical vector,
#'   a data frame/matrix with one or more columns,
#'   or character column name(s) when `data` is provided.
#' @param Y Measured outcomes. Can be:
#'   a data frame/matrix with 2+ outcome columns,
#'   or a character vector of outcome column names when `data` is provided.
#'   The first outcome is the reference loading (`lambda_1 = 1`).
#' @param data Optional data frame used to resolve character inputs in `Z` and `Y`,
#'   and to supply covariates referenced in `mod`.
#' @param method Estimation method: `"sem"` or `"gmm"` (default `"sem"`).
#' @param IV_Y Optional character vector of IV names for loading moments in GMM.
#'   Names can come from `data` (when provided), outcomes, treatment, and/or covariates implied by `mod`.
#'   Default is all treatment columns in `Z`. Regression moments always use covariates from `mod` as their own IVs.
#' @param opt Logical. For `method="gmm"`, `TRUE` (default) uses efficient two-step GMM;
#'   `FALSE` uses one-step GMM with identity weighting.
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

estlatent_ave <- function(mod=NULL,Z,Y,data=NULL,method="sem",IV_Y=NULL,opt=TRUE){

  if(!is.null(data)){
    if(!is.data.frame(data)){stop("data must be a data.frame")}
    if(is.character(Z)){
      if(any(!(Z %in% colnames(data)))){stop("Some Z columns are not found in data")}
      Z <- data[,Z,drop=FALSE]
    }
    if(is.character(Y)){
      if(any(!Y %in% colnames(data))){stop("Some Y columns are not found in data")}
      Y <- data[,Y,drop=FALSE]
    }
  }

  if(is.data.frame(Z)==FALSE){
    if(is.matrix(Z)){
      Z <- as.data.frame(Z, check.names = FALSE)
    }else{
      Z <- data.frame(Z=Z)
    }
  }
  if(is.data.frame(Y)==FALSE){
    Y <- as.data.frame(Y, check.names = FALSE)
  }

  if(is.null(colnames(Z))){colnames(Z) <- paste0("Z",seq_len(ncol(Z)))}
  if(is.null(colnames(Y))){colnames(Y) <- paste0("Y",seq_len(ncol(Y)))}

  if(sum(is.na(Z))>0 ){stop("Z has NAs; please remove or fill in NAs first")}

  if(nrow(Z)!=nrow(Y)){stop("Y and Z have different numbers")}
  if(!is.logical(opt) || length(opt)!=1 || is.na(opt)){stop("opt must be TRUE or FALSE")}
  if(!(method %in% c("sem","gmm"))){stop("method must be one of 'sem' or 'gmm'")}

  n_y <- ncol(Y)

  if(n_y==1){stop("The function needs more than 1 outcome measures")}

  for (j in 1:n_y) {
    if(sum(is.na(Y[,j]))>0 ){stop("Y has NAs; please remove or fill in NAs first")}
  }


  ### sem prepare (flexible model specification)
  cov_dat <- if(!is.null(data)) data else Z

  if(is.null(mod)){
    rhs_text <- paste(colnames(Z),collapse = "+")
  }else{
    if(!is.character(mod) || length(mod)!=1){stop("mod must be a single character string")}
    mod_txt <- trimws(mod)
    if(startsWith(mod_txt,"~")){
      rhs_text <- trimws(sub("^~","",mod_txt))
    }else if(grepl("~",mod_txt,fixed=TRUE)){
      parts <- strsplit(mod_txt,"~",fixed=TRUE)[[1]]
      if(length(parts)!=2){stop("mod must be in the form 'lhs ~ rhs'")}
      lhs <- trimws(parts[1])
      rhs_text <- trimws(parts[2])
      if(nchar(lhs)>0 && !(lhs %in% c("Y","eta"))){stop("mod lhs must be 'Y' or 'eta' (or omit lhs and use '~rhs')")}
    }else{
      rhs_text <- mod_txt
    }
    if(nchar(rhs_text)==0){stop("mod rhs is empty")}
  }

  mm <- tryCatch(
    model.matrix(as.formula(paste0("~",rhs_text)), data = cov_dat),
    error = function(e) stop("Variables in mod must be available in `data` (or in `Z` when `data` is NULL).")
  )
  if("(Intercept)" %in% colnames(mm)){
    mm <- mm[,colnames(mm)!="(Intercept)",drop=FALSE]
  }
  if(ncol(mm)==0){stop("mod produced no predictors after removing intercept")}
  colnames(mm) <- make.names(colnames(mm), unique = TRUE)
  cov_mm <- data.frame(mm, check.names = FALSE)

  z_names <- make.names(colnames(Z), unique = FALSE)
  missing_z <- setdiff(z_names, colnames(cov_mm))
  if(length(missing_z)>0){
    stop("mod must include all treatment variables in Z")
  }
  x_names <- setdiff(colnames(cov_mm), z_names)
  n_x <- length(x_names)
  w_names <- colnames(cov_mm)

  ## Build IV sets for GMM:
  ## iv_load (user-specified) for loading moments;
  ## iv_reg automatically uses covariates from mod
  iv_pool <- data.frame(check.names = FALSE)
  if(!is.null(data)){
    iv_pool <- data
  }
  for(nm in colnames(cov_mm)){
    if(!(nm %in% colnames(iv_pool))){iv_pool[[nm]] <- cov_mm[[nm]]}
  }
  for(nm in colnames(Y)){
    if(!(nm %in% colnames(iv_pool))){iv_pool[[nm]] <- Y[[nm]]}
  }
  for(nm in colnames(Z)){
    if(!(nm %in% colnames(iv_pool))){iv_pool[[nm]] <- Z[[nm]]}
  }
  if(is.null(IV_Y)){
    iv_names_load <- z_names
  }else{
    if(!is.character(IV_Y)){stop("IV_Y must be a character vector of variable names")}
    iv_names_load <- unique(IV_Y)
  }
  if(any(!iv_names_load %in% colnames(iv_pool))){
    stop("Some IV names are not found in available variables (data, outcomes, treatment, and mod covariates)")
  }
  iv_load_mm <- iv_pool[,iv_names_load,drop=FALSE]
  iv_names_reg <- colnames(cov_mm)
  iv_reg_mm <- cov_mm[,iv_names_reg,drop=FALSE]
  y_names <- colnames(Y)
  n_iv_load <- ncol(iv_load_mm)
  n_iv_reg <- ncol(iv_reg_mm)
  n_w <- ncol(cov_mm)

  reg_text <- paste0("eta ~ ",paste(colnames(cov_mm),collapse = "+"))
  dat <- cbind(cov_mm,Y)

  ld_text <-  paste0("eta =~ 1*", paste0(colnames(Y),collapse = "+"))
  var_text <- paste0(z_names,"~~",z_names, collapse = "\n")

  mod_c <- paste(ld_text,reg_text,var_text,sep="\n")

  sem_tmp <- sem(mod_c,data=dat)

  ### sem results

  pe <- summary(sem_tmp)$pe

  ld <- pe[pe$lhs=="eta" & pe$op=="=~",]
  ld_est <- setNames(ld$est,ld$rhs)
  ld_se <- setNames(ld$se,ld$rhs)
  ld_z <- setNames(ld$z,ld$rhs)
  ld_p <- setNames(ld$pvalue,ld$rhs)
  sem_lambda_est <- as.numeric(ld_est[colnames(Y)])
  sem_lambda_se <- as.numeric(ld_se[colnames(Y)])
  sem_lambda_z <- as.numeric(ld_z[colnames(Y)])
  sem_lambda_p <- as.numeric(ld_p[colnames(Y)])

  rg <- pe[pe$lhs=="eta" & pe$op=="~",]
  rg_est <- setNames(rg$est,rg$rhs)
  rg_se <- setNames(rg$se,rg$rhs)
  rg_z <- setNames(rg$z,rg$rhs)
  rg_p <- setNames(rg$pvalue,rg$rhs)

  if(any(!z_names %in% names(rg_est))){stop("SEM result does not contain all treatment coefficients")}
  sem_beta_est <- as.numeric(rg_est[z_names])
  sem_beta_se <- as.numeric(rg_se[z_names])
  sem_beta_z <- as.numeric(rg_z[z_names])
  sem_beta_p <- as.numeric(rg_p[z_names])

  if(n_x>0){
    sem_x_est <- as.numeric(rg_est[x_names])
    sem_x_se <- as.numeric(rg_se[x_names])
    sem_x_z <- as.numeric(rg_z[x_names])
    sem_x_p <- as.numeric(rg_p[x_names])
  }

  sem_eta_var <- pe$est[pe$lhs=="eta" & pe$op=="~~" & pe$rhs=="eta"][1]
  yv <- pe[pe$op=="~~" & pe$lhs==pe$rhs & pe$lhs %in% colnames(Y),]
  sem_var_y <- as.numeric(setNames(yv$est,yv$lhs)[colnames(Y)])

  if(method=="sem"){
    final_lambda_est <- sem_lambda_est
    final_lambda_se <- sem_lambda_se
    final_lambda_p <- sem_lambda_p
    final_lambda_z <- sem_lambda_z
    final_beta_est <- sem_beta_est
    final_beta_se <- sem_beta_se
    final_beta_z <- sem_beta_z
    final_beta_p <- sem_beta_p

    if(n_x >0){
      final_x_est <- sem_x_est
      final_x_se <- sem_x_se
      final_x_z <- sem_x_z
      final_x_p <- sem_x_p
    }
  }

  ### gmm (average-proxy moments via gmm_ave)
  if(method=="gmm"){
    dat_gmm <- cbind(Y, cov_mm, iv_load_mm, iv_reg_mm)
    w_idx <- (n_y+1):(n_y+n_w)
    iv_load_idx <- (n_y+n_w+1):(n_y+n_w+n_iv_load)
    iv_reg_idx <- (n_y+n_w+n_iv_load+1):(n_y+n_w+n_iv_load+n_iv_reg)
    iv_names_norm <- make.names(iv_names_load)
    y_names_norm <- make.names(y_names)
    n_load_mom <- 0
    for(j in 2:n_y){
      bad_j <- which(iv_names_norm %in% c(y_names_norm[1], y_names_norm[j]))
      n_load_mom <- n_load_mom + (n_iv_load - length(bad_j))
    }
    # moments: loading + cov(Y1,Yj; j=2..n_y) + var(epsilon/eta block) + E[v]=0 + E[v*W]=0
    n_mom <- n_load_mom + (n_y - 1) + n_y + 1 + n_w
    # params: sigma_eta + lambda_2..lambda_ny + psi_1..psi_ny + alpha + beta vector
    n_par <- 2*n_y + 1 + n_w
    if(n_mom < n_par){
      stop("gmm_ave under-identified: add more IVs via IV=...")
    }
    g_opt <- function(theta, dat_gmm) {
      gmm_ave(theta, dat_gmm, mod = mod, n_y = n_y, w_idx = w_idx,
              iv_load_idx = iv_load_idx, iv_reg_idx = iv_reg_idx,
              iv_names_load = iv_names_load, iv_names_reg = iv_names_reg, y_names = y_names)
    }
    # starting value
    {
      sem_alpha <- pe$est[pe$lhs=="eta" & pe$op=="~1"][1]
      if(length(sem_alpha)==0 || is.na(sem_alpha)){sem_alpha <- 0}
      theta_s <- c(sem_eta_var, sem_lambda_est[2:n_y], sem_var_y, sem_alpha, as.numeric(rg_est[w_names]))
      res_opt <-momentModel(g_opt,dat_gmm,theta0=theta_s,vcov="iid")
      fit_type <- if(opt) "twostep" else "onestep"
      rec_gmm_opt <- gmmFit(res_opt, type = fit_type)
      final_lambda_est <- c(1,coef(rec_gmm_opt)[2:n_y])
      final_lambda_se <- c(0,sqrt(diag(vcov(rec_gmm_opt))[2:n_y]))
      final_lambda_z <- c(NA,final_lambda_est[-1]/final_lambda_se[-1])
      final_lambda_p <- 2 * (1 - pnorm(abs(final_lambda_z)))   # two‑sided
      b_idx <- (2*n_y+2):(2*n_y+1+n_w)
      beta_pos <- b_idx[match(z_names,w_names)]
      if(any(is.na(beta_pos))){stop("Internal error: treatment variables not found in GMM coefficient map")}
      final_beta_est <- coef(rec_gmm_opt)[beta_pos]
      final_beta_se <- sqrt(diag(vcov(rec_gmm_opt))[beta_pos])
      final_beta_z <- final_beta_est/final_beta_se
      final_beta_p <- 2 * (1 - pnorm(abs(final_beta_z)))   # two‑sided
      if(n_x>0){
        x_pos <- b_idx[match(x_names,w_names)]
        final_x_est <- coef(rec_gmm_opt)[x_pos]
        final_x_se <- sqrt(diag(vcov(rec_gmm_opt))[x_pos])
        final_x_z <- final_x_est/final_x_se
        final_x_p <- 2 * (1 - pnorm(abs(final_x_z)))
      }
    }

  }


  ### ICW


  ### output

  lamdba_name <- rep(NA,n_y)

  for (j in 1:n_y) {
    lamdba_name[j] <- paste0("lambda_",j)
  }
  trt_name <- z_names

  if(n_x>0){
    coef_name <- c(lamdba_name,trt_name,x_names)
    coef_est <- c(final_lambda_est,final_beta_est,final_x_est)
    coef_se <- c(final_lambda_se,final_beta_se,final_x_se)
    coef_z <- c(final_lambda_z,final_beta_z,final_x_z)
    coef_p <- c(final_lambda_p,final_beta_p,final_x_p)
  }else{
    coef_name <- c(lamdba_name,trt_name)
    coef_est <- c(final_lambda_est,final_beta_est)
    coef_se <- c(final_lambda_se,final_beta_se)
    coef_z <- c(final_lambda_z,final_beta_z)
    coef_p <- c(final_lambda_p,final_beta_p)
  }

  out_n <- n_y + length(trt_name) + n_x
  output <- matrix(NA, nrow = out_n, ncol = 4)
  output[,1] <- coef_est
  output[,2] <- coef_se
  output[,3] <- coef_z
  output[,4] <- coef_p

  rownames(output) <- coef_name
  colnames(output) <- c("Estimate", "SE","Z","P-value")

  fit <- list(
    call = match.call(),
    method = method,
    mod = mod,
    IV = if(exists("iv_names_load")) iv_names_load else NULL,
    coefficients = output
  )
  class(fit) <- "estlatent"
  fit

}
