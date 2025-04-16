#' Estimate the Causal Effects with Latent Outcomes
#'
#'@param Z a vector of Treatment variable.
#'@param Y a data frame that includes the measured outcomes. The outcome in the first column will be set as lambda=1.
#'@param X an optional data frame that includes other covariates. The defauly is NULL.
#'@param method estimators options.
#'@param eta tbd
#'
#'@import lavaan
estlatent <- function(Z,Y,X=NULL,eta = 1,method="sem"){

  n_z <- ncol(Z)
  n_y <- ncol(Y)

  ### sem

  if(is.null(X)){
    reg_text <- paste0("eta ~ ",colnames(Z))
    dat <- cbind(Z,Y)
  }else{
    reg_text <- paste0("eta ~ ",colnames(Z),"+", paste0(colnames(X),collapse = "+"))
    dat <- cbind(Z,Y,X)
  }

  ld_text <-  paste0("eta =~ 1*", paste0(colnames(Y),collapse = "+"))
  # reg_text <- paste0("eta ~ ",colnames(Z))
  var_text <- paste0(colnames(Z),"~~",  colnames(Z))

  mod_c <- paste(ld_text,reg_text,var_text,sep="\n")

  tmp <- sem(mod_c,data=dat)

  tmp_sum <- summary(tmp)

  return(tmp_sum)

  ### optimal


  ### ICW


  ##### method of moment


}
