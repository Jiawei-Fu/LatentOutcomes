#' @export
print.estlatent <- function(x, ...) {
  cat("estlatent object. Use summary(x) to view the coefficient table.\n")
  invisible(x)
}

#' @export
summary.estlatent <- function(object, ...) {
  out <- list(
    call = object$call,
    method = object$method,
    mod = object$mod,
    IV = object$IV,
    coefficients = object$coefficients
  )
  class(out) <- "summary.estlatent"
  out
}

#' @export
print.summary.estlatent <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nMethod:", x$method, "\n")
  if(!is.null(x$mod)){
    cat("Model spec:", x$mod, "\n")
  }
  if(!is.null(x$IV)){
    cat("IV_Y:", paste(x$IV, collapse = ", "), "\n")
  }
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  invisible(x)
}
