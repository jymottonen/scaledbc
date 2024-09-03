#' summary.stepmixl
#'
#' summary.stepmixl is used to print the summary of
#' an object of class stepmixl.
#'
#' @param object an object of class stepmixl.
#' @param ... further arguments passed to or from other methods.
#' @param digits minimal number of significant digits.
#' @details
#' Here are the details of the function...
#' @method summary stepmixl
#' @export
summary.stepmixl<-function(object, ..., digits=6)
{
  res<-data.frame(K=object$K, conv=object$conv,
                  iter=object$iter, lambda=object$lambda, logL=object$logL,
                  AIC=object$AIC, BIC=object$BIC, ICL=object$ICL)
  print(res, digits=digits)
}


