#' summary.stepmixl
#'
#' summary.stepmixl is used to print the summary of
#' stepmixl.
#'
#' @param object an object of class stepmixl.
#' @param ... further arguments passed to or from other methods.
#' @param digits minimal number of significant digits.
#' @details
#' Here are the details of the function...
#' @export
summary.stepmixl<-function(object, ..., digits=6)
{
  res<-data.frame(K=object$K, conv=object$conv,
                  iter=object$iter, lambda=object$lambda, logL=object$logL,
                  AIC=object$AIC, BIC=object$BIC, ICL=object$ICL)
  print(res, digits=digits)
}

#' summary.stepmixl_orig
#'
#' summary.stepmixl_orig is used to print the summary of
#' stepmixl_orig.
#'
#' @param object an object of class stepmixl_orig.
#' @param ... further arguments passed to or from other methods.
#' @param digits minimal number of significant digits.
#' @details
#' Here are the details of the function...
#' @export
summary.stepmixl_orig<-function(object, ..., digits=6)
{
  res<-data.frame(K=object$K, conv=object$conv,
                  iter=object$iter, logL=object$logL,
                  AIC=object$AIC, BIC=object$BIC, ICL=object$ICL)
  print(res, digits=digits)
}


