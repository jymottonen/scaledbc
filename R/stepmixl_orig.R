#' Mixture model estimation without transformation in flexmix
#'
#' \code{stepmixl_orig} is used to ....
#'
#' @param y response variable (in long format)
#' @param x predictor(s) as a vector or a matrix (in long format)
#' @param id identifier unique to each subject (in long format)
#' @param K number of mixture components, single value or a vector
#' @param classes known classes of the subjects for computing cluster purity (optional)
#' @param data a data frame containing the variables y, x and id (optional)
#' @details
#' Fits the model onto the original, untransformed data.
#' @return A list with objects:
#' \describe{
#' \item{K}{number of mixture components}
#' \item{AIC}{Akaike information criterion (lower is better)}
#' \item{BIC}{Bayesian information criterion (lower is better)}
#' \item{ICL}{Integrated completed likelihood (lower is better)}
#' \item{purity}{cluster purity, higher is better (if classes are known)}
#' \item{pass}{proportion of clusters that pass Shapiro test for normality of residuals at alpha=0.05}
#' \item{conv}{1 if algorithm converged before reaching maximum number of iterations, 0 otherwise}
#' \item{iter}{number of iterations}
#' \item{logL}{log-likelihood of model}
#' \item{models}{flexmix model objects fitted with best lambda value for each K}
#' }
#' @references
#' Grün, B.  and Leisch, F. (2023), \emph{flexmix: Flexible Mixture Modeling}. R package version 2.3-19,
#' \url{https://CRAN.R-project.org/package=flexmix}.
#'
#' Leisch, F. (2004). FlexMix: A General Framework for Finite Mixture Models and Latent Class
#' Regression in R. \emph{Journal of Statistical Software}, \strong{11}(8), 1-18,
#' \url{https://doi.org/10.18637/jss.v011.i08}.
#'
#' Grün, B. and Leisch, F. (2007). Fitting Finite Mixtures of Generalized Linear Regressions in R,
#' \emph{Computational Statistics & Data Analysis}, \strong{51}, 5247–5252,
#' \url{https://doi.org/10.1016/j.csda.2006.08.014}.
#'
#' Grün, B. and Leisch, F. (2008). FlexMix Version 2: Finite Mixtures with Concomitant Variables
#' and Varying and Constant Parameters. \emph{Journal of Statistical Software}, \strong{28}(4), 1-35,
#' \url{https://doi.org/10.18637/jss.v028.i04}.
#'
#' @examples
#' \dontrun{
#' library(scaledbc)
#' res <- stepmixl_orig(y,x,id,K=1:5,data=ex.long)
#' plot(res)
#' summary(res)
#' summary(res,digits=7)
#' }
#' @export
#' @import flexmix
stepmixl_orig <- function(y,x,id,K, classes, data){
  #require(flexmix)
  if(hasArg(data)){y<-data$y; x<-data$x; id<-data$id}
  x<-as.data.frame(x)
  colnames(x)<-paste0("x",1:ncol(x))
  data<-cbind(y,x,id)
  (fmla <- as.formula(paste("y ~ ", paste(colnames(x), collapse= "+")," | id")))
  res <- matrix(nr=0, nc=9)
  models <- c()
  startTime <- Sys.time()
  for(k in K){
    best <- stepFlexmix(fmla, data=data, k=k, nrep=10, verbose=F)
    if(hasArg(classes)){pur <- purity(classes=classes, clusters=best@cluster[1:length(unique(id))])}
    else{pur <- NA}
    pass <- test.resid(best, y)
    conv <- as.logical(best@converged)
    iter <- best@iter
    logL <- best@logLik
    res <- rbind(res, c(k, AIC(best), BIC(best), ICL(best), pur, pass, conv, iter, logL))
    models <- c(models, best)
    nowTime <- Sys.time()
    print(paste("K=",k))
    print(nowTime-startTime)
    startTime<-nowTime
  }
  res <- as.data.frame(res)
  res$lambda<-NA
  names(res) <- c("K", "AIC", "BIC", "ICL", "purity", "pass", "conv", "iter", "logL","lambda")
  fit<-list(K=res$K,lambda=res$lambda,AIC=res$AIC,BIC=res$BIC,ICL=res$ICL,purity=res$purity,
            pass=res$pass,conv=res$conv,iter=res$iter,logL=res$logL,models=models)
  class(fit) <- "stepmixl"
  return(fit)
}

