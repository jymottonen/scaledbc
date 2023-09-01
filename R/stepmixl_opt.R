#' Optimized application of Box-Cox transformation in flexmix
#'
#' \code{stepmixl_opt} is used to ....
#'
#' @param y response variable (in long format)
#' @param x predictor(s) as a vector or a matrix (in long format)
#' @param id identifier unique to each subject (in long format)
#' @param K number of mixture components, single value or a vector
#' @param classes known classes of the subjects for computing cluster purity (optional)
#' @param data a data frame containing the variables y, x and id (optional)
#' @details
#' The method first checks 11 values of lambda
#' in range [-5,5] with interval 1.
#' The lambda with highest likelihood is selected
#' as a mid-point for the next iteration.
#' On the second iteration the range is of length 2
#' and the interval between values is 0.2.
#' The final iteration is done similarly with
#' a range of length 0.4 and interval 0.04.
#' @return A list with objects:
#' \describe{
#' \item{K}{number of mixture components}
#' \item{lambda}{best value of lambda found}
#' \item{AIC}{Akaike information criterion (lower is better)}
#' \item{BIC}{Bayesian information criterion (lower is better)}
#' \item{ICL}{Integrated completed likelihood (lower is better)}
#' \item{purity}{cluster purity, higher is better (if classes are known)}
#' \item{pass}{proportion of clusters that pass Shapiro test for normality of residuals at alpha=0.05}
#' \item{conv}{1 if algorithm converged before reaching maximum number of iterations, 0 otherwise}
#' \item{iter}{number of iterations}
#' \item{logL}{log-likelihood of model}
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
#' res <- stepmixl_opt(y,x,id,K=1:5,data=ex.long)
#' plot(res)
#' summary(res)
#' summary(res,digits=7)
#' }
#' @export
#' @import flexmix
stepmixl_opt <- function(y, x, id, K, classes, data){
  #require(flexmix)
  if(hasArg(data)){y<-data$y; x<-data$x; id<-data$id}
  x<-as.data.frame(x)
  colnames(x)<-paste0("x",1:ncol(x))
  data<-cbind(y,x,id)
  (fmla <- as.formula(paste("tran(y, lambda[l]) ~ ", paste(colnames(x), collapse= "+")," | id")))
  res <- matrix(nr=0, nc=10)
  models <- c()
  startTime <- Sys.time()
  for(k in K){
    lambda <- seq(-5, 5, 1)
    likelihood <- c()
    for(l in 1:11){
      likelihood <- c(likelihood,logLik(stepFlexmix(fmla,k=k, nrep=10, verbose=F, data=data)))
    }
    mid <- lambda[which.max(likelihood)]
    lambda <- round(seq(mid-1, mid+1, 0.2),2)
    likelihood <- c()
    for(l in 1:11){
      likelihood <- c(likelihood,logLik(stepFlexmix(fmla,k=k, nrep=10, verbose=F, data=data)))
    }
    mid <- lambda[which.max(likelihood)]
    lambda <- round(seq(mid-0.2, mid+0.2, 0.04),3)
    likelihood <- c()
    for(l in 1:11){
      likelihood <- c(likelihood,logLik(stepFlexmix(fmla,k=k, nrep=10, verbose=F, data=data)))
    }
    final_lambda <- lambda[which.max(likelihood)]
    (fmlb <- as.formula(paste("tran(y, final_lambda) ~ ", paste(colnames(x), collapse= "+")," | id")))
    best <- stepFlexmix(fmlb,k=k, nrep=10, verbose=F, data=data)
    if(hasArg(classes)){pur <- purity(classes=classes,clusters=best@cluster[1:length(classes)])}
    else{pur <- NA}
    pass <- test.resid(best, tran(data$y, final_lambda))
    conv <- best@converged
    iter <- best@iter
    logL <- best@logLik
    res <- rbind(res, c(k, final_lambda, AIC(best),BIC(best), ICL(best), pur, pass, conv, iter, logL))
    models <- c(models, best)
    nowTime <- Sys.time()
    print(paste("K=",k))
    print(nowTime-startTime)
    startTime<-nowTime
  }
  res <- as.data.frame(res)
  names(res) <- c("K","lambda","AIC","BIC","ICL","purity","pass","conv","iter","logL")
  fit<-list(K=res$K,lambda=res$lambda,AIC=res$AIC,BIC=res$BIC,ICL=res$ICL,purity=res$purity,
            pass=res$pass,conv=res$conv,iter=res$iter,logL=res$logL,models=models)
  class(fit) <- "stepmixl"
  return(fit)
}
