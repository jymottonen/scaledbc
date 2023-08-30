#' Brute-force application of scaled Box-Cox transformation in flexmix
#'
#' \code{stepmixl2} is used to ....
#'
#' @param y response variable (in long format).
#' @param x predictor(s) as a vector or a matrix (in long format).
#' @param id identifier unique to each subject (in long format).
#' @param lambdas sequence of lambda calculation points.
#' @param K number of mixture components, single value or a vector.
#' @param classes known classes of the subjects for computing cluster purity (optional).
#' @param data an optional data frame containing the variables y, x and id.
#'
#' @details
#' Here are the details of the function...
#' @return A list containing the following components:
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
#' \item{models}{flexmix model objects fitted with best lambda value for each K}
#' }
#' @references
#' ...
#' @seealso
#' ...
#' @examples
#' \dontrun{
#' library(scaledbc)
#' summary(ex)
#' res <- stepmixl(y,x,id,lambdas=seq(-2,2,0.05),K=1:5,data=ex.long)
#' plot(res)
#' summary(res)
#' summary(res,digits=7)
#' }
#' @export
stepmixl2 <- function(y, x, id, lambdas, K, classes, data){
  require(flexmix)
  nlambda <- length(lambdas)
  if(hasArg(data)){y<-data$y; x<-data$x; id<-data$id}
  x<-as.data.frame(x)
  colnames(x)<-paste0("x",1:ncol(x))
  data<-cbind(y,x,id)
  (fmla <- as.formula(paste("tran(y, lambdas[l]) ~ ", paste(colnames(x), collapse= "+")," | id")))
  print(fmla)
  (fmlb <- as.formula(paste("tran(y, best_lambda) ~ ", paste(colnames(x), collapse= "+")," | id")))
  res <- matrix(nr=0, nc=10)
  models <- c()
  startTime <- Sys.time()
  for(k in K){
    likelihood <- c()
    for(l in 1:nlambda){
      likelihood <- c(likelihood, logLik(stepFlexmix(fmla, data=data, k=k, nrep=10, verbose=F)))
    }
    best_lambda <- lambdas[which.max(likelihood)]
    best <- stepFlexmix(fmlb, data=data, k=k, nrep=10, verbose=F)
    if(hasArg(classes)){
      pur <- purity(classes=classes, clusters=best@cluster[1:length(unique(id))])
    }
    else{
      pur <- NA
    }
    pass <- test.resid(best, tran(y, best_lambda))
    conv <- as.logical(best@converged)
    iter <- best@iter
    logL <- best@logLik

    res <- rbind(res, c(k, best_lambda, AIC(best), BIC(best), ICL(best), pur, pass, conv, iter, logL))
    models <- c(models, best)
    nowTime <- Sys.time()
    print(paste("K=",k))
    print(nowTime-startTime)
    startTime<-nowTime
  }
  res <- as.data.frame(res)
  names(res) <- c("K", "lambda", "AIC", "BIC", "ICL", "purity", "pass", "conv", "iter", "logL")
  fit<-list(K=res$K,lambda=res$lambda,AIC=res$AIC,BIC=res$BIC,ICL=res$ICL,purity=res$purity,
            pass=res$pass,conv=res$conv,iter=res$iter,logL=res$logL,models=models)
  class(fit) <- "stepmixl"
  return(fit)
}


