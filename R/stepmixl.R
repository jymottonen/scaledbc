#' Application of Box-Cox transformation with flexmix
#'
#' \code{stepmixl} is used to ....
#'
#' @param y outcome vector
#' @param x predictor(s) as a vector or a matrix.
#' @param id identifier, unique to each subject.
#' @param K number of mixture components, single value or a vector.
#' @param seed RNG seed for reproducibility.
#' @param notrans if TRUE, analysis is done without transformation for comparison (TRUE by default)
#' @details
#' Here are the details of the function...
#' @return A list containing the following components:
#' \describe{
#' \item{models}{maximum likelihood solutions for different numbers of components. The numbers of components are given in K.}
#' \item{nt}{maximum likelihood solutions for different numbers of components without transformation, if notrans=TRUE}
#' \item{K}{number of mixture components}
#' \item{conv}{convergence of the stepFlexmix function. Equals 1 if converged.}
#' \item{iter}{number of iterations until convergence}
#' \item{lambda}{value of lambda}
#' \item{logL}{log-likelihood}
#' \item{AIC}{Akaike information criterion}
#' \item{BIC}{Bayesian information criterion}
#' \item{ICL}{Integrated completed likelihoo}
#' }
#' @references
#' ...
#' @seealso
#' ...
#' @examples
#' \dontrun{
#' library(scaledbc)
#' summary(ex)
#' res <- stepmixl(ex$y,ex$x,ex$id,K=1:5,seed=1)
#' plot(res)
#' summary(res)
#' summary(res,digits=7)
#' }
#' @export
stepmixl <- function(y, x, id, K, seed=.Random.seed, notrans=TRUE){
  #-----------------------------------------------
  # Scaled power transformation
  tran <- function(y, lambda){
    gm <- exp(mean(log(y)))
    if (lambda==0){ yt <- gm*log(y) }
    else{ yt <- ((y^(lambda))-1)/(lambda*(gm^(lambda-1))) }
    yt
  }
  #-----------------------------------------------
  require(flexmix)
  if(any(y <= 0, na.rm=TRUE)){
    stop("All values of the response variable y must be positive ",
          "in the version of the transformation applied here.")
  }
  res <- matrix(nr=0, nc=8)
  data <- list(y=y, x=x, id=id)
  models <- c()
  set.seed(seed)
  startTime <- Sys.time()
  for(k in K){
    # The method first checks 11 values of lambda
    # in range [-5,5] with interval 1.
    # The lambda with highest likelihood is selected
    # as a mid-point for the next iteration.
    # On the second iteration the range is of length 2
    # and the interval between values is 0.2.
    # The final iteration is done similarly with
    # a range of length 0.4 and interval 0.04.
    lambda <- seq(-5, 5, 1)
    likelihood <- c()
    for(l in 1:11){
      likelihood <- c(likelihood,
                        logLik(stepFlexmix(tran(y, lambda[l])~x|id,
                        k=k, nrep=10, verbose=F, data=data)))
    }
    mid <- lambda[which.max(likelihood)]

    lambda <- round(seq(mid-1, mid+1, 0.2),2)
    likelihood <- c()
    for(l in 1:11){
      likelihood <- c(likelihood,
                        logLik(stepFlexmix(tran(y, lambda[l])~x|id,
                        k=k, nrep=10, verbose=F, data=data)))
    }
    mid <- lambda[which.max(likelihood)]

    lambda <- round(seq(mid-0.2, mid+0.2, 0.04),3)
    likelihood <- c()
    for(l in 1:11){
      likelihood <- c(likelihood,
                        logLik(stepFlexmix(tran(y, lambda[l])~x|id,
                        k=k, nrep=10, verbose=F, data=data)))
    }
    final_lambda <- lambda[which.max(likelihood)]

    best <- stepFlexmix(tran(y, final_lambda)~x|id,
                        k=k, nrep=10, verbose=F, data=data)
    conv <- best@converged
    iter <- best@iter
    logL <- best@logLik

    res <- rbind(res, c(k, conv, iter, final_lambda, logL,
                        AIC(best), BIC(best), ICL(best)))
    models <- c(models, best)
    nowTime <- Sys.time()
    print(paste("K=",k))
    print(nowTime-startTime)
    startTime<-nowTime
  }
  res <- as.data.frame(res)
  names(res) <- c("K", "conv", "iter", "lambda", "logL",
                  "AIC", "BIC", "ICL")
  nt <- NA
  if(notrans){
    nt <- stepFlexmix(y~x|id, k=K, nrep=10, verbose=F, data=data)
  }

 fit<-list(models=models, nt=nt, K=res$K, conv=res$conv,
           iter=res$iter, lambda=res$lambda, logL=res$logL,
           AIC=res$AIC, BIC=res$BIC, ICL=res$ICL )
 class(fit) <- "stepmixl"
 return(fit)
}


