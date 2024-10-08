#' plot.stepmixl
#'
#' plot.stepmixl is used to plot the results of
#' an object of class stepmixl.
#'
#' @param x an object of class stepmixl.
#' @param ... further arguments passed to or from other methods.
#' @details
#' Here are the details of the function...
#' @importFrom graphics plot points legend axis
#' @method plot stepmixl
#' @export
plot.stepmixl<-function(x, ...)
{
    uplim <- max(c(x$AIC, x$BIC, x$ICL))
    lowlim <- min(c(x$AIC, x$BIC, x$ICL))
    plot(x$AIC, xaxt="n", type="b", col="blue",
         xlab="Mixture components", ylab="Score",
         ylim=c(lowlim, uplim),
         pch=ifelse(x$AIC == min(x$AIC),16,1))
    points(x$BIC, type="b", col="red",
           pch=ifelse(x$BIC == min(x$BIC),16,1))
    points(x$ICL, type="b", col="green",
           pch=ifelse(x$ICL == min(x$ICL),16,1))
    legend(x="topright", legend=c("AIC","BIC","ICL"),
           col=c("blue","red","green"), pch=1)
    axis(1, at=x$K)
}
