#' Data generation (this function uses normal distribution)
#'
#' \code{getSample} is used to ....
#'
#' @param n number of individuals
#' @param p number of measurement points
#' @param a intercept (for all classes)
#' @param b slope (for classes 2 and 3)
#' @param c quadratic growth term (for class 3)
#' @param pi prior probabilities
#' @param v variance of random error
#' @details
#' Here are the details of the function...
#' @return A list containing the following components:
#' \describe{
#' \item{1}{data in long format (id, y, time, time^2)}
#' \item{2}{data in wide format (mostly for plotting purposes)}
#' \item{3}{classes}
#' }
#' @noRd
#' @export
getSample <- function(n=100, p=5, a=1, b=0.3, c=0.3, pi=c(0.65, 0.35, 0), v=1){
  obs <- function(subpop, p, a, b, c, v){
    x <- seq(0,1,1/(p-1))
    if (subpop == 1){
      ts(exp(a+rnorm(p, 0, v)), start = 0, end=1, frequency = p-1)
    }
    else if (subpop == 2){
      ts(exp(a+b*x+rnorm(p, 0, v)), start = 0, end=1, frequency = p-1)
    }
    else if (subpop == 3){
      ts(exp(a+b*x+c*(x^2)+rnorm(p, 0, v)), start = 0, end=1, frequency = p-1)
    }
  }
  data <- matrix(nr=n, nc=p, data=NA)
  subpops <- numeric(n)
  for(i in 1:n){
    subpop <- sample(1:3, 1, prob=pi)
    subpops[i] <- subpop
    data[i,] <- obs(subpop, p, a, b, c, v)
  }
  ts <- as.data.frame(data)
  y <- c(as.matrix(ts))
  id <- rep(1:n,p)
  time <- floor((0:(n*p-1))/n)/(p-1)
  time2 <- time^2
  data <- as.data.frame(cbind(id, y, time, time2))
  list(data, ts, subpops)
}

#' Data generation (this function uses gamma distribution)
#'
#' \code{getSample2} is used to ....
#'
#' @param n number of individuals
#' @param p number of measurement points
#' @param a intercept (for all classes)
#' @param b slope (for class 2)
#' @param pi prior probabilities
#' @details
#' Here are the details of the function...
#' @return A list containing the following components:
#' \describe{
#' \item{1}{data in long format (id, y, time, time^2)}
#' \item{2}{data in wide format (mostly for plotting purposes)}
#' \item{3}{classes}
#' }
#' @noRd
#' @export
getSample2 <- function(n=200, p=6, a=1, b=2, pi=c(0.5, 0.5)){
  obs2 <- function(subpop, a, b, p){
    x <- seq(0,1,1/(p-1))
    if (subpop == 1){
      ts(a+rgamma(p, 2, 1), start = 1, end=p, frequency = 1)
    }
    else if (subpop == 2){
      ts(a+b*x+rgamma(p, 2, 1), start = 1, end=p, frequency = 1)
    }
  }
  data <- matrix(nr=n, nc=p, data=NA)
  subpops <- sample(1:2, n, prob=pi, replace=T)
  for(i in 1:n){
    data[i,] <- obs2(subpops[i], a, b, p)
  }
  ts <- as.data.frame(data)
  y <- c(as.matrix(ts))
  id <- rep(1:n,p)
  time <- floor((0:(n*p-1))/n)+1
  data <- as.data.frame(cbind(id, y, time))
  list(data, ts, subpops)
}

#' Cluster purity
#'
#' \code{purity} is used to ....
#'
#' @param classes known classes
#' @param clusters clusters
#' @details
#' Here are the details of the function...
#' @noRd
#' @export
purity <- function(classes, clusters){
  tab <- table(clusters, classes)
  n <- length(classes)
  pure <- 0
  for(i in 1:length(tab[,1])){
    pure <- pure + max(tab[i,])
  }
  pure/n
}

#' Performs a Shapiro test of normality for each cluster in a model and
#' returns the proportion of clusters that "pass" the test at alpha = 0.05
#'
#' \code{test.resid} is used to ....
#'
#' @param model model
#' @param y y
#' @details
#' Here are the details of the function...
#' @noRd
#' @export
test.resid <- function(model, y){
  preds <- c()
  for(i in 1:length(y)){
    preds <- c(preds, fitted(model)[i,model@cluster[i]])
  }
  resids <- y-preds
  ps <- c()
  for(i in unique(model@cluster)){
    ps <- c(ps, shapiro.test(resids[model@cluster == i])$p)
  }
  sum(ps >= 0.05)/length(unique(model@cluster))
}

#' Scaled power transformation
#'
#' \code{tran} is used to ....
#'
#' @param y y
#' @param lambda lambda
#' @details
#' Here are the details of the function...
#' @noRd
#' @export
tran <- function(y, lambda){
  gm <- exp(mean(log(y)))
  if (lambda==0){ yt <- gm*log(y) }
  else{ yt <- ((y^(lambda))-1)/(lambda*(gm^(lambda-1))) }
  yt
}


