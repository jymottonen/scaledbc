#' Simulated data in long format
#'
#' The simulated dataset  has n=100 longitudinal measurements of length p=10.
#'
#' @format A list with objects:
#' \describe{
#'   \item{y}{outcome vector (length n*p=1000)}
#'   \item{x}{predictor vector (length n*p=1000)}
#'   \item{id}{vector of the subject ids (length n*p=1000)}
#'   \item{class}{vector of the known classes of the subjects (length n=100)}
#' }
"ex.long"

#' Simulated data in wide format
#'
#' The simulated dataset  has n=100 longitudinal measurements of length p=10.
#'
#' @format A list with objects:
#' \describe{
#'   \item{y}{outcome matrix (dimension 100x10)}
#'   \item{x}{predictor matrix (dimension 100x10)}
#'   \item{id}{vector of the subject ids (length n=100)}
#'   \item{class}{vector of the known classes of the subjects (length n=100)}
#' }
"ex.wide"
