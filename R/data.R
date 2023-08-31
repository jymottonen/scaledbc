#' A simulated data in long format
#'
#' A simulated dataset with n=100 longitudinal measurements of length p=10.
#'
#' @format A list with objects:
#' \describe{
#'   \item{y}{outcome vector (length n*p=1000)}
#'   \item{x}{predictor vector (length n*p=1000)}
#'   \item{id}{vector of the subject ids (length n*p=1000)}
#'   \item{class}{vector of the known classes of the subjects (length n=100)}
#' }
"ex.long"

#' A simulated data in wide format
#'
#' A simulated dataset with n=100 longitudinal measurements of length p=10.
#'
#' @format A list with objects:
#' \describe{
#'   \item{y}{outcome matrix (dimension 100x10)}
#'   \item{x}{predictor matrix (dimension 100x10)}
#'   \item{id}{vector of the subject ids (length 100)}
#'   \item{class}{vector of the known classes of the subjects (length 100)}
#' }
"ex.wide"
