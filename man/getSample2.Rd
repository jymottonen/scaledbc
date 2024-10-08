% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{getSample2}
\alias{getSample2}
\title{Data generation (this function uses gamma distribution)}
\usage{
getSample2(n = 200, p = 6, a = 1, b = 2, pi = c(0.5, 0.5))
}
\arguments{
\item{n}{number of individuals}

\item{p}{number of measurement points}

\item{a}{intercept (for all classes)}

\item{b}{slope (for class 2)}

\item{pi}{prior probabilities}
}
\value{
A list containing the following components:
\describe{
\item{1}{data in long format (id, y, time, time^2)}
\item{2}{data in wide format (mostly for plotting purposes)}
\item{3}{classes}
}
}
\description{
\code{getSample2} is used to ....
}
\details{
Here are the details of the function...
}
