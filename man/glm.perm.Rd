% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glm.perm.R
\name{glm.perm}
\alias{glm.perm}
\title{Fit GLM/Cox to a given model matrix x and response y}
\usage{
glm.perm(y, x, Family)
}
\arguments{
\item{y}{Response vector}

\item{x}{Model matrix with "resid" column for residuals}

\item{Family}{GLM family object}
}
\value{
a vector with (deviance, dispersion)
}
\description{
This is part of the implementation in the original glmperm package (GPL-3)
}
