\name{lga-package}
\alias{lga-package}
\docType{package}

\title{The lga Package}
\description{
  The lga package is an implementation of the algorithms described in
  \cite{Van Aelst et al (2006)} and \cite{Garcia-Escudero et al (2008)}.  It has three main functions (with accompanying print and plot methods).
  \itemize{
    \item{\code{\link{lga}}}{The core linear grouping analysis}
    \item{\code{\link{rlga}}}{The robust version of linear grouping analysis}
    \item{\code{\link{gap}}}{Performs a gap analysis to find the number of clusters}
  }
  For more details refer to the relevant help files.
}

\references{
  Van Aelst, S. and Wang, X. and Zamar, R. and Zhu, R. (2006)
  \sQuote{Linear Grouping Using Orthogonal Regression},
  \emph{Computational Statistics \& Data Analysis} \bold{50},
  1287--1312.

  Garcia-Escudero, L. and Gordaliz, A. and San Martin, R. and Van Aelst, S. and Zamar, R.H. (2008) \sQuote{Robust linear clustering}. Unpublished Manuscript.
}

\author{Justin Harrington \email{harringt@stat.ubc.ca}}
\keyword{multivariate}
\keyword{cluster}