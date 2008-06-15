\name{nhl94}
\alias{nhl94}
\docType{data}
\title{Player Performance in NHL for 1994-1995}
\description{
  This data set gives four variables for each of 871 players for the NHL
  1994-1995 season.  From \cite{Van Aelst et al (2006)}
}
\usage{data(nhl94)}
\format{
  A matrix with 871 observations on the following 3 variables:
  \describe{
    \item{\code{PTS}}{Points scored, a numeric vector}
    \item{\code{PM}}{plus/minus average rating, a numeric vector}
    \item{\code{PIM}}{Total penalty time in minutes, a numeric vector}
    \item{\code{PP}}{Power play goals, a numeric vector}
  }
}

\references{
  Van Aelst, S. and Wang, X. and Zamar, R. and Zhu, R. (2006)
  \sQuote{Linear Grouping Using Orthogonal Regression},
  \emph{Computational Statistics \& Data Analysis} \bold{50}, 1287--1312.
}

\keyword{datasets}
  