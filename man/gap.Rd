\name{gap}
\alias{gap}
\title{Perform gap analysis}
\description{
  Performs the gap analysis using lga to estimate the number of clusters.
}

\usage{
gap(x, K, B, criteria=c("tibshirani", "DandF","none"),
    nnode=NULL, scale=TRUE)
}

\arguments{
  \item{x}{a numeric matrix.}
  \item{K}{an integer giving the maximum number of clusters to
    consider.}
  \item{B}{an integer giving the number of bootstraps.}
  \item{criteria}{a character string indicating which criteria to
    evaluate the gap data. One of \sQuote{\dQuote{tibshirani}}
    (default),\sQuote{\dQuote{DandF}} or \sQuote{\dQuote{none}}.  Can be
    abbreviated.}
  \item{nnode}{an integer of many CPUS to use for parallel
    processing.  Defaults to NULL i.e. no parallel processing.}
  \item{scale}{logical.  Should the data be scaled?}
}

\details{
  This code performs the gap analysis using lga.  The gap statistic is
  defined as the difference between the log of the Residual Orthogonal
  Sum of Squared Distances (denoted \eqn{log(W_k)}) and its expected
  value derived using bootstrapping under the null hypothesis that there
  is only one cluster.  In this implementation, the reference
  distribution used for the bootstrapping is a random uniform hypercube,
  transformed by the principal components of the underlying data set.
  For further details see \cite{Tibshirani et al (2001)}.

  For different criteria, different rules apply. With
  \sQuote{\dQuote{tibshirani}} (\cite{ibid}) we calculate the gap
  statistic for
  \eqn{k = 1, \ldots, K}, stopping when
  \deqn{\mbox{gap}(k) \ge \mbox{gap}(k+1) - s_{k+1}}{gap(k) >= gap(k+1)
    - s_(k+1)}
  where \eqn{s_{k+1}}{s_(k+1)} is a function of standard deviation of
  the bootstrapped estimates.
  
  With the \sQuote{\dQuote{DandF}} criteria from \cite{Dudoit et al
    (2002)}, we calculate the gap statistic for
  all values of \eqn{k = 1, \ldots, K}, selecting the number of clusters
  as
  \deqn{\hat{k} = \mbox{ smallest } k \ge 1 \mbox{such that gap}(k) \ge
    \mbox{gap}(k^*) - s_{k*}}{khat = smallest k >= 1 such that gap(k) >=
    gap(kstar) - s_(kstar)}
  where \eqn{k^* = \arg \max_{k \ge 1}\mbox{gap}(k)}{kstar = argmax_(k
    >= 1) gap(k)}.
  
  Finally, for the criteria \dQuote{none}, no rules are applied, and
  just the gap data is returned.
  
  As lga is ostensibly unsupervised in this case, the parameter niter
  is set to 20 to ensure convergence.

  This function is parallel computing aware via the \code{nnode}
  argument, and works with the package \code{snow}.  In order to
  use parallel computing, one of MPI (e.g. lamboot) or PVM is necessary.
  For further details, see the documentation for \code{snow}.

}

\value{
  An object of class \sQuote{\dQuote{gap}} with components
  \item{finished}{a logical.  For the \dQuote{tibshirani}, was there a
    solution found?}
  \item{nclust}{a integer for the number of clusters estimated.  Returns
    NA if nothing conclusive is found.}
  \item{data}{the original data set, scaled if specified in the
    arguments.}
  \item{criteria}{the criteria used.}
}

\references{
  Tibshirani, R. and Walther, G. and Hastie, T. (2001)
  \sQuote{Estimating the number of clusters in a data set via the gap
    statistic}, \emph{J. R. Statist. Soc. B} \bold{63}, 411--423.
  
  Dudoit, S. and Fridlyand, J. (2002) \sQuote{A prediction-based
    resampling method for estimating the number of clusters in a
    dataset}, \emph{Genome Biology} \bold{3}.
  
  Van Aelst, S. and Wang, X. and Zamar, R. and Zhu, R. (2006)
  \sQuote{Linear Grouping Using Orthogonal Regression},
  \emph{Computational Statistics \& Data Analysis} \bold{50},
  1287--1312.  }

\author{Justin Harrington \email{harringt@stat.ubc.ca}}

\seealso{\code{\link{lga}}}

\examples{

## Synthetic example
## Make a dataset with 2 clusters in 2 dimensions

library(MASS)
set.seed(1234)
X <- rbind(mvrnorm(n=100, mu=c(1, -2), Sigma=diag(0.1, 2) + 0.9),
           mvrnorm(n=100, mu=c(1, 1), Sigma=diag(0.1, 2) + 0.9))

gap(X, K=4, B=20)

## to run this using parallel processing with 4 nodes, the equivalent
## code would be

\dontrun{gap(X, K=4, B=20, nnode=4)}


## Quakes data (from package:datasets)
## Including the first two dimensions versus three dimensions
## yields different results

set.seed(1234)
\dontrun{
gap(quakes[,1:2], K=4, B=20)
gap(quakes[,1:3], K=4, B=20)
}

library(maps)
lgaout1 <- lga(quakes[,1:2], k=3)
plot(lgaout1)

lgaout2 <- lga(quakes[,1:3], k=2)
plot(lgaout2)

## Let's put this in context
par(mfrow=c(1,2))
map("world", xlim=range(quakes[,2]), ylim=range(quakes[,1])); box()
points(quakes[,2], quakes[,1], pch=lgaout1$cluster, col=lgaout1$cluster)

map("world", xlim=range(quakes[,2]), ylim=range(quakes[,1])); box()
points(quakes[,2], quakes[,1], pch=lgaout2$cluster, col=lgaout2$cluster)

}

\keyword{multivariate}

