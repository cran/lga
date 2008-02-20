\name{lga}
\alias{lga}
\alias{lga.default}
\alias{rlga}
\alias{rlga.default}
\title{Perform LGA/RLGA}
\description{
  Linear Grouping Analysis
}

\usage{
\method{lga}{default}(x, k, biter = NULL, niter = 10, showall = FALSE, scale = TRUE,
    nnode=NULL, silent=FALSE, ...)
\method{rlga}{default}(x, k, alpha=0.5, biter = NULL, niter = 10, showall = FALSE, scale = TRUE,
    nnode=NULL, silent=FALSE, ...)
}

\arguments{
  \item{x}{a numeric matrix.}
  \item{k}{an integer for the number of clusters.}
  \item{alpha}{a numeric value between 0.5 and 1. For the robust
  estimate of LGA, specifying the percentage of points in the best subset.}
  \item{biter}{an integer for the number of different starting
    hyperplanes to try.}
  \item{niter}{an integer for the number of iterations to attempt for
    convergence.}
  \item{showall}{logical.  If TRUE then display all the outcomes, not just
    the best one.}
  \item{scale}{logical.  Allows you to scale the data, dividing each
    column by its standard deviation, before fitting.}
  \item{nnode}{an integer of many CPUS to use for parallel
    processing.  Defaults to NULL i.e. no parallel processing.}
  \item{silent}{logical.  If TRUE, produces no text output during
    processing.}
  \item{...}{For any other arguments passed from the generic function.}
}

\details{
  This code tries to find k clusters using the lga algorithm described
  in \cite{Van Aelst et al (2006)}.  For each attempt, it has up to
  \code{niter} steps to get to convergence, and it does this from
  \code{biter} different starting hyperplanes.  It then selects the
  clustering with the smallest Residual Orthoganal Sum of Squareds.
  
  If \code{biter} is left as NULL, then it is selected via the equation
  given in \cite{Van Aeslt et al (2006)}.

  The function \code{rlga} is the robust equivalent to LGA, and is
  introduced in \cite{Garcia-Escudero et al (2008)}.

  Both functions are parallel computing aware via the \code{nnode}
  argument, and works with the package \code{snow}.  In order to use
  parallel computing, one of MPI (e.g. lamboot) or PVM is necessary.
  For further details, see the documentation for \code{snow}.
  
  Associated with the lga and rlga functions are a print method and a
  plot method (see the examples).  In the plot method, the fitted
  hyperplanes are also shown as dashed-lines when there are only two
  dimensions.  }

\value{
  An object of class \sQuote{\dQuote{lga}} with components
  \item{cluster}{a vector containing the cluster memberships.}
  \item{ROSS}{the Residual Orthogonal Sum of Squares for the solution.}
  \item{converged}{a logical. True if at least one solution has converged.}
  \item{biter}{the biter setting used.}
  \item{niter}{the niter setting used.}
  \item{nconverg}{the number of converged solutions (out of biter starts).}
  \item{scaled}{logical. Is the data scaled?}
  \item{k}{the number of clusters to be found.}
  \item{x}{the (scaled if selected) dataset.}
}

\references{
    Van Aelst, S. and Wang, X. and Zamar, R. and Zhu, R. (2006)
  \sQuote{Linear Grouping Using Orthogonal Regression},
  \emph{Computational Statistics \& Data Analysis} \bold{50}, 1287--1312.

    Garcia-Escudero, L. and Gordaliz, A. and San Martin, R. and Van Aelst, S. and Zamar, R.H. (2008) \sQuote{Robust linear clustering}. Unpublished Manuscript.
}

\author{Justin Harrington \email{harringt@stat.ubc.ca}}

\seealso{\code{\link{gap}}}

\examples{

## Synthetic Data
## Make a dataset with 2 clusters in 2 dimensions

library(MASS)
set.seed(1234)
X <- rbind(mvrnorm(n=100, mu=c(1,-1), Sigma=diag(0.1,2)+0.9),
            mvrnorm(n=100, mu=c(1,1), Sigma=diag(0.1,2)+0.9))

lgaout <- lga(X,2)
plot(lgaout)
print(lgaout)

## Robust equivalent

rlgaout <- rlga(X,2, alpha=0.75)
plot(rlgaout)
print(rlgaout)


## nhl94 data set

data(nhl94)
plot(lga(nhl94, k=3, niter=30))


## Allometry data set
data(brain)
plot(lga(log(brain, base=10), k=3))


## Second Allometry data set
data(ob)
plot(lga(log(ob[,2:3]), k=3), pch=as.character(ob[,1]))


## Parallel processing case
## In this example, running using 4 nodes. 

\dontrun{
set.seed(1234)
X <- rbind(mvrnorm(n=1e6, mu=c(1,-1), Sigma=diag(0.1,2)+0.9),
            mvrnorm(n=1e6, mu=c(1,1), Sigma=diag(0.1,2)+0.9))
abc <- lga(X, k=2, nnode=4)
}

}

\keyword{multivariate}
\keyword{cluster}
