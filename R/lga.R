"lga" <- function(x, k, biter=NULL, niter=10, showall=FALSE, scale=TRUE,
                  nnode=NULL, silent=FALSE){
    ## Scale the dataset (if required), otherwise parse to a matrix
    x <- as.matrix(x)
    if (any(is.na(x))) stop("Missing data in x")

    if(scale)
        x <- scale(x, center=FALSE, scale=TRUE)
    if (!is.numeric(x)) stop("Data does not appear to be numeric.\n")
    n <- nrow(x); d <- ncol(x)

    ## Set the number of random starts (if not provided)
    if (is.null(biter))
        if(is.na(biter <- lga.NumberStarts(n, d, k, p=0.95)))
            stop("NumberStarts ill-specified. Rerun, setting biter \n")

    if (!silent)
        cat("LGA Algorithm \nk =", k, "\tBiter =", biter, "\tNiter =", niter, "\n")

    ## Choose the starting hyperplane coefficients for each biter
    hpcoef <- list()
    for (j in 1:biter){
        ## Choose starting clusters
        clindex <- matrix(sample(1:n, size=k*d, replace=FALSE), nrow=k)
        ## form initial hyperplanes - each row is a hyperplane
        hpcoef[[j]] <- matrix(NA, nrow=k, ncol=(d+1))
        for (i in 1:k)
            hpcoef[[j]][i,] <- lga.orthreg(x[clindex[i,],])
    }

    if(is.null(nnode)){
        ## not parallel
        outputsl <- lapply(hpcoef, lga.iterate, x, k, d, n, niter)
    }
    else {
        ## parallel
        if (!silent) cat("Setting up nodes \n")
        if (!(require(snow))) stop("Can't find required packages: snow")
        cl <- makeCluster(nnode)
        clusterEvalQ(cl, library(lga))
        outputsl <- clusterApplyLB(cl, hpcoef, lga.iterate, x, k, d, n, niter)
        stopCluster(cl)
        if (!silent) cat("Nodes closed \n")
    }

    outputs <- matrix(unlist(outputsl), ncol = biter, nrow = n + 2)

    ## Find the number of converged results
    nconverg <- sum(outputs[n+1,])
    if (nconverg == 0)
        warning("LGA failed to converge for any iteration\n")
    if (!showall){
        ## remove any columns with NAs
        outputs <- outputs[,complete.cases(t(outputs)), drop=FALSE]
        if (nconverg != 0)
            outputs <- outputs[,outputs[n+1,]==1, drop=FALSE]
        outputs <- outputs[, which.min(outputs[n+2,]), drop=FALSE]
        if(ncol(outputs) > 1)
            outputs <- lga.CheckUnique(outputs)
    }
    if (!silent) {
        cat("\nFinished.\n")
        if(showall) cat("\nReturning all outputs \n", "")
    }
    fout <- list(cluster=outputs[1:n,], ROSS=outputs[n+2,], converged=outputs[n+1,], biter=biter, niter=niter, nconverg=nconverg, scaled=scale, k=k, x=x)
    class(fout) <- "lga"
    return(fout)
}


"print.lga" <- function(x, ...) { ## S3method for printing LGA class
    cat("Output from LGA:\n\nDataset scaled =", x$scale, "\tk =",x$k,"\tniter =",x$niter,"\tbiter =",
        x$biter,"\nNumber of converged =",x$nconverg,"\nBest clustering has ROSS =",x$ROSS,"\n")
}

"plot.lga" <- function(x, ...) { ## S3method for plotting LGA class
    nclust <- x$k;   d <- ncol(x$x)
    hp <- matrix(NA, ncol=nclust, nrow=(d+1))

    ## Fit the hyperplanes
    for (i in 1:nclust)
        hp[,i] <- lga.orthreg(x$x[x$cluster == i,])

    ## 2-d case
    if (d == 2){
        plot(x$x, col=x$cluster, type="p", xlab="",ylab="", ...)
        for (i in 1:nclust)
            abline(a= hp[3,i]/hp[2,i], b= -hp[1,i]/hp[2,i], lty=2, col=i)
    }
    else{
        ## We do a rather complicated trick here using frames, as we need to know what pair of axes we are plotting
        ## in order to do the correct hyperplane intersections.
	lgaPlotCounter <- 0  ## this increments for each plot
	sysFrame <- sys.frame(sys.nframe()) ## this is the pointer to the current frame

	refMatrix<- matrix(NA, ncol=2, nrow=0) ## this matrix tells us which pair of axes
	for (j in 2:d) for (i in 1:(j-1)) refMatrix <- rbind(refMatrix, c(i,j))

        ## data.frame makes splom compatible with earlier versions of R
        splom(data.frame(x$x), lower.panel=function(x,y,...){}, upper.panel=function(x, y, hp, clusters, sysframe, refMatrix, ...)
          {
              int <- nrow(hp); nclust <- ncol(hp)
              assign("lgaPlotCounter", get("lgaPlotCounter", envir=sysframe) + 1, envir=sysframe)
              refV <- refMatrix[get("lgaPlotCounter", envir=sysframe),]
              panel.xyplot(x, y, col=clusters, ...)
              for (i in 1:nclust)
                  panel.abline(a= hp[int,i]/hp[ refV[2], i], b= -hp[ refV[1], i]/hp[ refV[2], i], lty=2, col=i)
          }, hp=hp, clusters=x$cluster, sysframe=sysFrame, refMatrix=refMatrix)
    }
}

"lga.calculateROSS" <- function(hpcoef, xsc, n, d, groups){
    ## This function calculates the total Residual Orthogonal Sum of Squares for a given grouping
    dist <- (sweep(xsc %*% t(hpcoef[,1:d, drop=FALSE]), 2, hpcoef[,d+1, drop=FALSE], FUN="-"))^2
    return(sum(dist[cbind(1:n, groups)]))
}

"lga.CheckUnique" <- function(x){
    "CheckUnique.Rand" <- function(z) {
        sum(z^2)-0.5*(sum(apply(z,1,sum)^2)+ sum(apply(z,2,sum)^2))
    }
    d <- dim(x)[2]
    index <- rep(TRUE,d)
    for (i in 1:(d-1)){
        for (j in (i+1):d){
            y <- table(x[,i], x[,j])
            z <- CheckUnique.Rand(y)
            if (z==0) index[j] <- FALSE
        }
    }
    if (sum(index) > 1) { ## In the incredibly unlikely situation....
        warning("Two unique solutions with identical ROSS.")
    }
    return(x[,index, drop=FALSE])
}

"lga.dodist" <- function(y, coeff, k, d, n){
    ## This function calculates the (orthogonal) Residuals for different hyerplanes,
    ## and returns the closest for each observation
    dist <- (y %*% t(coeff[,1:d])- matrix(coeff[,d+1], ncol=k, nrow=n, byrow=TRUE))^2
    return(max.col(-dist, ties.method="first"))
}

"lga.iterate" <-  function(hpcoef, xsc, k, d, n, niter){
    ## give the function the inital set of hyperplanes (in hpcoef)
    groups <- lga.dodist(xsc, hpcoef, k, d, n)
    iter <- 0
    if (any(table(groups) < d) | (length(unique(groups)) < k))
        iter <- (niter+10)
    oldgroups <- 1
    converged <- FALSE
    while (!converged & iter < niter) {
        iter <- iter+1
        oldgroups <- groups
        for (i in 1:k){
            hpcoef[i,] <- lga.orthreg(xsc[groups==i,])
        }
        groups <- lga.dodist(xsc, hpcoef, k, d, n)
        if (any(table(groups) < d) | (length(unique(groups)) < k))
            iter <- (niter+10) # if there aren't enough obs in a group
        if (all(oldgroups==groups)) converged <- TRUE
    }
    return(c(groups, converged, lga.calculateROSS(hpcoef, xsc, n, d, groups)))
}

"lga.NumberStarts" <- function(n, d, k, p=0.95){
    ## Calculate the number of starting positions (from article)
    n1 <- ceiling(n/k)
    return(ceiling(log(1-p)/log(1-choose(n1, d)^k/choose(n1*k, k*d))))
}

"lga.orthreg" <- function(x){
    ## Perform orthogonal regression.
    ## We use svd rather than eigen for numerical stability (esp when eigenvalues are small)
    y <- scale(x, scale=FALSE)
    emat <- svd(y)$v[,dim(y)[2]]
    return(c(emat, emat %*% attr(y, 'scaled:center')))
}

