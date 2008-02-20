rlga <- function(x, k, alpha=0.5, ...)  UseMethod("rlga")

"rlga.default" <- function(x, k, alpha=0.5, biter=NULL, niter=10, showall=FALSE, scale=TRUE,
                  nnode=NULL, silent=FALSE, ...){
    ## Scale the dataset (if required), otherwise parse to a matrix
    x <- as.matrix(x)
    if (any(is.na(x))) stop("Missing data in x")

    if(scale)
        x <- scale(x, center=FALSE, scale=sqrt(apply(x,2,var)))
    if (!is.numeric(x)) stop("Data does not appear to be numeric.\n")
    n <- nrow(x); d <- ncol(x)

    ## Check alpha
    if (alpha < 0.5 || alpha > 1)
      stop("Alpha must be in [0.5, 1].")
    
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
        outputsl <- lapply(hpcoef, rlga.iterate, x, k, d, n, niter, alpha)
    }
    else {
        ## parallel
        if (!silent) cat("Setting up nodes \n")
        if (!(require(snow))) stop("Can't find required packages: snow")
        cl <- makeCluster(nnode)
        clusterEvalQ(cl, library(lga))
        outputsl <- clusterApplyLB(cl, hpcoef, rlga.iterate, x, k, d, n, niter, alpha)
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
    fout <- list(cluster=outputs[1:n,], ROSS=outputs[n+2,], converged=outputs[n+1,],
                 biter=biter, niter=niter, nconverg=nconverg, scaled=scale, k=k, x=x, alpha=alpha)
    class(fout) <- c("rlga","lga")
    return(fout)
}

"rlga.iterate" <-  function(hpcoef, xsc, k, d, n, niter, alpha){
    ## give the function the inital set of hyperplanes (in hpcoef)
    groups <- rlga.dodist(xsc, hpcoef, k, d, n, alpha)
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
        groups <- rlga.dodist(xsc, hpcoef, k, d, n, alpha)
        if (any(table(groups) < d) | (length(unique(groups)) < k))
            iter <- (niter+10) # if there aren't enough obs in a group
        if (all(oldgroups==groups)) converged <- TRUE
    }
    return(c(groups, converged, lga.calculateROSS(hpcoef, xsc, n, d, groups)))
}

"rlga.dodist" <- function(y, coeff, k, d, n, alpha){
    ## This function calculates the (orthogonal) Residuals for different hyerplanes,
    ## calculates which hyperplane each observation is closest to, and then takes the smallest alpha of them (setting the rest to zero)
    dist <- (y %*% t(coeff[,1:d])- matrix(coeff[,d+1], ncol=k, nrow=n, byrow=TRUE))^2
    mindist <- max.col(-dist, ties.method="first")
    mindist2 <- dist[cbind(1:n, mindist)]
    mindist[ mindist2 > quantile(mindist2, alpha) ] <- 0
    return(mindist)
}


"print.rlga" <- function(x, ...) { ## S3method for printing LGA class
    cat("Output from RLGA:\n\nDataset scaled =", x$scale, "\tk =",x$k,"\tniter =",x$niter,"\tbiter =", 
        x$biter,"\nAlpha =", x$alpha, "\nNumber of converged =",x$nconverg,"\nBest clustering has ROSS =",x$ROSS,"\n")
}

