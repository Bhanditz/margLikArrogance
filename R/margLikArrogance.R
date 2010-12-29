#########################################################
# Build and query a sparse histogram

.GetID <- function(point, center, width) {
  # Return the ID string of the histogram given a point
  return(paste(floor((point - center) / width), collapse=" "))
}

.PopulateHist <- function(points, heights, width) {
  # Return a histogram with those points added
  # Each bucket is associated with a tuple (height, count)
  # where count is the number of points in that bucket
  max.point <- points[heights == max(heights), , drop=FALSE][1, ]
  center <- max.point - width / 2
  
  hist.env <- new.env(hash=TRUE, parent=emptyenv())
  for (i in seq(along=heights)) {
    # Update the histogram with a single point
    id <- .GetID(points[i, ], center, width)
    if (exists(id, envir=hist.env))
      cur.tuple <- get(id, envir=hist.env)
    else cur.tuple <- c(Inf, 0)

    assign(id, c(min(heights[i], cur.tuple[1]), cur.tuple[2] + 1),
           envir=hist.env)
  }
  return(list(env=hist.env, center=center, width=width, dim=length(center)))
}

.Profile <- function(hist) {
  # Return number of histogram basis points in each bucket (descending order)
  ids <- ls(hist$env)
  counts <- sapply(ids, function(id) get(id, envir=hist$env)[2])
  names(counts) <- ids
  return(sort(counts, decreasing=TRUE))
}

.Lookup <- function(point, hist) {
  # Return the height of the histogram at that point
  id <- .GetID(point, hist$center, hist$width)
  return(ifelse(exists(id, envir=hist$env), get(id, envir=hist$env)[1], -Inf))
}

.SetHistNorm <- function(hist) {
  # Return histogram with the log of integrated area underneath set in hist$norm
  ids <- ls(hist$env)
  heights <- sapply(ids, function(id) get(id, envir=hist$env)[1])
  # max.height is separated out below due computer floating point limitations
  max.height <- max(heights) 
  hist$norm <- (log(hist$width ^ hist$dim * sum(exp(heights - max.height)))
                + max.height)
  return(hist)
}

.Coverage <- function(points, hist) {
  # Return proportion number of points that fall under defined buckets in hist
  heights <- apply(points, 1, function(point) .Lookup(point, hist))
  return(length(heights[heights > -Inf]) / length(heights))
}

.Dist <- function(p1, p2) {
  # Distance using the Linf norm
  return(max(abs(p1 - p2)))
}

.GetWidth <- function(theta.hist, theta.width, ll.hist, opt.prob=0.5) {
  # Find the "optimum" width of the histogram
  # "Optimum" here just means that about prob of the samples from theta.width
  # will fall under the kernels.
  #
  # To get bounds on the density:
  # max density would give us a small volume around a low ll point
  # min density would give us a large volume around a high ll point
  low.point <- theta.hist[which(ll.hist == min(ll.hist))[1], , drop=FALSE]
  high.point <- theta.hist[which(ll.hist == max(ll.hist))[1], , drop=FALSE]
  
  minll.dists <- apply(theta.hist, 1, function(theta) .Dist(theta, low.point))
  maxll.dists <- apply(theta.hist, 1, function(theta) .Dist(theta, high.point))
  small.dist <- 0.5 * min(maxll.dists[maxll.dists > 0])
  big.dist <- 2 * max(minll.dists)
  
  F <- function(width) {
    hist <- .PopulateHist(theta.hist, ll.hist, width)
    return(.Coverage(theta.width, hist) - opt.prob)
  }
  solve.out <- uniroot(F, lower=small.dist, upper=big.dist,
                       tol=min(.05, 2/nrow(theta.width)))
  return(solve.out$root)
}

.GetLimits <- function(hist) {
  # Return a matrix of upper and lower bounds for each dimension of the hist
  split.ids <- strsplit(ls(hist$env), " ")
  coords <- matrix(NA, ncol=hist$dim, nrow=length(split.ids))
  for (i in seq(along=split.ids))
    coords[i, ] <- as.integer(split.ids[[i]])

  max.coords <- apply(coords, 2, max)
  min.coords <- apply(coords, 2, min)
  
  result <- rbind(apply(coords, 2, min) * hist$width,
                  (apply(coords, 2, max) + 1) * hist$width)
  rownames(result) <- c("min", "max")
  return(result)
}

#########################################################
# Use the histogram to compute marginal likelihood

.MargLL <- function(hist, theta.imp, ll.imp) {
  # Return marginal log likelihood given normalized histogram
  n <- length(ll.imp)
  samples <- rep(0, n)
  for (i in seq(length=n))
    samples[i] <- exp(.Lookup(theta.imp[i, ], hist) - ll.imp[i])
  return(list(samples=samples, mll=-log(mean(samples))))
}

.Assert <- function(boolean, err.msg=NULL) {
  # Abort and display traceback and err.msg if boolean is false
  if (isTRUE(boolean))
    return(invisible(NULL))
  call.list <- as.list(match.call())
  cat("----------------------------------------------------\n", file=stderr())
  cat("Error: Expression \"", as.character(as.expression(call.list[[2]])),
      "\"\n    evaluates to \"", boolean, "\" instead of TRUE.\n",
      sep="", file=stderr())
  cat("----------------------------------------------------\n", file=stderr())
  stop(err.msg)
}

.CheckArgs <- function(theta, ll, cov.prob, bounds) {
  # Do initial validation of arguments
  .Assert(is.matrix(theta), "theta should be a matrix of posterior samples")
  .Assert(length(dim(theta)) == 2,
         "theta needs to be two dimensional: each row is a sample, and each column is a dimension of the parameter space")
  .Assert(is.numeric(ll), "ll should be a vector of log likelihoods")
  sample.count <- length(ll)
  .Assert(sample.count == nrow(theta),
          "There needs to be one log likelihood for each row in theta")
  .Assert(sample.count > 300, "Give at least 300 data points")
  .Assert(0 < cov.prob && cov.prob < 1,
          "Coverage probability must be between 0 and 1")
  .Assert(is.matrix(bounds) && nrow(bounds) == 2,
          "bounds should be a matrix with two rows")
  .Assert(ncol(bounds) == ncol(theta),
          "bounds needs to have the same number of columns as theta")
}

.SplitComponents <- function(theta, ll) {
  # Split the log-likelihood and parameter values into histogram, width,
  # and importance buckets.  Return a list with these separated.
  tot.count <- length(ll)
  # this is not very scientific!
  hist.count <- min(0.2 * tot.count, round(sqrt(tot.count) * 2))
  width.count <- 40 # this many are used to determine width
  imp.count <- tot.count - hist.count - width.count
  return(list(theta.hist=theta[1:hist.count, , drop=FALSE],
              theta.bw=theta[(hist.count + 1):(hist.count + width.count), ,
                             drop=FALSE],
              theta.imp=theta[(tot.count - imp.count + 1):tot.count, ,
                              drop=FALSE],
              ll.hist=ll[1:hist.count],
              ll.imp=ll[(tot.count - imp.count + 1):tot.count]))
}

.CheckBounds <- function(hist, bounds) {
  # Make sure the histogram is within the hypercube specified by bounds
  # also sets the limits of the histogram
  limits <- hist$limits <- .GetLimits(hist)
  if (!all((bounds[1, ] == -Inf) | (bounds[1, ] < limits[1, ])))
    stop(paste("Bounds error: histogram lower limits (",
               paste(limits[1, ], collapse=" "),
               ") lower than specified bounds (",
               paste(bounds[1, ], collapse=" "), ").", sep=""))
  if (!all((bounds[2, ] == Inf) | (bounds[2, ] > limits[2, ])))
    stop(paste("Bounds error: histogram upper limits (",
               paste(limits[2, ], collapse=" "),
               ") greater than specified bounds (",
               paste(bounds[2, ], collapse=" "), ").", sep=""))
  return(hist)
}

MarginalLikelihood <- function(theta, ll, cov.prob=0.5,
                               bounds=matrix(c(-Inf, Inf),
                                             ncol=ncol(theta), nrow=2)) {
  # Return the marginal log likelihood
  .CheckArgs(theta, ll, cov.prob, bounds)
  comps <- .SplitComponents(theta, ll)
  width <- .GetWidth(comps$theta.hist, comps$theta.bw,
                     comps$ll.hist, cov.prob)
  hist <- .PopulateHist(comps$theta.hist, comps$ll.hist, width)
  hist <- .SetHistNorm(hist)
  hist <- .CheckBounds(hist, bounds)
  ml.out <- .MargLL(hist, comps$theta.imp, comps$ll.imp)
  sd.samples <- sd(ml.out$samples) / sqrt(length(ml.out$samples))
  conf.interval <- -log(qnorm(c(.975, .025)) * sd.samples + mean(ml.out$samples))
  return(list(mll=ml.out$mll + hist$norm,
              width=width, conf.interval=conf.interval + hist$norm, hist=hist))
}

