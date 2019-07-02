##   R functions used in:
##   Bremnes(2019). Constrained quantile regression splines for ensemble postprocessing.
##   Monthly Weather Review.
##
##


##  function for fitting quantiles
##
##  ens            matrix or data.frame of ensemble forecasts. Dimension: #forecasts, #ensemble members
##  obs            vector of observations
##  knots          vector of knots. If not specified equidistant knots with number of interior knots
##                 set equal to argument knots.interior
##  knots.interior number of interior knots. Not used if knots is specified
##  degree         degree of B-splines
##  reorder        if TRUE reorder B-spline coefficients in the probability direction
##  increasing     if TRUE reorder B-spline coefficients in ascending order in the covariate direction
##  nboot          number of bootstraps
##  aggregation    if TRUE B-spline coefficients are averaged over the bootstraps
##  c0             upper and lower limits for the spline fits
##  c1             upper and lower limits for the first order differences (after scaling covariates to [0,1])
##  c2             upper and lower limits for the second order differences (after scaling covariates to [0,1])
##
cqrs_ensemble.fit <- function(ens, obs, knots = NULL, knots.interior = 1, degree = 3, reorder = TRUE,
                              increasing = TRUE, nboot = 1, aggregation = TRUE,
                              c0 = c(0, NA), c1 = c(0, NA), c2 = c(NA, NA),
                              c2.first = c(NA, NA), c2.last = c(NA, NA)) {

    ##  define knots
    if (is.null(knots)) {
        knots <- seq(0, 1, length=knots.interior + 2)
        dx    <- knots[2] - knots[1]
        knots <- c(knots[1] - dx*degree:1, knots, knots[length(knots)] + dx*1:degree)
    }
    nbasis <- length(knots) - degree - 1  # number of B-splines
    
    ##  define quantile levels
    nqts  <- ncol(ens)
    prob  <- 1:nqts / (nqts + 1)
    qts   <- as.matrix(t(apply(ens, 1, sort)))  
     
    ##  organize data
    lim           <- apply(qts, 2, range)
    rownames(lim) <- c("lower", "upper")
    for (i in 1:nqts)
        qts[, i] <- (qts[, i] - lim[1, i]) / (lim[2, i] - lim[1, i])
    
    ##  define general constraints
    nbasis <- length(knots) - degree - 1
    R      <- r <- NULL
    constr <- rbind(c0, c1, c2)
    for (i in 1:nrow(constr)) {
        for (j in 1:ncol(constr)) {
            if (!is.na(constr[i,j])) {
                gt <- if (j == 1) 1 else -1
                R  <- rbind(R, if (i==1) gt*diag(nbasis) else gt*diff(diag(nbasis), differences=i-1))
                d  <- if (i == 1) 1 else dx*(i-1)
                r  <- c(r, gt * rep(constr[i,j], nbasis-i+1) * d)
            }
        }
    }
    
    ##  start and end point constraints (2nd order differences)
    R0 <- diff(diag(nbasis), differences=2)
    if (!is.na(c2.first[1])) {
        R <- rbind(R, R0[1,])
        r <- c(r, c2.start[1]*2*dx)
    }
    if (!is.na(c2.first[2])) {
        R <- rbind(R, -R0[1,])
        r <- c(r, -c2.start[2]*2*dx)
    }    
    if (!is.na(c2.last[1])) {
        R <- rbind(R, R0[nrow(R0),])
        r <- c(r, c2.end[1]*2*dx)
    }
    if (!is.na(c2.last[2])) {
        R <- rbind(R, -R0[nrow(R0),])
        r <- c(r, -c2.end[2]*2*dx)
    }  
    rm(R0)
    
    ##  fit splines
    coef <- array(NA, dim=c(length(prob), nbasis, nboot), dimnames=list(prob=prob, basis=NULL, boot=NULL))
    for (j in 1:nboot) {
        kb <- if (j == 1) 1:nrow(qts) else sample(nrow(qts), replace=TRUE)
        for (i in 1:nqts) {
            xt   <- splines::splineDesign(knots=knots, x=qts[kb, i], ord=degree+1, outer.ok=TRUE)
            if (is.null(R)) 
                try(coef[i,,j] <- quantreg::rq.fit.br(x=xt, y=obs[kb], tau=prob[i])$coefficients)
            else 
                try(coef[i,,j] <- quantreg::rq.fit.fnc(x=xt, y=obs[kb], tau=prob[i], R=R, r=r)$coefficients)
        }
    }
    
    ##  reorder spline coefficients in the probability direction
    if (reorder) {
        if (increasing & all(is.na(c1))) {
            coef <- apply(coef, c(1,3), sort, decreasing=FALSE)
            coef <- aperm(coef, c(2,1,3))
        }        
        coef <- apply(coef, 2:3, sort, decreasing=FALSE)
    }
        
    if (aggregation)
        coef <- apply(coef, 1:2, median, na.rm=TRUE)  ## some bootstrap samples may cause singularity/NA
    
    return( list(coefficients = coef, knots = knots, degree = degree, prob = prob, limits = lim,
                 obs.limits = range(obs), nboot = nboot, R = R, r = r) )
}



##  function for prediction
##
cqrs_ensemble.predict <- function(fit, ens, reorder=TRUE, extrapolation=FALSE) {
  nqts <- nrow(fit$coefficients)
  qts  <- as.matrix(t(apply(ens, 1, sort)))
  for (i in 1:nqts)
      qts[, i] <- (qts[, i] - fit$lim["lower", i]) / (fit$lim["upper", i] - fit$lim["lower", i])
  if (!extrapolation) {
      qts[qts>1] <- 1
      qts[qts<0] <- 0
  }
  out  <- matrix(NA, nrow=nrow(ens), ncol=ncol(ens), dimnames=list(NULL, fit$prob))
  for (i in 1:nqts) {
      xp       <- splines::splineDesign(knots=fit$knots, x=qts[, i], ord=fit$degree+1, outer.ok=TRUE)
      out[, i] <- xp %*% fit$coefficients[i, ]
  }

  if (reorder)
      out <- t(apply(out, 1, sort))
  
  return(out)
}
