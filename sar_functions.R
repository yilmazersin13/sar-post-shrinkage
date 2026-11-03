# Core functions for post-shrinkage estimation in the SAR model with specified W.

library(Matrix)

# Row-standardized W on a lattice (queen) or a 4-nearest-neighbour graph.

make_weights <- function(n, kind = "queen", seed = 27011) {
  if (kind == "queen") {
    layouts <- list("100" = c(10, 10), "250" = c(10, 25),
                    "500" = c(20, 25), "64" = c(8, 8))
    dims <- layouts[[as.character(n)]]
    if (is.null(dims)) stop("No lattice layout specified for this n.")
    xy <- as.matrix(expand.grid(x = seq_len(dims[1]), y = seq_len(dims[2])))
    dx <- abs(outer(xy[, 1], xy[, 1], "-"))
    dy <- abs(outer(xy[, 2], xy[, 2], "-"))
    adjacency <- 1 * (pmax(dx, dy) == 1)
  } else {
    set.seed(seed)
    xy <- cbind(runif(n), runif(n))
    distance <- as.matrix(dist(xy))
    diag(distance) <- Inf
    adjacency <- matrix(0, n, n)
    for (i in seq_len(n)) adjacency[i, order(distance[i, ])[1:4]] <- 1
    adjacency <- 1 * ((adjacency + t(adjacency)) > 0)
    if (kind == "distance") adjacency <- adjacency / distance
    if (!kind %in% c("knn", "distance")) stop("Unknown weights construction.")
  }
  diag(adjacency) <- 0
  degree <- rowSums(adjacency)
  stopifnot(all(degree > 0), max(abs(adjacency - t(adjacency))) < 1e-12)
  W <- adjacency / degree
  es <- eigen(adjacency / sqrt(outer(degree, degree)), symmetric = TRUE)
  vectors <- es$vectors
  left <- vectors / sqrt(degree)
  right <- t(vectors) * rep(sqrt(degree), each = n)
  # tr(G'G) = g' K g for g_i = eigenvalue_i / (1-rho*eigenvalue_i).
  K <- crossprod(vectors, degree * vectors) *
    crossprod(vectors, vectors / degree)
  list(W = W, values = es$values, left = left, right = right,
       K = K, degree = degree, xy = xy, n = n, kind = kind)
}

# S(rho)^-1 z
spatial_solve <- function(weights, z, rho) {
  as.vector(weights$left %*% ((weights$right %*% z) /
                               (1 - rho * weights$values)))
}

# log|S(rho)| = sum log(1 - rho * lambda_i)
log_jacobian <- function(rho, weights) {
  sum(log1p(-rho * weights$values))
}

make_design <- function(n, p, rho_x, distribution = "gaussian") {
  covariance <- rho_x ^ abs(outer(seq_len(p), seq_len(p), "-"))
  z <- if (distribution == "bounded") {
    matrix(runif(n * p, -sqrt(3), sqrt(3)), n, p)
  } else if (distribution == "gaussian") {
    matrix(rnorm(n * p), n, p)
  } else stop("Unknown design distribution.")
  X <- z %*% chol(covariance)
  X <- sweep(X, 2, colMeans(X), "-")

  X <- sweep(X, 2, sqrt(colMeans(X^2)), "/")
  colnames(X) <- paste0("x", seq_len(p))
  X
}

.
sar_data <- function(y, X, weights) {
  wy <- as.vector(weights$W %*% y)
  list(y = y, X = X, wy = wy, n = length(y), p = ncol(X),
       Q = crossprod(X), xy = as.vector(crossprod(X, y)),
       xwy = as.vector(crossprod(X, wy)),
       yy = sum(y^2), ywy = sum(y * wy), wywy = sum(wy^2),
       weights = weights)
}


sar_rss <- function(beta, rho, data) {
  value <- data$yy - 2 * rho * data$ywy + rho^2 * data$wywy
  if (length(beta)) value <- value - 2 * sum(beta * (data$xy - rho * data$xwy)) +
    sum(beta * (data$Q %*% beta))
  if (!is.finite(value) || value <= 0) stop("Nonpositive SAR residual sum of squares.")
  as.numeric(value)
}

scalar_minimum <- function(fn, interval = c(-0.98, 0.98), grid_size = 31) {
  grid <- seq(interval[1], interval[2], length.out = grid_size)
  val <- vapply(grid, fn, numeric(1))
  candidates <- which(val[2:(grid_size - 1)] <= val[1:(grid_size - 2)] &
                        val[2:(grid_size - 1)] <= val[3:grid_size]) + 1
  points <- c(grid[1], grid[grid_size])
  for (j in candidates) {
    fit <- optimize(fn, grid[c(j - 1, j + 1)], tol = 1e-9)
    points <- c(points, fit$minimum)
  }
  points[which.min(vapply(points, fn, numeric(1)))]
}

# SAR-QMLE (Lee, 2004). 

sar_qmle <- function(y, X, weights, interval = c(-0.98, 0.98)) {
  data <- sar_data(y, X, weights)
  if (data$p && qr(X)$rank < data$p) stop("Rank-deficient regression design.")
  by <- if (data$p) qr.coef(qr(X), y) else numeric(0)
  bw <- if (data$p) qr.coef(qr(X), data$wy) else numeric(0)
  ey <- y - as.vector(X %*% by)
  ew <- data$wy - as.vector(X %*% bw)
  aa <- sum(ey^2); bb <- sum(ey * ew); cc <- sum(ew^2)
  profile <- function(rho) {
    rss <- aa - 2 * rho * bb + rho^2 * cc
    if (rss <= 0) return(Inf)
    0.5 * data$n * log(rss / data$n) - log_jacobian(rho, weights)
  }
  rho <- scalar_minimum(profile, interval)
  beta <- as.vector(by - rho * bw)
  sigma2 <- sum((y - rho * data$wy - as.vector(X %*% beta))^2) / data$n
  loglik <- log_jacobian(rho, weights) -
    0.5 * data$n * (log(2 * pi * sigma2) + 1)
  list(beta = beta, rho = rho, sigma2 = sigma2, loglik = loglik,
       converged = TRUE, boundary = min(abs(rho - interval)) < 1e-5)
}

sar_information <- function(X, beta, rho, sigma2, weights) {
  n <- nrow(X); p <- ncol(X)
  g <- weights$values / (1 - rho * weights$values)
  gx <- as.vector(weights$left %*% (g * (weights$right %*% (X %*% beta))))
  bb <- crossprod(X) / sigma2
  br <- as.vector(crossprod(X, gx)) / sigma2
  rr <- sum(gx^2) / sigma2 + sum(g^2) + as.numeric(crossprod(g, weights$K %*% g))
  ee <- matrix(c(rr, sum(g) / sigma2, sum(g) / sigma2,
                 n / (2 * sigma2^2)), 2, 2)
  be <- cbind(br, rep(0, p))
  C <- (bb - be %*% solve(ee, t(be))) / n
  C <- (C + t(C)) / 2
  if (min(eigen(C, symmetric = TRUE, only.values = TRUE)$values) <= 0)
    stop("Nonpositive effective information.")
  list(C = C, information = rbind(cbind(bb, be), cbind(t(be), ee)) / n)
}

soft_threshold <- function(x, threshold) sign(x) * pmax(abs(x) - threshold, 0)

# Coordinate descent for the weighted Lasso step.
weighted_lasso <- function(Q, z, penalty, beta, tolerance = 1e-9, maxit = 3000) {
  for (iteration in seq_len(maxit)) {
    old <- beta
    for (j in seq_along(beta)) {
      partial <- z[j] - sum(Q[j, ] * beta) + Q[j, j] * beta[j]
      beta[j] <- soft_threshold(partial, penalty[j]) / Q[j, j]
    }
    if (max(abs(beta - old)) < tolerance * (1 + max(abs(beta)))) break
  }
  list(beta = beta, converged = iteration < maxit, iterations = iteration)
}

penalty_kkt <- function(beta, rho, sigma2, data, lambda, penalty, adaptive_weights) {
  gradient <- as.vector(data$Q %*% beta - data$xy + rho * data$xwy) /
    (data$n * sigma2)
  if (penalty == "AL") {
    active <- abs(beta) > 1e-9
    violation <- pmax(abs(gradient) - lambda * adaptive_weights, 0)
    violation[active] <- abs(gradient[active] +
                              lambda * adaptive_weights[active] * sign(beta[active]))
  } else violation <- abs(gradient + lambda * beta)
  grad_rho <- sum(data$weights$values / (1 - rho * data$weights$values)) / data$n -
    (data$ywy - rho * data$wywy - sum(beta * data$xwy)) / (data$n * sigma2)
  if (rho <= -0.98 + 1e-6) grad_rho <- min(grad_rho, 0)
  if (rho >= 0.98 - 1e-6) grad_rho <- max(grad_rho, 0)
  max(c(violation, abs(grad_rho)))
}


sar_penalized <- function(data, lambda, penalty, adaptive_weights, start,
                          tolerance = 2e-6, maxit = 300) {
  beta <- start$beta; rho <- start$rho; sigma2 <- start$sigma2
  old_objective <- Inf
  for (iteration in seq_len(maxit)) {
    z <- data$xy - rho * data$xwy
    if (penalty == "AL") {
      update <- weighted_lasso(data$Q, z, data$n * sigma2 * lambda * adaptive_weights, beta)
      beta <- update$beta
    } else {
      beta <- as.vector(solve(data$Q + diag(data$n * sigma2 * lambda, data$p), z))
    }
    conditional <- function(r) {
      -log_jacobian(r, data$weights) / data$n +
        sar_rss(beta, r, data) / (2 * data$n * sigma2)
    }
    rho <- optimize(conditional, c(-0.98, 0.98), tol = 1e-9)$minimum
    sigma2 <- sar_rss(beta, rho, data) / data$n
    pen <- if (penalty == "AL") sum(adaptive_weights * abs(beta)) else sum(beta^2) / 2
    objective <- 0.5 * log(sigma2) - log_jacobian(rho, data$weights) / data$n + lambda * pen
    kkt <- penalty_kkt(beta, rho, sigma2, data, lambda, penalty, adaptive_weights)
    if (kkt < tolerance && abs(old_objective - objective) < tolerance) break
    old_objective <- objective
  }
  loglik <- log_jacobian(rho, data$weights) -
    0.5 * data$n * (log(2 * pi * sigma2) + 1)
  list(beta = beta, rho = rho, sigma2 = sigma2, loglik = loglik,
       objective = objective, lambda = lambda, kkt = kkt,
       converged = kkt < tolerance, iterations = iteration,
       boundary = min(abs(rho - c(-0.98, 0.98))) < 1e-5)
}


tune_penalty <- function(y, X, weights, full, penalty, nlambda = 20) {
  data <- sar_data(y, X, weights)
  adaptive_weights <- 1 / (abs(full$beta) + 1 / sqrt(data$n))
  null <- sar_qmle(y, X[, FALSE, drop = FALSE], weights)
  null$beta <- rep(0, data$p)
  if (penalty == "AL") {
    lambda_max <- max(abs(data$xy - null$rho * data$xwy) /
                        (data$n * null$sigma2 * adaptive_weights))
    lambdas <- lambda_max * 10^seq(0, -4, length.out = nlambda)
  } else {
    scale <- median(diag(data$Q)) / (data$n * full$sigma2)
    lambdas <- scale * 10^seq(2, -5, length.out = nlambda)
  }
  fits <- vector("list", length(lambdas) + 1)
  path <- vector("list", length(fits))
  previous <- null
  for (j in seq_along(lambdas)) {
    a <- sar_penalized(data, lambdas[j], penalty, adaptive_weights, previous)
    b <- sar_penalized(data, lambdas[j], penalty, adaptive_weights, full)
    eligible <- Filter(function(z) z$converged, list(a, b))
    if (!length(eligible)) {
      fit <- if (a$objective < b$objective) a else b
    } else fit <- eligible[[which.min(vapply(eligible, function(z) z$objective, numeric(1)))]]
    fit$start_gap <- abs(a$objective - b$objective)
    df <- if (penalty == "AL") sum(abs(fit$beta) > 1e-8) else {
      sum(diag(solve(data$Q + diag(data$n * fit$sigma2 * lambdas[j], data$p), data$Q)))
    }
    fit$bic <- -2 * fit$loglik + (df + 2) * log(data$n)
    fit$df <- df
    fits[[j]] <- fit
    path[[j]] <- data.frame(lambda = lambdas[j], bic = fit$bic, df = df,
                           kkt = fit$kkt, converged = fit$converged,
                           start_gap = fit$start_gap)
    if (fit$converged) previous <- fit
  }
  full$lambda <- 0; full$kkt <- 0; full$df <- data$p
  full$bic <- -2 * full$loglik + (data$p + 2) * log(data$n)
  full$start_gap <- 0
  fits[[length(fits)]] <- full
  path[[length(path)]] <- data.frame(lambda = 0, bic = full$bic, df = data$p,
                                    kkt = 0, converged = TRUE, start_gap = 0)
  path <- do.call(rbind, path)
  # Do not silently select a path containing failed fits.
  if (any(!path$converged)) stop(paste(penalty, "path has nonconverged candidates."))
  chosen <- which.min(path$bic)
  fit <- fits[[chosen]]
  fit$selected <- which(abs(fit$beta) > 1e-8)
  fit$path <- path
  fit$edge <- chosen %in% c(1, length(lambdas))
  fit
}

# Post-shrinkage estimators of beta1 (paper Section 3):

post_shrinkage <- function(full, restricted, retained, C, n, alpha = 0.05) {
  p <- length(full$beta); excluded <- setdiff(seq_len(p), retained); q <- length(excluded)
  uf <- rep(0, p); uf[retained] <- restricted$beta
  estimates <- list(OF = full$beta, UF = uf, PT = full$beta, S = full$beta, PS = full$beta)
  if (!q) return(list(estimates = estimates, T = 0, q = q, fallback = FALSE))
  C22 <- C[excluded, excluded, drop = FALSE]
  if (length(retained)) C22 <- C22 - C[excluded, retained, drop = FALSE] %*%
    solve(C[retained, retained, drop = FALSE], C[retained, excluded, drop = FALSE])
  T <- as.numeric(n * crossprod(full$beta[excluded], C22 %*% full$beta[excluded]))
  if (!is.finite(T) || T < 0) stop("Invalid Wald statistic.")
  difference <- uf[retained] - full$beta[retained]
  estimates$PT[retained] <- full$beta[retained] +
    difference * (T <= qchisq(1 - alpha, q))
  if (q >= 3) {
    if (T == 0) stop("Stein rule undefined at exactly zero T.")
    estimates$S[retained] <- full$beta[retained] + difference * (q - 2) / T
    estimates$PS[retained] <- full$beta[retained] + difference * min(1, (q - 2) / T)
  }

  list(estimates = estimates, T = T, q = q, fallback = q > 0 && q < 3)
}

# Average direct, total and indirect impact factors (LeSage and Pace papers).
impact_factors <- function(rho, weights) {
  rho <- unname(rho)
  direct <- mean(1 / (1 - rho * weights$values))
  total <- 1 / (1 - rho)
  c(ADE = direct, ATE = total, AIE = total - direct)
}

inverse_chisq_moment <- function(df, ncp, power = 1, cutoff = NULL) {
  stopifnot(df > 2 * power, ncp >= 0)
  k <- 0:max(40, ceiling(qpois(1 - 1e-13, ncp / 2)))
  degrees <- df + 2 * k
  values <- exp(-power * log(2) + lgamma(degrees / 2 - power) - lgamma(degrees / 2))
  if (!is.null(cutoff)) values <- values *
    pchisq(cutoff, df = degrees - 2 * power, lower.tail = FALSE)
  sum(dpois(k, ncp / 2) * values)
}

# Asymptotic distributional risk of OF, UF, PT, S and PS for beta1 (Theorem 4),
adr_formula <- function(C, kappa, p1, alpha = 0.05) {
  p <- nrow(C); q <- p - p1; i <- seq_len(p1); j <- (p1 + 1):p
  A <- solve(C[i, i, drop = FALSE], C[i, j, drop = FALSE])
  C22 <- C[j, j, drop = FALSE] - C[j, i, drop = FALSE] %*% A
  Phi <- A %*% solve(C22, t(A)); mu <- as.vector(A %*% kappa)
  Rof <- sum(diag(solve(C)[i, i, drop = FALSE]))
  Rphi <- sum(diag(Phi)); Rmu <- sum(mu^2)
  Delta <- as.numeric(crossprod(kappa, C22 %*% kappa))
  a <- q - 2; critical <- qchisq(1 - alpha, q)
  h2 <- pchisq(critical, q + 2, ncp = Delta)
  h4 <- pchisq(critical, q + 4, ncp = Delta)
  m2 <- inverse_chisq_moment(q + 2, Delta, 1)
  m4 <- inverse_chisq_moment(q + 4, Delta, 1)
  v2 <- inverse_chisq_moment(q + 2, Delta, 2)
  v4 <- inverse_chisq_moment(q + 4, Delta, 2)
  g2 <- pchisq(a, q + 2, ncp = Delta) + a * inverse_chisq_moment(q + 2, Delta, 1, a)
  g4 <- pchisq(a, q + 4, ncp = Delta) + a * inverse_chisq_moment(q + 4, Delta, 1, a)
  gg2 <- pchisq(a, q + 2, ncp = Delta) + a^2 * inverse_chisq_moment(q + 2, Delta, 2, a)
  gg4 <- pchisq(a, q + 4, ncp = Delta) + a^2 * inverse_chisq_moment(q + 4, Delta, 2, a)
  risk <- c(OF = Rof, UF = Rof - Rphi + Rmu,
            PT = Rof - Rphi * h2 + Rmu * (2 * h2 - h4),
            S = Rof - 2 * a * Rphi * m2 - 2 * a * Rmu * (m4 - m2) +
              a^2 * Rphi * v2 + a^2 * Rmu * v4,
            PS = Rof - 2 * Rphi * g2 - 2 * Rmu * (g4 - g2) + Rphi * gg2 + Rmu * gg4)
  list(risk = risk, bias = list(OF = mu * 0, UF = mu, PT = mu * h2,
                                S = mu * a * m2, PS = mu * g2),
       Delta = Delta, A = A, C22 = C22, Phi = Phi, mu = mu)
}

# Missing-response functions, used for geographic validation in the applications.
# Held-out responses are integrated out of the Gaussian likelihood (Suesse, 2018).

marginal_data <- function(y_train,X,w,train) {
  n <- nrow(X); held <- setdiff(seq_len(n),train)
  stopifnot(length(y_train)==length(train),!anyDuplicated(train),
    all(is.finite(y_train)),qr(X[train,,drop=FALSE])$rank==ncol(X))
  W <- Matrix(w$W,sparse=TRUE)
  y0 <- numeric(n); y0[train] <- y_train # zeros are algebraic placeholders, not imputation
  wy0 <- as.vector(W%*%y0)
  xx <- crossprod(X); xy <- as.vector(crossprod(X,y0)); xwy <- as.vector(crossprod(X,wy0))
  yy <- sum(y0^2); ywy <- sum(y0*wy0); wywy <- sum(wy0^2)
  if(length(held)) {
    B1 <- as.matrix(W[held,held,drop=FALSE]+t(W[held,held,drop=FALSE]))
    B2 <- as.matrix(crossprod(W[,held,drop=FALSE]))
    mh0 <- X[held,,drop=FALSE]; mh1 <- as.matrix(crossprod(W[,held,drop=FALSE],X))
    v1 <- as.vector((W%*%y0+t(W)%*%y0)[held])
    v2 <- as.vector(crossprod(W,wy0)[held])
  }
  stats <- function(rho,derivatives=FALSE) {
    R <- xx; b <- xy-rho*xwy; a <- yy-2*rho*ywy+rho^2*wywy
    ld <- log_jacobian(rho,w)
    if(length(held)) {
      Qh <- diag(length(held))-rho*B1+rho^2*B2
      ch <- chol(Qh)
      F <- chol2inv(ch)
      M <- mh0-rho*mh1; v <- -rho*v1+rho^2*v2
      FM <- F%*%M; Fv <- as.vector(F%*%v)
      R <- R-crossprod(M,FM); b <- b-as.vector(crossprod(M,Fv))
      a <- a-sum(v*Fv); ld <- ld-sum(log(diag(ch)))
    }
    ans <- list(R=(R+t(R))/2,b=b,a=a,ld=ld)
    if(derivatives) {
      dR <- matrix(0,ncol(X),ncol(X)); db <- -xwy; da <- -2*ywy+2*rho*wywy
      dld <- -sum(w$values/(1-rho*w$values))
      if(length(held)) {
        dQ <- -B1+2*rho*B2; dM <- -mh1; dv <- -v1+2*rho*v2
        cross <- crossprod(dM,FM)
        dR <- -cross-t(cross)+crossprod(FM,dQ%*%FM)
        db <- db-as.vector(crossprod(dM,Fv)+crossprod(FM,dv)-crossprod(FM,dQ%*%Fv))
        da <- da-2*sum(dv*Fv)+as.numeric(crossprod(Fv,dQ%*%Fv))
        dld <- dld-.5*sum(F*t(dQ))
      }
      ans$dR <- dR; ans$db <- db; ans$da <- da; ans$dld <- dld
    }
    ans
  }
  list(stats=stats,n=length(train),p=ncol(X),X=X,w=w,train=train,held=held,y_train=y_train)
}
marginal_rss <- function(g,beta) {
  rss <- as.numeric(g$a-2*sum(beta*g$b)+crossprod(beta,g$R%*%beta))
  if(!is.finite(rss)||rss<=0) stop("Nonpositive marginal RSS.")
  rss
}
marginal_qmle <- function(data,cols=seq_len(data$p)) {
  profile <- function(rho,return_fit=FALSE) {
    g <- data$stats(rho)
    beta <- numeric(data$p)
    beta[cols] <- solve(g$R[cols,cols,drop=FALSE],g$b[cols])
    sigma2 <- marginal_rss(g,beta)/data$n
    ll <- g$ld-.5*data$n*(log(2*pi*sigma2)+1)
    if(return_fit) return(list(beta=beta[cols],rho=rho,sigma2=sigma2,loglik=ll,
      converged=TRUE,boundary=min(abs(rho-c(-.98,.98)))<1e-5))
    -ll
  }
  rho <- scalar_minimum(profile)
  profile(rho,TRUE)
}

profile_intercept <- function(g,slopes) {
  (g$b[1]-sum(g$R[1,-1]*slopes))/g$R[1,1]
}
slope_system <- function(g) {
  list(Q=g$R[-1,-1,drop=FALSE]-outer(g$R[-1,1],g$R[1,-1])/g$R[1,1],
    z=g$b[-1]-g$R[-1,1]*g$b[1]/g$R[1,1])
}
marginal_penalized <- function(data,lambda,penalty,aw,start,tolerance=2e-6,maxit=300L) {
  beta <- start$beta; rho <- start$rho; sigma2 <- start$sigma2
  old <- Inf
  for(iteration in seq_len(maxit)) {
    g <- data$stats(rho); sys <- slope_system(g)
    if(penalty=="AL") {
      update <- weighted_lasso(sys$Q,sys$z,data$n*sigma2*lambda*aw,beta[-1])
      slopes <- update$beta
    } else slopes <- as.vector(solve(sys$Q+diag(data$n*sigma2*lambda,length(aw)),sys$z))
    objective_rho <- function(r) {
      h <- data$stats(r)
      b <- c(profile_intercept(h,slopes),slopes)
      -h$ld/data$n+marginal_rss(h,b)/(2*data$n*sigma2)
    }
    rho <- optimize(objective_rho,c(-.98,.98),tol=1e-9)$minimum
    g <- data$stats(rho,TRUE)
    beta <- c(profile_intercept(g,slopes),slopes)
    sigma2 <- marginal_rss(g,beta)/data$n
    penalty_value <- if(penalty=="AL") sum(aw*abs(slopes)) else sum(slopes^2)/2
    objective <- .5*log(sigma2)-g$ld/data$n+lambda*penalty_value
    grad <- as.vector(g$R%*%beta-g$b)/(data$n*sigma2)
    if(penalty=="AL") {
      active <- abs(slopes)>1e-9
      violation <- pmax(abs(grad[-1])-lambda*aw,0)
      violation[active] <- abs(grad[-1][active]+lambda*aw[active]*sign(slopes[active]))
    } else violation <- abs(grad[-1]+lambda*slopes)
    drss <- g$da-2*sum(beta*g$db)+as.numeric(crossprod(beta,g$dR%*%beta))
    grho <- -g$dld/data$n+drss/(2*data$n*sigma2)
    if(rho<=-.98+1e-6) grho <- min(grho,0)
    if(rho>=.98-1e-6) grho <- max(grho,0)
    kkt <- max(abs(grad[1]),violation,abs(grho))
    if(kkt<tolerance && abs(old-objective)<tolerance) break
    old <- objective
  }
  list(beta=beta,rho=rho,sigma2=sigma2,
    loglik=g$ld-.5*data$n*(log(2*pi*sigma2)+1),lambda=lambda,objective=objective,
    kkt=kkt,iterations=iteration,converged=kkt<tolerance && iteration<maxit,
    boundary=min(abs(rho-c(-.98,.98)))<1e-5)
}
marginal_tune <- function(data,full,penalty,nlambda=16L) {
  aw <- 1/(abs(full$beta[-1])+data$n^(-.5))
  null <- marginal_qmle(data,1L)
  null$beta <- c(null$beta,rep(0,data$p-1L))
  if(penalty=="AL") {
    g <- data$stats(null$rho)
    gradient <- as.vector(g$R%*%null$beta-g$b)[-1]
    lambda_max <- max(abs(gradient)/(data$n*null$sigma2*aw))
    lambdas <- lambda_max*10^seq(0,-4,length.out=nlambda)
  } else {
    Q <- slope_system(data$stats(full$rho))$Q
    lambda_scale <- median(diag(Q))/(data$n*full$sigma2)
    lambdas <- lambda_scale*10^seq(2,-5,length.out=nlambda)
  }
  fits <- vector("list",nlambda+1L); previous <- null
  for(j in seq_along(lambdas)) {
    f1 <- marginal_penalized(data,lambdas[j],penalty,aw,previous)
    f2 <- marginal_penalized(data,lambdas[j],penalty,aw,full)
    candidates <- Filter(function(z) z$converged,list(f1,f2))
    if(!length(candidates)) stop(sprintf("%s candidate %d failed: KKT %.3g, %.3g",penalty,j,f1$kkt,f2$kkt))
    f <- candidates[[which.min(vapply(candidates,function(z) z$objective,numeric(1)))]]
    Q <- slope_system(data$stats(f$rho))$Q
    slope_df <- if(penalty=="AL") sum(abs(f$beta[-1])>1e-8) else
      sum(diag(solve(Q+diag(data$n*f$sigma2*lambdas[j],data$p-1L),Q)))
    f$df <- 1+slope_df; f$bic <- -2*f$loglik+(f$df+2)*log(data$n)
    f$start_gap <- abs(f1$objective-f2$objective)
    fits[[j]] <- f; previous <- f
    if(j %% 4L==0L) cat(sprintf("    %s %d/%d\n",penalty,j,nlambda))
  }
  f0 <- full; f0$lambda <- 0; f0$kkt <- 0; f0$df <- data$p; f0$start_gap <- 0
  f0$bic <- -2*f0$loglik+(data$p+2)*log(data$n)
  fits[[nlambda+1L]] <- f0
  chosen <- which.min(vapply(fits,function(z) z$bic,numeric(1)))
  f <- fits[[chosen]]
  f$retained <- c(1L,which(abs(f$beta[-1])>1e-8)+1L)
  f$edge <- chosen %in% c(1L,nlambda)
  f$path <- do.call(rbind,lapply(fits,function(z) data.frame(lambda=z$lambda,BIC=z$bic,
    df=z$df,KKT=z$kkt,converged=z$converged,start_gap=z$start_gap,boundary=z$boundary)))
  f
}
marginal_precision <- function(w,train,rho) {
  W <- Matrix(w$W,sparse=TRUE); n <- nrow(W); held <- setdiff(seq_len(n),train)
  S <- Diagonal(n)-rho*W
  Q <- crossprod(S); dQ <- -crossprod(W,S)-crossprod(S,W)
  K <- as.matrix(Q[train,train,drop=FALSE]); dK <- as.matrix(dQ[train,train,drop=FALSE])
  if(length(held)) {
    Qh <- as.matrix(Q[held,held,drop=FALSE]); F <- chol2inv(chol(Qh))
    C <- as.matrix(Q[held,train,drop=FALSE]); dC <- as.matrix(dQ[held,train,drop=FALSE])
    FC <- F%*%C
    K <- K-crossprod(C,FC)
    B <- crossprod(dC,FC)
    dK <- dK-B-t(B)+crossprod(FC,as.matrix(dQ[held,held,drop=FALSE])%*%FC)
  }
  list(K=(K+t(K))/2,dK=(dK+t(dK))/2,S=S,Q=Q,held=held)
}
marginal_information <- function(data,beta,rho,sigma2) {
  m <- marginal_precision(data$w,data$train,rho)
  H <- as.matrix(solve(m$S,data$X))
  mu <- as.vector(H%*%beta)
  dm <- as.vector(solve(m$S,data$w$W%*%mu))[data$train]
  HT <- H[data$train,,drop=FALSE]
  VdK <- solve(m$K,m$dK)
  bb <- crossprod(HT,m$K%*%HT)/sigma2
  br <- as.vector(crossprod(HT,m$K%*%dm))/sigma2
  rr <- as.numeric(crossprod(dm,m$K%*%dm))/sigma2+.5*sum(VdK*t(VdK))
  # V=K^-1, hence tr(V^-1 V')=-tr(V K').
  rs <- -.5*sum(diag(VdK))/sigma2
  ee <- matrix(c(rr,rs,rs,data$n/(2*sigma2^2)),2,2)
  be <- cbind(br,rep(0,data$p))
  C <- (bb-be%*%solve(ee,t(be)))/data$n
  C <- (C+t(C))/2
  stopifnot(min(eigen(C,symmetric=TRUE,only.values=TRUE)$values)>0,
    max(abs(bb-data$stats(rho)$R/sigma2))<1e-6)
  list(C=C,information=rbind(cbind(bb,be),cbind(t(be),ee))/data$n)
}

marginal_predict <- function(data,beta,rho) {
  # E[y_H | y_T,X,W]=Q_HH^-1{(S'X beta)_H - Q_HT y_T}.
  n <- nrow(data$X); S <- Diagonal(n)-rho*Matrix(data$w$W,sparse=TRUE)
  Q <- crossprod(S); H <- data$held
  rhs <- as.vector(crossprod(S,data$X%*%beta))[H]-
    as.vector(Q[H,data$train,drop=FALSE]%*%data$y_train)
  as.vector(solve(Q[H,H,drop=FALSE],rhs))
}

# Prediction rule for PT, S and PS 
shrinkage_weights <- function(post, alpha = 0.05) {
  q <- post$q; Tn <- post$T
  c(PT = as.numeric(q > 0 && Tn <= qchisq(1 - alpha, q)),
    S = if (q >= 3) (q - 2) / Tn else 0,
    PS = if (q >= 3) min(1, (q - 2) / Tn) else 0)
}
shrink_theta <- function(full, restricted, retained, weight) {
  uf <- numeric(length(full$beta)); uf[retained] <- restricted$beta
  list(beta = full$beta + weight * (uf - full$beta),
       rho = full$rho + weight * (restricted$rho - full$rho))
}
