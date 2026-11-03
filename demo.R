# Small example: one simulated SAR data set, all seven estimators. 

source("sar_functions.R")

# Settings: 10 x 10 queen lattice, 4 target and 8 candidate covariates
set.seed(2026)
n <- 100; rho <- 0.5; p1 <- 4; p2 <- 8
beta <- c(1.5, -1, 2, 0.8, rep(0, p2))
w <- make_weights(n, "queen")
X <- make_design(n, p1 + p2, rho_x = 0.3)
y <- spatial_solve(w, as.vector(X %*% beta) + rnorm(n), rho)

# OF and the effective information C
full <- sar_qmle(y, X, w)
C <- sar_information(X, full$beta, full$rho, full$sigma2, w)$C

# Adaptive Lasso selects the submodel; UF refits it
al <- tune_penalty(y, X, w, full, "AL", nlambda = 16)
ridge <- tune_penalty(y, X, w, full, "RIDGE", nlambda = 16)
restricted <- sar_qmle(y, X[, al$selected, drop = FALSE], w)

# PT, S and PS for the retained coefficients
post <- post_shrinkage(full, restricted, al$selected, C, n)
est <- c(post$estimates, list(AL = al$beta, RIDGE = ridge$beta))

cat(sprintf("rho_OF = %.3f, selected = %s, q = %d, T_n = %.3f\n",
            full$rho, paste(al$selected, collapse = ","), post$q, post$T))
print(shrinkage_weights(post))
tab <- sapply(est, function(b) b[1:p1])
rownames(tab) <- paste0("beta", 1:p1)
tab <- rbind(tab, "Squared error" = colSums((tab - beta[1:p1])^2))
print(round(cbind(truth = c(beta[1:p1], NA), tab), 3))
