# Run this script after first running map.qtls.gemma.R up to the QTL
# mapping part.
source("polygenic.R")
X <- geno
rm(geno)

# Get the phenotype data.
y <- pheno[[phenotype]]
n <- length(y)

# Get the covariate data.
Z <- pheno[covariates]
for (col in covariates)
  if (is.factor(Z[[col]]))
    Z[[col]] <- binfactor2num(Z[[col]])
Z <- as.matrix(cbind(data.frame(intercept = rep(1,n)),Z))

# Adjust the genotypes and phenotypes so that the linear effects of
# the covariates are removed. This is equivalent to integrating out
# the regression coefficients corresponding to the covariates with
# respect to an improper, uniform prior; see Chipman, George and
# McCulloch, "The Practical Implementation of Bayesian Model
# Selection," 2001. The equivalent expressions in MATLAB are  
#
#   y = y - Z*((Z'*Z)\(Z'*y))
#   X = X - Z*((Z'*Z)\(Z'*X))  
#
# Note that this should give the same result as centering the
# columns of X and subtracting the mean from y when we have only
# one covariate, the intercept.
y <- y - c(Z %*% solve(crossprod(Z),t(y %*% Z)))
X <- X - Z %*% solve(crossprod(Z),t(Z) %*% X)

# COMPUTE POLYGENIC MODEL ESTIMATES
# ---------------------------------
# Give summary of the analysis.
h     <- seq(0.01,0.99,0.01)
out   <- polygenic.model(X,y,h)
logw  <- out$logw
sigma <- out$sigma
rm(out)

# Compute the normalized importance weights.
w <- normalizelogweights(logw)


