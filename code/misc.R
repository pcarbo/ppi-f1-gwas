# This file contains miscellaneous functions

# Combines cat and paste into one function
cp0 <- function(...) {
  cat(paste0(...))
}

# Converts a factor vector to a numeric vector
as.numeric.factor <- function(x) { as.numeric(levels(x)[x]) }

# Prints specified symbol specified times in the console (to be used to separate steps)
print.separation <- function(symbol = "_", n = 80) {
  cat(paste(rep(symbol, n), collapse = ""), "\n")
}

# Adds a row into a data frame
insert.row <- function(original, new.row, row.ind) {
  original[seq(row.ind + 1, nrow(original) + 1), ] <- original[seq(row.ind, nrow(original)), ]
  original[row.ind, ] <- new.row
  return (original)
}

# Adds a column into a data frame
insert.col <- function(original, new.col, col.ind, col.name) {
  orig.col <- ncol(original)
  original <- cbind(original, rep(-1, nrow(original)))
  colnames(original)[seq(col.ind + 1, orig.col + 1)] <- colnames(original)[seq(col.ind, orig.col)]
  original[, seq(col.ind + 1, orig.col + 1)] <- original[, seq(col.ind, orig.col)]
  original[, col.ind] <- new.col
  colnames(original)[col.ind] <- col.name
  return (original)
}

# Creates binary variables from categorical variables
binary.from.categorical <- function(categorical.vector, col.names) {
  # Creates binary factor for each 'level' of the categorical variable
  l <- list()
  for(level in levels(categorical.vector))
    l[[level]] <- factor(as.integer(categorical.vector == level))
  
  # Returns data.frame
  dat <- data.frame(l, check.names=FALSE)
  colnames(dat) <- col.names
  return(dat)
}

# Returns true if every element of given vector is numeric or integer
all.num.or.int <- function(x) {
  return(class(x) == "numeric" || class(x) == "integer")
}

# Fits a histogram to inputted data, with the number of bins chosen
# by cross validation
cross.validation.hist <- function(x, verbose = FALSE, ...) {
  n     <- length(x)
  nbins <- 1:n
  h     <- (max(x) - min(x)) / nbins
  risks <- c()
  
  for(i in 1:n) {
    currentBreaks <- seq(min(x), max(x), length=nbins[i] + 1)
    counts <- hist(x, breaks=currentBreaks, plot=FALSE)$counts
    risks[i] <- ((2 / (h[i] * (n-1))) - (((n+1) / (h[i]*(n-1))) * sum((counts/n)^2)))
  }
  loc.of.best <- which(risks == min(risks))
  
  min.risk   <- risks[loc.of.best]
  hbest      <- h[loc.of.best]
  nbins.best <- nbins[loc.of.best]
  
  if (verbose) cp0("Minimum risk: ", round(min.risk, 3), "\nOptimal h: ",
                   round(hbest, 3), "\nOptimal nbins: ", round(nbins.best, 3), "\n")
    
  h = hist(x, breaks=nbins.best)
  h$density = h$counts / sum(h$counts)
  plot(h, freq = FALSE, ...)
  
  return(invisible(0))
}

# Plots an estimate of the density of x, along with:
#   (a) a histogram of x
#   (b) a normal distribution with x's mean and sd
plot.density.estimate <- function(x, ...) {
  # Get mean and sd
  u   <- mean(x, na.rm = TRUE)
  std <- sd(x, na.rm = TRUE)
  
  # Get density estimate of x
  density.est <- density(x, kernel = "epanechnikov", na.rm = TRUE)
  
  # Find bounds to avoid cutting off graphs
  max.y <- 1.1 * max(dnorm(0, mean = u, sd = std),
                      max(density.est$y))
  
  min.x <- 1.15 * min(min(x), min(density.est$x))
  max.x <- 1.15 * max(max(x), max(density.est$x))
  
  # Make plot
  plot(density.est, ylim = c(0, max.y), xlim = c(min.x, max.x), ...)
  xbounds <- seq(1.25 * min(x, na.rm = TRUE), 1.25 * max(x, na.rm = TRUE), length=500)
  lines(xbounds, dnorm(xbounds, mean = u, sd = std), lty = "dashed", col = "blue")
  legend("topright", lty = c("solid", "dashed"), col = c("black", "blue"), 
         legend = c("Density Estimate", paste0("N(", round(u, 2), ", sd=", round(std, 2), ")")))
  
  return(invisible(0))
}

# Plot a graph of residuals vs. fitted values
plot.heteroscedastisity <- function(model) {
  # Extract information
  resid     <- model$residuals
  fitted    <- model$fitted.values
  phenotype <- as.character(model$terms[[2]])
  
  # Make plot
  plot(fitted, resid, pch = "*", xlab = "Fitted Values", ylab = "Residuals",
       main = paste0("Test for Heteroscedastisity: ", phenotype))
  
  lines(locfit(resid ~ fitted), lty = "dashed", col = "red")
  abline(0, 0, lty = 3)
  
  return(invisible(0))
}

# Plot a graph of observed values vs. fitted values
plot.linearity <- function(fitted, y.vals, phenotype) {
  plot(fitted, y.vals, pch = "*", xlab = "Fitted Values", ylab = "Observed Values",
     main = paste0("Observed vs Fitted Values: ", phenotype))
  abline(lm(y.vals ~ fitted))
  lines(locfit(y.vals ~ fitted), lty = "dashed", col = "red")
  
  return(invisible(0))
}

# Gets the phenotypes and covariates used in the analysis
get.pheno.and.covs <- function(model.info) {

  # First, get a list of phenotypes and covariates for each analysis
  n <- length(model.info)
  COVARIATES <- c()
  PHENOTYPES <- c()
  for(i in 1:n){
    current.covs  <- union(model.info[[i]]$covariates, "CACNA1C")
    current.pheno <- names(model.info)[i]
    
    PHENOTYPES[i] <- current.pheno
    COVARIATES[i] <- paste(current.covs, collapse = "+")
  }

  return(list(phenotypes = PHENOTYPES, covariates = COVARIATES))
}

# Adds significance symbols based on pvalues in a data.frame
add.sig.symbols <- function(anovaData, col){ 
  SigLevel <- c()
  for(i in 1:dim(anovaData)[1]){
    if(is.na(anovaData[i, col])){
      SigLevel[i] <- ""
    } else if(anovaData[i, col] < 0.001){ 
      SigLevel[i] <- "***"
    } else if(anovaData[i, col] < 0.01){
      SigLevel[i] <- "**"
    } else if(anovaData[i, col] < 0.05){
      SigLevel[i] <- "*"
    } else if(anovaData[i, col] < 0.1){
      SigLevel[i] <- "."
    } else {
      SigLevel[i] <- ""
    }
  }
  final <- data.frame(anovaData, SigLevel)
  colnames(final) <- c("Pr(>F)", "Sig")
  return(final)
}

# Computes the length of a vector with the option of removing NA's
len <- function(x, na.rm = TRUE) {
  if (na.rm) {
    return(sum(!is.na(x)))
  }
  return(length(x))
}

# Peter's miscellaneous functions
# Output the string using 'cat', then move the cursor back to the
# beginning of the string so that subsequent output will overwrite
# this string.
caterase <- function (s)
  cat(s,rep("\b",nchar(s)),sep="")

# ----------------------------------------------------------------------
# Create a formula for the linear model of the phenotype given the
# specified covariates.
create.lm.formula <- function (phenotype, covariates)
  return(formula(paste(phenotype,"~",paste(covariates,collapse = "+"))))

# ----------------------------------------------------------------------
# Returns TRUE if the input vector contains at least one element that
# is different from others.
is.poly <- function (x)
  return(length(unique(x[!is.na(x)])) > 1)

# ----------------------------------------------------------------------
# For each row of the matrix or data frame, returns true if all the
# entries in the row are provided (not missing).
none.missing.row <- function (x)
  rowSums(is.na(x)) == 0

# ----------------------------------------------------------------------
# Returns TRUE if x is a factor with exactly 2 levels.
is.binary.factor <- function (x)
  return(is.factor(x) & nlevels(x) == 2)

# ----------------------------------------------------------------------
# Returns a data frame with one column for each level of factor x. The
# columns of the data frame encode the categorical variable with n
# levels as n binary variables.
binary.from.categorical <- function (x, col.names = NULL) {
  
  # Create a binary factor for each value of the categorical variable.
  d <- list()
  for (value in levels(x))
    d[[value]] <- factor(as.integer(x == value))
  
  # Convert the list to a data frame, and adjust the column names, if
  # requested.
  d <- data.frame(d,check.names = FALSE)
  if (!is.null(col.names))
    names(d) <- col.names
  
  # Output the newly created data frame.
  return(d)
}

# ----------------------------------------------------------------------
# Convert a factor with exactly 2 levels to a numeric vector with
# values 0 or 1.
binfactor2num <- function (x) {
  if (!is.binary.factor(x))
    stop("Factor must have exactly 2 levels")
  return(as.numeric(x) - 1)
}

# ----------------------------------------------------------------------
# Does the same thing as repmat(A,m,n) in MATLAB.
repmat <- function (A,m,n)
  return(kronecker(matrix(1,m,n),A))

# ----------------------------------------------------------------------
# This is the same as VAR(X,1)' in MATLAB.
var1 <- function (X) {
  n <- nrow(X)
  return(apply(X,2,function(x) (n-1)/n*var(x)))
}

# ----------------------------------------------------------------------
# Returns the quadratic norm (2-norm) of vector x.
norm2 <- function (x)
  sqrt(sum(x^2))

# ----------------------------------------------------------------------
# Returns the matrix product X*X'.
matrix.square <- function (X)
  return(X %*% t(X))

# ----------------------------------------------------------------------
# Centers the columns of matrix X so that the entries in each column
# of X add up to zero.
center.columns <- function (X) {
  mu <- matrix(colMeans(X),1,ncol(X))
  X  <- X - repmat(mu,nrow(X),1)
  return(X)
}

# ----------------------------------------------------------------------
# Check whether the observed quantiles match what we would expect
# under the normal distribution.
check.normal.quantiles <- function (x) {
  
  # Discard the missing values.
  x <- x[!is.na(x)]
  
  # Transform the observations so that they have zero mean and unit
  # variance.
  x <- (x - mean(x))/sqrt(var(x))
  
  # Create a data frame giving the observed and expected proportion of
  # samples within 1, 2 and 3 standard deviations of the mean.
  return(data.frame(exp = c(pnorm(1) - pnorm(-1),
                            pnorm(2) - pnorm(-2),
                            pnorm(3) - pnorm(-3)),
                    obs = c(mean(-1 < x & x < 1),
                            mean(-2 < x & x < 2),
                            mean(-3 < x & x < 3)),
                    row.names = c("sd1","sd2","sd3")))
}

# ----------------------------------------------------------------------
# ROOTS2(A,B,C) returns solutions X to quadratic A*X^2 + B*X + C = 0.
roots2 <- function (a, b, c) {
  q <- -(b + sign(b)*sqrt(b^2 - 4*a*c))/2
  return(c(q/a,c/q))
}











