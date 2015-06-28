# This file contains miscellaneous functions.

# ----------------------------------------------------------------------
# Combines cat and paste into one function.
cp0 <- function (...)
  cat(paste0(...))

# ----------------------------------------------------------------------
# For each row of the matrix or data frame, returns true if all the
# entries in the row are provided (not missing).
none.missing.row <- function (x)
  rowSums(is.na(x)) == 0

# ----------------------------------------------------------------------
# Does the same thing as repmat(A,m,n) in MATLAB.
repmat <- function (A,m,n)
      return(kronecker(matrix(1,m,n),A))

# ----------------------------------------------------------------------
# Centers the columns of matrix X so that the entries in each column
# of X add up to zero.
center.columns <- function (X) {
  mu <- matrix(colMeans(X),1,ncol(X))
  X  <- X - repmat(mu,nrow(X),1)
  return(X)
}







