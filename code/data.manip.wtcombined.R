# This file contains functions for manipulating data in the analysis

# Transforms the given phenotype with the given transformation function
# Returns an updated prepared.pheno with the phenotype's old
# column replaced by a transformed version
transform.pheno <- function(prepared.pheno, phenotype, transformation, verbose = TRUE) {
  
  if (verbose) cp0("\n*Starting transform.pheno on ", phenotype, ".\n")
  
  # Check if there is a transformation to do
  # Return the original prepared.pheno if there isn't
  if (is.null(transformation)) {
    if (verbose) cp0("Did not transform ", phenotype, " because there is no specified transformation.\n")
    return(prepared.pheno)
  }
  # Perform the transformation
  if (verbose) cp0("Transformed ", phenotype, ".\n")
  prepared.pheno[[phenotype]] <- transformation(prepared.pheno[[phenotype]])
  
  # Return the updated phenotype data
  return(prepared.pheno)
  
}

# Removes outliers in a given phenotype based on the given outlier function
# Note that the function assumes transformed data is given
# This function removes outliers based on the outlier function and based on outliers given by ID
# Returns an updated transformed.pheno with a NA assigned to each outlying observation
remove.outliers <- function(transformed.pheno, phenotype, covariates,
                            outlier.function, outliers, verbose = TRUE) {
    
  if (verbose) cp0("\n*Starting remove.outliers on ", phenotype, ".\n")
  
  # Check if there is an outlier function, if there isn't then don't do anything
  # If there is, apply it and set outliers to NA
  if (is.null(outlier.function)) {
    if (verbose) cp0("Did not apply an outlier removal function on ", phenotype, " (none given).\n")
  } else { # There is an outlier function
    # Create a function that tells if x is an outlier (takes into account NAs)
    is.outlier <- function(x) {
      y <- outlier.function(x)
      y[is.na(y)] <- FALSE
      return(y)
    }
    
    # Now fit a model and save the residuals as r
    if (length(covariates) > 0) {
      f <- formula(paste(phenotype, "~", paste(covariates, collapse = " + ")))
      r <- resid(lm(f, transformed.pheno, na.action = na.exclude))
    } else {
      r <- transformed.pheno[[phenotype]]
    }
    
    # Check the residuals for outliers and report the number that will be removed
    if (verbose) {
      n <- sum(is.outlier(r))
      if (n == 0) {
        cp0("No outliers for ", phenotype)
      } else {
        cp0("Removed ", n, " outliers for ", phenotype)
      }
      if (length(covariates) > 0) {
        cp0(" conditioned on [", paste(covariates, collapse = " + "), "]")
      }
      cp0(" (outlier removal function).\n")    
    }
    
    # Assign a "NA" to outliers
    transformed.pheno[is.outlier(r), phenotype] <- NA
  }
  
  # Check if there are other outliers specified to be removed
  if (length(outliers) > 0) {
    to.remove <- match(outliers, transformed.pheno$id)
    transformed.pheno[to.remove, phenotype] <- NA
    
    # Report amount removed
    if (verbose) {
      cp0("Removed ", length(to.remove), " more outliers for ", phenotype)
      if (length(covariates) > 0) {
        cp0(" conditioned on [" ,paste(covariates, collapse = " + "), "]")
      }
      cp0(" (ID specified).\n")
    }
    
  } else { # No ID specified outliers to remove
    if (verbose) cp0("Did not remove any ID specified outliers (none given).\n")
  }

  return(transformed.pheno)
}



