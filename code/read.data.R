# This file contains the functions to read in data.

# Read TCF7L2 or CACNA1C phenotype data.
read.pheno <- function (pheno.filename, GENE, lines.to.skip = 0,
                        verbose = TRUE) {
  
  if (GENE == "TCF7L2") {
    if (verbose)
      cat(paste0("Reading raw phenotype data from: ",pheno.filename,".\n"))
    pheno.data <- read.csv(pheno.filename,comment.char = "#",
                           skip = lines.to.skip,header = TRUE,
                           check.names = FALSE,quote = "",
                           as.is = c("TCF7L2","OFTtimeofday", 
                                     "FCtimeofday","FSTtimeofday"))
    
    if (verbose)
      cat(paste0("Converting necessary variables to factors.\n"))
    pheno.data <- transform(pheno.data,
                            TCF7L2 = factor(TCF7L2, c("WT", "HET")),
                            OFTtimeofday = factor(OFTtimeofday, c("AM", "PM")),
                            FCtimeofday = factor(FCtimeofday, c("AM", "PM")),
                            FSTtimeofday = factor(FSTtimeofday, c("AM", "PM")),
                            oftbox = factor(oftbox),
                            fcbox = factor(fcbox),
                            PPIbox = factor(PPIbox),
                            FSTbucket = factor(FSTbucket))
    
    return(pheno.data)
  
  } else if (GENE == "CACNA1C") {
    
    if (verbose) cat(paste0("Reading raw phenotype data from: ", pheno.filename, ".\n"))
    pheno.data <- read.csv(pheno.filename, header = TRUE, check.names = FALSE, 
                           as.is = c("CACNA1C"))
    
    if (verbose) cat(paste0("Converting necessary variables to factors.\n"))
    pheno.data <- transform(pheno.data,
                            CACNA1C = factor(CACNA1C, c("WT", "HET")),
                            oftbox = factor(oftbox),
                            ppibox = factor(ppibox))    
    
    return(pheno.data)
    
  } else {
    stop("GENE must be set to TCF7L2 or CACNA1C.\n")
  }
  
}

# Read residuals.
read.residual <- function (resid.filename, GENE, verbose = TRUE) {
  if (verbose)
    cat(paste0("Reading in residual data from: ",resid.filename,".\n",
               "Converting genotype column to factor.\n"))
  resid.data <- read.csv(resid.filename, header = TRUE)
  resid.data[, which(colnames(resid.data) == GENE)] <- factor(resid.data[[GENE]], c("WT", "HET"))
  return(resid.data)
}

# Read combined wild-type data.
read.combined.wt <- function(combined.filename, verbose = TRUE) {
  if (verbose)
    cat(paste0("Reading WT combined data from: ", combined.filename, ".\n"))
  combined.data <- read.csv(combined.filename, header = TRUE, check.names = FALSE, 
                            as.is = c("genotype"))
  
  if (verbose) cat(paste0("Converting necessary variables to factors.\n"))
  combined.data <- transform(combined.data,
                             genotype = factor(genotype, "WT"),
                             oftbox = factor(oftbox),
                             ppibox = factor(ppibox))
  
  return(combined.data)
  
}

# Reads the Mouse Diversity Array genotypes from the CSV file, and
# returns a list containing two list elements: "geno", the n x p
# matrix of genotypes (which are all homozygous since these are
# genotypes of inbred lab strains), where n is the number of strains,
# and p is the number of SNPs; "map", a data frame giving the
# chromosome, base-pair position and alleles for each SNP.
read.mda <- function (file) {
  chromosomes <- c(1:19,"X","Y","M")
  genotypes   <- c("A","T","G","C")
  
  # Read the genotype information from the CSV file. Here I'm assuming
  # that the first 22 lines of the file are comments (lines beginning
  # with #).
  mda <- fread(file,sep = ",",header = TRUE,stringsAsFactors = FALSE,
               skip = 22,na.strings = "NA")
  
  # Discard the data.table attributes.
  class(mda) <- "data.frame"
  
  # Set the row names to the SNP ids.
  rownames(mda) <- mda$id
  mda           <- mda[-1]
  
  # Split the data frame into an n x p matrix of genotypes (where n is
  # the number of strains, and p is the number of markers), and the
  # remaining SNP information ("map").
  map  <- mda[c("chr","pos","A","B")]
  geno <- t(mda[-(1:4)])
  rm(mda)
  
  # I convert chromosome to a factor manually so that I can control
  # the order of the chromosomes in the factor. I also convert the
  # allele columns ("A" and "B") to factors.
  map <- transform(map,
                   chr = factor(chr,paste0("chr",chromosomes)),
                   A   = factor(A,genotypes),
                   B   = factor(B,genotypes))
  levels(map$chr) <- chromosomes
  
  # Set the heterozygous genotypes to missing. Note that the only
  # strain that has heterozygous genotypes is the KK/HlJ ("KK")
  # strain.
  markers <- which(geno["KK",] == "H")
  geno["KK",markers] <- NA
  
  # Return a list containing two list elements: the genotype matrix
  # ("geno"), and the SNP information ("map").
  return(list(map = map,geno = geno))
}

# ----------------------------------------------------------------------
# Reads the genotypes of the F1 mice derived from the Mouse Diversity
# Array panel. This function returns a list containing two list
# elements: "geno", the n x p matrix of genotypes, where n is the
# number of strains, and p is the number of SNPs; "map", a data frame
# giving the chromosome, base-pair position and alleles for each SNP.
read.mda.F1 <- function (file) {
  chromosomes <- c(1:19,"X","Y","M")
  genotypes   <- c("A","T","G","C")
  
  # Read the genotype information from the CSV file. Here I'm assuming
  # that the first 25 lines of the file are comments (lines beginning
  # with #).
  #   d <- fread(file,sep = ",",header = TRUE,stringsAsFactors = FALSE,
  #              skip = 25, na.strings = "NA")
  d <- fread(file,sep = ",",header = TRUE,stringsAsFactors = FALSE,
             na.strings = "NA")
  
  # Discard the data.table attributes.
  class(d) <- "data.frame"
  
  # Set the row names to the SNP ids.
  rownames(d) <- d$id
  d           <- d[-1]
  
  # Split the data frame into an n x p matrix of genotypes (where n is
  # the number of strains, and p is the number of SNPs), and the
  # remaining SNP information ("map").
  map  <- d[c("chr","pos","A1","A2")]
  geno <- t(d[-(1:4)])
  storage.mode(geno) <- "double"
  rm(d)
  
  # I convert chromosome to a factor manually so that I can control
  # the order of the chromosomes in the factor. I also convert the
  # allele columns ("A1" and "A2") to factors.
  map <- transform(map,
                   chr = factor(chr,chromosomes),
                   A1  = factor(A1,genotypes),
                   A2  = factor(A2,genotypes))
  
  # Return a list containing two list elements: the genotype matrix
  # ("geno"), and the SNP information ("map").
  return(list(map = map,geno = geno))
}









