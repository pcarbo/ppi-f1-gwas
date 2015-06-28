# This file contains the functions to read in data.

# ----------------------------------------------------------------------
# Read combined wild-type phenotype data.
read.combined.wt <- function (file.name, verbose = TRUE) {
  if (verbose)
    cat(paste0("Reading WT combined data from: ",file.name,".\n"))
  combined.data <- read.csv(file.name,header = TRUE,check.names = FALSE, 
                            as.is = "genotype",comment.char = "#")
  if (verbose)
    cat(paste0("Converting necessary variables to factors.\n"))
  return(transform(combined.data,
                   genotype = factor(genotype,"WT"),
                   oftbox   = factor(oftbox),
                   ppibox   = factor(ppibox)))
}

# ----------------------------------------------------------------------
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
read.mda.F1 <- function (file, skip = 25) {
  chromosomes <- c(1:19,"X","Y","M")
  genotypes   <- c("A","T","G","C")
  
  # Read the genotype information from the CSV file. Here I'm assuming
  # that the first 25 lines of the file are comments (lines beginning
  # with #).
  d <- fread(file,sep = ",",header = TRUE,stringsAsFactors = FALSE,
             na.strings = "NA",skip = skip)
  
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









