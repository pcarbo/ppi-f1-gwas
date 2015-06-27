# Maps QTLs for phenotypes measured in wild-type F1 mice using a linear
# mixed model (GEMMA) that corrects for possible confounding to due
# relatedness

# Load packages
require(data.table)
require(qtl)

 # Set location of dropbox files
dropbox.dir <- "/home/palmerlab/Dropbox (Palmer Lab)/Palmer Lab/Kyle Engel/data_analysis/"
setwd(dropbox.dir)

# Load some extra functions
source("code/read.data.R")
source("code/data.manip.wtcombined.R")
source("code/transformation.functions.R")
source("code/wtcombined.gemma.functions.R")
source("code/misc.R")

# Set location of data files, gemma executable, and the default model information
pheno.filename <- "data/combined.wts.csv"
geno.filename  <- "data/mda.F1.lg.csv"
gemma.exe      <- "/usr/local/src/gemma/bin/gemma"
source("defaults/wtdefaults.version3.R")

system("mkdir ~/Desktop/gemma_output")

for (which.analysis in names(model.info)) {
  
  # Since this is on wtcombined data, we only want WTs
  choose.samples <- "wt"
  
  # Create folders for GEMMA and specify some model information
  system(paste0("mkdir ~/Desktop/gemma_output/", which.analysis))
  gemmadir       <- paste0("~/Desktop/gemma_output/", which.analysis, "/", choose.samples)
  system(paste0("mkdir ", gemmadir))
  
  # Name of final result RData file
  resultsfile <- paste0("~/Dropbox (Palmer Lab)/Palmer Lab/Kyle Engel/wt_combined_gemma/",
                        which.analysis, ".wtcombined.RData")
  
  # A few more specifications
  use.kinship <- TRUE
  seed        <- 1
  
  # Print out status
  cp0("On ", which.analysis, " and ", choose.samples, "\n")
  
  # Extract information out of the model.info list
  analysis         <- model.info[[which.analysis]]
  phenotype        <- which.analysis
  transformation   <- analysis$transformation
  covariates       <- unique(analysis$covariates)
  outliers         <- analysis$outlier.function
  strain.outliers  <- unique(analysis$outliers)

  # LOAD PHENOTYPE DATA
  # -------------------
  cat("Loading phenotype data.\n")
  
  # Read in the data, transform it if necessary, and remove outliers
  pheno <- read.combined.wt(pheno.filename)
  pheno <- transform.pheno(pheno, phenotype, transformation)
  pheno <- remove.outliers(pheno, phenotype, covariates, outliers, strain.outliers)
  
  # Remove MA strain (no genotype data is available for MA)
  pheno        <- subset(pheno, strain != "MA")
  pheno        <- pheno[order(pheno$strain), ]
  pheno$strain <- droplevels(pheno$strain)
  
  # Only include samples in the analysis for which the phenotype and all
  # covariates are observed
  cols  <- c(phenotype,covariates)
  rows  <- which(none.missing.row(pheno[cols]))
  pheno <- pheno[rows, ]
  
  # Print out some more information  
  cat("Including all (", nrow(pheno), ") WT samples in analysis.\n", sep="")
  
  # Convert the sex, time of day, and study columns to binary
  levels(pheno$sex)         <- 0:1
  levels(pheno$fctimeofday) <- 0:1
  levels(pheno$study)       <- 0:1
  
  # LOAD GENOTYPE DATA
  # ------------------
  cat("Loading genotype data.\n")
  d    <- read.mda.F1(geno.filename)
  map  <- d$map
  geno <- d$geno
  geno <- geno[order(rownames(geno)), ]
  map  <- cbind(data.frame(snp = rownames(map)),map)
  rownames(map) <- NULL
  rm(d)
  
  # Drop sex-linked and mitochondrial DNA genotypes
  markers <- which(is.element(map$chr,1:19))
  map     <- transform(map[markers, ], chr = droplevels(chr))
  geno    <- geno[, markers]
  
  # Also, drop SNPs that are polymorphic in 5 strains or less, because
  # it is unlikely that we will discover phenotype associations with
  # these SNPs
  markers <- which(pmin(colSums(geno), colSums(1 - geno)) > 5)
  map     <- map[markers, ]
  geno    <- geno[, markers]
  
  # Align the phenotypes and genotypes
  geno <- geno[match(pheno$strain, rownames(geno)), ]
  
  # Initialize the random number generator
  set.seed(seed)

  # MAP QTLs USING GEMMA
  # --------------------
  if (use.kinship) {

    # Map QTLs using a separate kinship matrix for each chromosome
    out <- run.gemma(phenotype, covariates, pheno, geno, map, gemmadir, gemma.exe)
    gwscan.gemma <- out$gwscan
    pve.gemma    <- out$pve
    
  } else {

    # Maps QTLs without fully accounting for population structure
    pve.gemma    <- NULL
    gwscan.gemma <- run.gemma.norr(phenotype, covariates, pheno, geno, map,
                                   gemmadir, gemma.exe)
    
  }

  # Save results to file
  cat("Saving results to file.\n")
  save(list = c("analysis", "choose.samples", "map", "gwscan.gemma", "pve.gemma"),
       file = resultsfile)

}
