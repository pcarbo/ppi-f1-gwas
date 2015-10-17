# Maps QTLs for phenotypes measured in wild-type F1 mice using a linear
# mixed model (GEMMA) that corrects for possible confounding to due
# relatedness.

# Load packages.
require(data.table)
require(qtl)

# Load some extra functions definitions.
source("read.data.R")
source("data.manip.R")
source("functions.transform.R")
source("functions.gemma.R")
source("misc.R")

# Load the model defaults.
source("wtdefaults.v3.R")

# SCRIPT PARAMETERS
# -----------------
which.analysis <- "ppi12"
choose.samples <- "wt"
use.kinship    <- TRUE
results.file   <- "results.gemma.RData"
seed           <- 1

# Set location of the data files and the gemma executable.
pheno.filename <- "../data/pheno.csv"
geno.filename  <- "../data/mda.F1.csv"
gemmadir       <- "~/gemma_output"
gemma.exe      <- "~/bin/gemma"

# Print out the status.
cp0("On ",which.analysis," and ",choose.samples,"\n")
  
# Extract information out of model.info list.
analysis        <- model.info[[which.analysis]]
phenotype       <- which.analysis
transformation  <- analysis$transformation
covariates      <- unique(analysis$covariates)
outliers        <- analysis$outlier.function
strain.outliers <- unique(analysis$outliers)

# LOAD PHENOTYPE DATA
# -------------------
cat("Loading phenotype data.\n")
  
# Read in the data, transform it if necessary, and remove outliers
pheno <- read.combined.wt(pheno.filename)
pheno <- transform.pheno(pheno,phenotype,transformation)
pheno <- remove.outliers(pheno,phenotype,covariates,outliers,strain.outliers)

# Remove MA strain (since no genotype data is available for MA).
pheno        <- subset(pheno,strain != "MA")
pheno        <- pheno[order(pheno$strain),]
pheno$strain <- droplevels(pheno$strain)

# Only include samples in the analysis for which the phenotype and all
# covariates are observed.
cols  <- c(phenotype,covariates)
rows  <- which(none.missing.row(pheno[cols]))
pheno <- pheno[rows,]
  
# Print out some more information  
cat("Including all (",nrow(pheno),") WT samples in analysis.\n",sep="")
  
# Convert the sex, time of day, and study columns to binary values.
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
  
# Drop sex-linked and mitochondrial DNA genotypes.
markers <- which(is.element(map$chr,1:19))
map     <- transform(map[markers,],chr = droplevels(chr))
geno    <- geno[,markers]

# Also, drop SNPs that are polymorphic in 5 strains or less, because
# it is unlikely that we will discover phenotype associations with
# these SNPs.
markers <- which(pmin(colSums(geno),colSums(1 - geno)) > 5)
map     <- map[markers,]
geno    <- geno[,markers]
  
# Align the phenotypes and genotypes.
geno <- geno[match(pheno$strain,rownames(geno)),]

# Initialize the random number generator.
set.seed(seed)

# MAP QTLs USING GEMMA
# --------------------
if (use.kinship) {

  # Map QTLs using a separate kinship matrix for each chromosome
  out <- run.gemma(phenotype,covariates,pheno,geno,map,gemmadir,gemma.exe)
  gwscan.gemma <- out$gwscan
  pve.gemma    <- out$pve
} else {

  # Map QTLs without fully accounting for population structure.
  pve.gemma    <- NULL
  gwscan.gemma <- run.gemma.norr(phenotype,covariates,pheno,geno,map,
                                 gemmadir,gemma.exe)
}

# Save results to file.
cat("Saving results to file.\n")
save(list = c("analysis","choose.samples","map","gwscan.gemma","pve.gemma"),
     file = results.file)
