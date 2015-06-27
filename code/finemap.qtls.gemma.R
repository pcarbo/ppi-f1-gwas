# Fine-map QTLs for phenotypes measured in wild-type F1 mice using a
# linear mixed model (GEMMA) that corrects for possible confounding to
# due relatedness

# Load packages
library(data.table)
library(qtl)

# Set location of dropbox files
dropbox.dir <- "~/Dropbox (Palmer Lab)/Palmer Lab/Kyle Engel/data_analysis/"
setwd(dropbox.dir)

# Load some extra functions
source("code/read.data.R")
source("code/data.manip.wtcombined.R")
source("code/transformation.functions.R")
source("code/wtcombined.gemma.functions.R")
source("code/misc.R")

# Set location of data files, gemma executable, and the default model information
pheno.filename <- "data/combined.wts.csv"
mda.filename   <- "data/mda.F1.lg.csv"
unc.filename   <- "data/unc.F1.lg.fixed.csv"
gemma.exe      <- "/usr/local/src/gemma/bin/gemma"
source("defaults/wtdefaults.version3.R")

# Specify the phenotype and chromosome to run finemapping on
which.analysis <- "ppi12"
chromosome     <- 2

# Should we remove markers that are polymorphic in only one strain?
extra.filter   <- TRUE

# Name of final result RData file and a directory to store some GEMMA results
resultsfile <- paste0("~/Dropbox (Palmer Lab)/Palmer Lab/Kyle Engel/wt_combined_gemma/",
                      "doublecheck.2.", which.analysis, ".chr", chromosome, ".RData")
gemmadir    <- "~/Desktop/finemap_gemma_output/"

# This is added to the diagonal of the kinship matrix to make sure that
# calculations involving this matrix are stable
delta <- 0.001

# Get the phenotype and covariates used in the QTL mapping, and the
# file for saving the results.
analysis        <- model.info[[which.analysis]]
phenotype       <- which.analysis
transformation  <- analysis$transformation
covariates      <- analysis$cov
outliers        <- analysis$outlier.function
strain.outliers <- analysis$outliers

# LOAD PHENOTYPE DATA
# -------------------
cat("Loading phenotype data.\n")
pheno <- read.combined.wt(pheno.filename)
pheno <- transform.pheno(pheno, phenotype, transformation)
pheno <- remove.outliers(pheno, phenotype, covariates, outliers, strain.outliers)

# Remove MA strain (no genotype data is available for MA)
pheno        <- subset(pheno,strain != "MA")
pheno        <- pheno[order(pheno$strain), ]
pheno$strain <- droplevels(pheno$strain)

# Only include samples in the analysis for which the phenotype and all
# covariates are observed
cols  <- c(phenotype,covariates)
rows  <- which(none.missing.row(pheno[cols]))
pheno <- pheno[rows, ]

# Print out some more information
cat("Including all (",nrow(pheno),") samples in analysis.\n",sep="")

# Convert the sex, time of day, and study columns to binary
levels(pheno$sex)         <- 0:1
levels(pheno$fctimeofday) <- 0:1
levels(pheno$study)       <- 0:1

# LOAD MDA GENOTYPE DATA
# ----------------------
cat("Loading MDA genotype data for computing kinship matrix.\n")
d    <- read.mda.F1(mda.filename)
map  <- d$map
geno <- d$geno
geno <- geno[order(rownames(geno)), ]
map  <- cbind(data.frame(snp = rownames(map)), map)
rownames(map) <- NULL
rm(d)

# Drop sex-linked and mitochondrial DNA genotypes
markers <- which(is.element(map$chr,1:19))
map     <- transform(map[markers, ], chr = droplevels(chr))
geno    <- geno[, markers]

# Keep SNPs on all chromosomes except the selected chromosome
markers <- which(map$chr != chromosome)
map     <- transform(map[markers,], chr = droplevels(chr))
geno    <- geno[, markers]

# Also, drop SNPs that are polymorphic in 5 strains or less (this was done in wtcombined.gemma.R)
markers <- which(pmin(colSums(geno), colSums(1 - geno)) > 5)
map     <- map[markers, ]
geno    <- geno[, markers]

# Align the phenotypes and MDA genotypes.
geno <- geno[match(pheno$strain, rownames(geno)), ]

# Compute the kinship matrix
cat("Computing kinship matrix.\n")
n <- nrow(pheno)
K <- tcrossprod(center.columns(geno))/length(markers)
K <- K + diag(delta, n)

# Clear unused structures
rm(map, geno, markers, n)

# LOAD UNC GENOTYPE DATA
# ----------------------
cat("Loading UNC genotype data.\n")

# Note that the read.mda.F1 function can also be used on UNC data
d    <- read.mda.F1(unc.filename)
map  <- d$map
geno <- d$geno
map  <- cbind(data.frame(snp = rownames(map)), map)
rownames(map) <- NULL
rm(d)

# Keep only the markers on the selected chromosome.
markers <- which(map$chr == chromosome)
map     <- transform(map[markers,],chr = droplevels(chr))
geno    <- geno[,markers]

# Remove markers that are polymorphic in only one strain if specified
if (extra.filter) {
  markers <- which(apply(geno, 2, function (x) sum(x) > 1))
  map <- map[markers, ]
  geno <- geno[, markers]
}

# Align the phenotypes and UNC genotypes.
geno <- geno[match(pheno$strain,rownames(geno)),]
rm(markers)

# Map QTLs using a kinship matrix estimated from the MDA SNP panel.
out          <- run.finemap(phenotype,covariates,pheno,geno,map,
                            gemmadir,gemma.exe,K)
gwscan.gemma <- out$gwscan
pve.gemma    <- out$pve

# Save the results to file.
cat("Saving results to file.\n")
save(list = c("analysis","map","gwscan.gemma","pve.gemma"),
     file = resultsfile)
