# Fine-map QTLs for phenotypes measured in wild-type F1 mice using a
# linear mixed model (GEMMA) that corrects for possible confounding to
# due relatedness

# Load packages.
library(data.table)
library(qtl)

# Load some extra functions definitions.
source("read.data.R")
source("data.manip.R")
source("functions.transform.R")
source("functions.gemma.R")
source("misc.R")

# Load the model defaults.
source("wtdefaults.v3.R")

# Set location of the data files and the gemma executable.
pheno.filename <- "../data/pheno.csv"
mda.filename   <- "../data/mda.F1.csv"
unc.filename   <- "../data/unc.F1.csv"
gemmadir       <- "~/gemma_output"
gemma.exe      <- "~/bin/gemma"

# SCRIPT PARAMETERS
# -----------------
chromosome     <- 7
which.analysis <- "ppi12"
resultsfile    <- "results.finemap.gemma.RData"

# Should we remove markers that are polymorphic in only one strain?
extra.filter <- TRUE

# This is added to the diagonal of the kinship matrix to make sure that
# calculations involving this matrix are stable.
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
pheno <- transform.pheno(pheno,phenotype,transformation)
pheno <- remove.outliers(pheno,phenotype,covariates,outliers,strain.outliers)

# Remove MA strain (no genotype data is available for MA).
pheno        <- subset(pheno,strain != "MA")
pheno        <- pheno[order(pheno$strain),]
pheno$strain <- droplevels(pheno$strain)

# Only include samples in the analysis for which the phenotype and all
# covariates are observed.
cols  <- c(phenotype,covariates)
rows  <- which(none.missing.row(pheno[cols]))
pheno <- pheno[rows,]

# Print out some more information.
cat("Including all (",nrow(pheno),") samples in analysis.\n",sep="")

# Convert the sex, time of day, and study columns to binary.
levels(pheno$sex)         <- 0:1
levels(pheno$fctimeofday) <- 0:1
levels(pheno$study)       <- 0:1

# LOAD MDA GENOTYPE DATA
# ----------------------
cat("Loading MDA genotype data for computing kinship matrix.\n")
d    <- read.mda.F1(mda.filename)
map  <- d$map
geno <- d$geno
geno <- geno[order(rownames(geno)),]
map  <- cbind(data.frame(snp = rownames(map)),map)
rownames(map) <- NULL
rm(d)

# Drop sex-linked and mitochondrial DNA genotypes.
markers <- which(is.element(map$chr,1:19))
map     <- transform(map[markers,],chr = droplevels(chr))
geno    <- geno[,markers]

# Keep SNPs on all chromosomes except the selected chromosome.
markers <- which(map$chr != chromosome)
map     <- transform(map[markers,],chr = droplevels(chr))
geno    <- geno[,markers]

# Also, drop SNPs that are polymorphic in 5 strains or less.
markers <- which(pmin(colSums(geno),colSums(1 - geno)) > 5)
map     <- map[markers,]
geno    <- geno[,markers]

# Align the phenotypes and MDA genotypes.
geno <- geno[match(pheno$strain,rownames(geno)),]

# Compute the kinship matrix.
cat("Computing kinship matrix.\n")
n <- nrow(pheno)
K <- tcrossprod(center.columns(geno)) / length(markers)
K <- K + diag(delta,n)

# Clear unused variables.
rm(map,geno,markers,n)

# LOAD UNC GENOTYPE DATA
# ----------------------
cat("Loading UNC genotype data.\n")

# Note that the same read.mda.F1 function can also be used to read in
# the UNC genotype data.
d    <- read.mda.F1(unc.filename,skip = 32)
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
  markers <- which(apply(geno,2,function (x) sum(x) > 1))
  map     <- map[markers,]
  geno    <- geno[,markers]
}

# Align the phenotypes and UNC genotypes.
geno <- geno[match(pheno$strain,rownames(geno)),]
rm(markers)

# Map QTLs using a kinship matrix estimated from the MDA SNP panel.
out          <- run.gemma(phenotype,covariates,pheno,geno,map,
                          gemmadir,gemma.exe,K = K)
gwscan.gemma <- out$gwscan
pve.gemma    <- out$pve

# Save the results to file.
cat("Saving results to file.\n")
save(list = c("analysis","map","gwscan.gemma","pve.gemma"),
     file = resultsfile)
