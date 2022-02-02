## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  install.packages("polyqtlR")

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("Rcpp")
#  install.packages("foreach")
#  install.packages("doParallel")

## -----------------------------------------------------------------------------
library(polyqtlR)

## -----------------------------------------------------------------------------
data("phased_maplist.4x", "SNP_dosages.4x", "Phenotypes_4x")

## ----eval = FALSE-------------------------------------------------------------
#  IBD_4x <- estimate_IBD(phased_maplist = phased_maplist.4x,
#                         genotypes = SNP_dosages.4x,
#                         ploidy = 4,
#                         bivalent_decoding = FALSE,
#                         ncores = 4)

## ----eval = FALSE-------------------------------------------------------------
#  nc <- parallel::detectCores() - 1

## ----echo = FALSE-------------------------------------------------------------
nc <- 2 #to pass CRAN checks

## ---- eval = FALSE------------------------------------------------------------
#  IBD_4x <- import_IBD(folder = "TetraOrigin",
#                       filename.vec = paste0("TetraOrigin_Output_bivs_LinkageGroup",1:5,"_Summary"),
#                       bivalent_decoding = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  Importing map data under description inferTetraOrigin-Summary,Genetic map of biallelic markers
#  Importing parental phasing under description inferTetraOrigin-Summary,MAP of parental haplotypes
#  Importing IBD data under description inferTetraOrigin-Summary,Conditonal genotype probability

## ----eval=FALSE---------------------------------------------------------------
#  IBD_4x <- estimate_IBD(phased_maplist = phased_maplist.4x,
#                         genotypes = SNP_dosages.4x,
#                         method = "heur",
#                         ploidy = 4)

## ---- eval = FALSE------------------------------------------------------------
#  thinned_maplist.4x <- thinmap(maplist = phased_maplist.4x,
#                                dosage_matrix = SNP_dosages.4x)

## ---- echo = FALSE------------------------------------------------------------
cat("87 markers from a possible 93 on LG1 were included.

89 markers from a possible 93 on LG2 were included.")

## ---- out.width = "780px", echo = FALSE---------------------------------------
knitr::include_graphics("figures/thinmap.png")

## ----eval = FALSE-------------------------------------------------------------
#  IBD_4x.spl <- spline_IBD(IBD_list = IBD_4x,
#                           gap = 1) #grid at 1 cM spacing

## ----echo = FALSE-------------------------------------------------------------
IBD_4x <- readRDS(file = "IBD_4x_DR.RDS")

## ---- fig.width = 8, fig.height = 7-------------------------------------------
visualiseHaplo(IBD_list = IBD_4x,
               display_by = "name",
               linkage_group = 1,
               select_offspring = colnames(SNP_dosages.4x)[3:11], #cols 1+2 are parents
               multiplot = c(3,3)) #plot layout in 3x3 grid

## ---- echo = FALSE------------------------------------------------------------
IBDnoisy <- readRDS("noisyIBD.RDS")

## ---- fig.width = 8, fig.height = 7, echo = FALSE-----------------------------
visualiseHaplo(IBDnoisy, 
               display_by = "name", 
               select_offspring = "all",
               multiplot = c(3,3))

## ----echo = FALSE-------------------------------------------------------------
data("IBD_4x")

## -----------------------------------------------------------------------------
GIC_4x <- estimate_GIC(IBD_list = IBD_4x)

## -----------------------------------------------------------------------------
visualiseGIC(GIC_list = GIC_4x)

## -----------------------------------------------------------------------------
dim(Phenotypes_4x)

head(Phenotypes_4x)

## ----eval = FALSE-------------------------------------------------------------
#  qtl_LODs.4x <- QTLscan(IBD_list = IBD_4x,
#                         Phenotype.df = Phenotypes_4x,
#                         genotype.ID = "geno",
#                         trait.ID = "pheno",
#                         block = "year")

## ----echo = FALSE-------------------------------------------------------------
qtl_LODs.4x <- readRDS("qtl_LODs.4x.RDS")

## -----------------------------------------------------------------------------
plotQTL(LOD_data = qtl_LODs.4x,
        multiplot = c(1,2), #1 row, 2 columns
        col = "darkgreen")

## -----------------------------------------------------------------------------
plotLinearQTL(LOD_data = qtl_LODs.4x,
              col = "dodgerblue")

## ----eval = FALSE-------------------------------------------------------------
#  qtl_LODs.4x <- QTLscan(IBD_list = IBD_4x,
#                         Phenotype.df = Phenotypes_4x,
#                         genotype.ID = "geno",
#                         trait.ID = "pheno",
#                         block = "year",
#                         perm_test = TRUE,
#                         ncores = nc)

## ----echo = FALSE-------------------------------------------------------------
qtl_LODs.4x <- readRDS("qtl_LODs.4x_perm.RDS")

## -----------------------------------------------------------------------------
plotLinearQTL(LOD_data = qtl_LODs.4x,
              col = "dodgerblue")

## -----------------------------------------------------------------------------
blues <- BLUE(data = Phenotypes_4x,
              model = pheno~geno,
              random = ~1|year,
              genotype.ID = "geno")

## ----eval = FALSE-------------------------------------------------------------
#  qtl_LODs.4x <- QTLscan(IBD_list = IBD_4x,
#                         Phenotype.df = blues,
#                         genotype.ID = "geno",
#                         trait.ID = "blue",
#                         perm_test = TRUE,
#                         ncores = nc)

## ----echo = FALSE-------------------------------------------------------------
qtl_LODs.4x <- readRDS("qtl_LODs.4x_fastpermute.RDS")

## -----------------------------------------------------------------------------
plotLinearQTL(LOD_data = qtl_LODs.4x,
              col = "dodgerblue")

## -----------------------------------------------------------------------------
findPeak(qtl_LODs.4x, linkage_group = 1)

## ----eval = FALSE-------------------------------------------------------------
#  qtl_LODs.4x_cofactor <- QTLscan(IBD_list = IBD_4x,
#                                  Phenotype.df = Phenotypes_4x,
#                                  genotype.ID = "geno",
#                                  trait.ID = "pheno",
#                                  block = "year",
#                                  cofactor_df = data.frame("LG" = 1,
#                                                           "cM" = 12.3),
#                                  perm_test = FALSE,
#                                  ncores = nc)#nc is the number of cores, defined earlier

## ----echo = FALSE-------------------------------------------------------------
qtl_LODs.4x_cofactor <- readRDS("qtl_LODs.4x_cofactor.RDS")

## -----------------------------------------------------------------------------
new_pheno <- qtl_LODs.4x_cofactor$Residuals
head(new_pheno)
colnames(new_pheno)

## -----------------------------------------------------------------------------
blues <- BLUE(data = new_pheno,
              model = Pheno~geno,
              random = ~1|Block,
              genotype.ID = "geno")

## ----eval = FALSE-------------------------------------------------------------
#  qtl_LODs.4x_cofactor <- QTLscan(IBD_list = IBD_4x,
#                                  Phenotype.df = blues,
#                                  genotype.ID = "geno",
#                                  trait.ID = "blue",
#                                  perm_test = TRUE,
#                                  ncores = nc)

## ----echo = FALSE-------------------------------------------------------------
qtl_LODs.4x_cofactor <- readRDS("qtl_LODs.4x_cof_fastpermute.RDS")

## -----------------------------------------------------------------------------
plotLinearQTL(LOD_data = qtl_LODs.4x_cofactor,
              col = "red")

## ----eval = FALSE-------------------------------------------------------------
#  blues <- BLUE(data = Phenotypes_4x,
#                model = pheno~geno,
#                random = ~1|year,
#                genotype.ID = "geno")
#  
#  QTLmodel <- check_cofactors(IBD_list = IBD_4x,
#                              Phenotype.df = blues,
#                              genotype.ID = "geno",
#                              trait.ID = "blue",
#                              LOD_data = qtl_LODs.4x,
#                              ncores = nc)
#  
#  
#  QTLmodel

## ----echo = FALSE-------------------------------------------------------------
blues <- BLUE(data = Phenotypes_4x,
              model = pheno~geno,
              random = ~1|year,
              genotype.ID = "geno")
QTLmodel <- readRDS(file = "QTLmodel.RDS")
QTLmodel

## -----------------------------------------------------------------------------
PVE(IBD_list = IBD_4x,
    Phenotype.df = blues,
    genotype.ID = "geno",
    trait.ID = "blue",
    QTL_df = QTLmodel)

## -----------------------------------------------------------------------------
plotLinearQTL_list(LOD_data.ls = list(qtl_LODs.4x,
                                      qtl_LODs.4x_cofactor),
                   plot_type = "lines")

## -----------------------------------------------------------------------------
blues <- BLUE(data = Phenotypes_4x,
              model = pheno~geno,
              random = ~1|year,
              genotype.ID = "geno")

## -----------------------------------------------------------------------------
qtl_explored <- exploreQTL(IBD_list = IBD_4x,
                           Phenotype.df = blues,
                           genotype.ID = "geno",
                           trait.ID = "blue",
                           linkage_group = 1,
                           LOD_data = qtl_LODs.4x)

## -----------------------------------------------------------------------------
default_QTLconfig <- get("segList_4x",envir = getNamespace("polyqtlR"))
default_QTLconfig[1:2]

length(default_QTLconfig)

## -----------------------------------------------------------------------------
default_QTLconfig[[1]] #replace 1 with 27 to see the second most-likely model

## ----eval = FALSE-------------------------------------------------------------
#  visualiseQTLeffects(IBD_list = IBD_4x,
#                      Phenotype.df = blues,
#                      genotype.ID = "geno",
#                      trait.ID = "blue",
#                      linkage_group = 1,
#                      LOD_data = qtl_LODs.4x)

## ---- out.width = "780px", echo = FALSE---------------------------------------
knitr::include_graphics("figures/qtl_effects.png")

## -----------------------------------------------------------------------------
hist(x = blues$blue)

## ----eval = FALSE-------------------------------------------------------------
#  visualiseHaplo(IBD_list = IBD_4x,
#                 display_by = "phenotype",
#                 Phenotype.df = blues,
#                 genotype.ID = "geno",
#                 trait.ID = "blue",
#                 pheno_range = c(44,max(blues$blue)),
#                 linkage_group = 1,
#                 multiplot = c(3,3))

## ---- out.width = "780px", echo = FALSE---------------------------------------
knitr::include_graphics("figures/fish_hom1.png")

## ----eval = FALSE-------------------------------------------------------------
#  qtl_SMR.4x <- singleMarkerRegression(dosage_matrix = SNP_dosages.4x,
#                                       Phenotype.df = blues,
#                                       genotype.ID = "geno",
#                                       trait.ID = "blue",
#                                       return_R2 = TRUE,
#                                       perm_test = TRUE,
#                                       ncores = nc,
#                                       maplist = phased_maplist.4x)

## ----echo = FALSE-------------------------------------------------------------
qtl_SMR.4x <- readRDS("qtl_SMR.4x.RDS")

## -----------------------------------------------------------------------------
plotLinearQTL_list(LOD_data.ls =  list(qtl_LODs.4x,
                                       qtl_SMR.4x),
                   plot_type = c("lines","points"), #IBD results as lines, SMR results as points..
                   pch = 19,
                   colours = c("darkgreen","navyblue"))

## -----------------------------------------------------------------------------
PVE(IBD_list = IBD_4x,
    Phenotype.df = Phenotypes_4x,
    genotype.ID = "geno",
    trait.ID = "pheno",
    block = "year",
    QTL_df = data.frame("LG" = 1,"cM" = 12.3))

## ----eval = FALSE-------------------------------------------------------------
#  mr <- meiosis_report(IBD_list = IBD_4x)

## ---- out.width = "500px", echo = FALSE---------------------------------------
knitr::include_graphics("figures/mr1.png")

## ---- out.width = "500px", echo = FALSE---------------------------------------
knitr::include_graphics("figures/mr2.png")

## ----echo = FALSE-------------------------------------------------------------
mr <- readRDS("mr.RDS")

## -----------------------------------------------------------------------------
par(mfrow = c(1,2))
visualisePairing(meiosis_report.ls = mr,
                 parent = "P1",
                 datawidemax = 5)

## -----------------------------------------------------------------------------
recom.ls <- count_recombinations(IBD_4x)

## -----------------------------------------------------------------------------
layout(matrix(c(1,3,2,3),nrow=2)) #make tidy layout
RLS_summary <- plotRecLS(recom.ls) #capture the function output as well

## ---- fig.width = 8, fig.height = 7-------------------------------------------
visualiseHaplo(IBD_list = IBD_4x,
               display_by = "name",
               linkage_group = 1,
               select_offspring = colnames(SNP_dosages.4x)[3:11], #cols 1+2 are parents
               multiplot = c(3,3), #plot layout in 3x3 grid
               recombination_data = recom.ls) 

## ----eval = FALSE-------------------------------------------------------------
#  visHap1 <- visualiseHaplo(IBD_list = IBD_4x,
#                            display_by = "name",
#                            select_offspring = "all",
#                            linkage_group = 1,
#                            cM_range = c(1.56,50),
#                            recombinant_scan = c(1,3),
#                            multiplot = c(4,2))

## -----------------------------------------------------------------------------
visHap1

## ---- fig.width = 8, fig.height = 7-------------------------------------------
visualiseHaplo(IBD_list = IBD_4x,
               display_by = "name",
               select_offspring = visHap1$recombinants,
               linkage_group = 1,
               cM_range = c(1.56,50),#note: start position = first marker @1.56 cM
               recombinant_scan = c(1,3),
               multiplot = c(3,2),
               recombination_data = recom.ls) # we can add this track for clarity

## ----eval = FALSE-------------------------------------------------------------
#  visHap2 <- visualiseHaplo(IBD_list = IBD_4x,
#                            display_by = "name",
#                            select_offspring = "all",
#                            linkage_group = 1,
#                            cM_range = c(1.56,15),
#                            allele_fish = 1,
#                            multiplot = c(4,2))

## -----------------------------------------------------------------------------
visHap2

## ----fig.width = 7, fig.height = 7--------------------------------------------
visualiseHaplo(IBD_list = IBD_4x,
               display_by = "name",
               select_offspring = visHap2$allele_fish$h1,
               linkage_group = 1,
               cM_range = c(1.56,15),
               multiplot = c(4,4))

## ---- eval = FALSE------------------------------------------------------------
#  IBD_multiError <- maxL_IBD(phased_maplist = phased_maplist.4x,
#                             genotypes = SNP_dosages.4x,
#                             ploidy = 4,
#                             bivalent_decoding = FALSE,
#                             errors = c(0.01,0.02,0.05,0.1,0.2),
#                             ncores = 4)

## ---- eval = FALSE------------------------------------------------------------
#  IBD_4x <- IBD_multiError$maxL_IBD

## ---- eval = FALSE------------------------------------------------------------
#  # Copy our dosage matrix to a new matrix:
#  error_dosages <- SNP_dosages.4x
#  
#  # Only work with offspring errors for now:
#  F1cols <- 3:ncol(error_dosages)
#  
#  # Count number of (offspring) dosage scores:
#  ndose <- length(error_dosages[,F1cols])
#  
#  # Generate a vector of random positions (10% of total nr):
#  set.seed(42)
#  error.pos <- sample(1:ndose,round(ndose*0.1))
#  
#  # Replace real values with these random scores (simple error simulation):
#  error_dosages[,F1cols][error.pos] <- sapply(error_dosages[,F1cols][error.pos], function(x) sample(setdiff(0:4,x),1))
#  
#  # Re-estimate IBD probabilities, using the maxL_IBD function:
#  IBD_multiError <- maxL_IBD(phased_maplist = phased_maplist.4x,
#                             genotypes = error_dosages,
#                             ploidy = 4,
#                             errors = c(0.01, 0.02, 0.05, 0.1, 0.3),
#                             ncores = nc)
#  
#  IBD_4x.err <- IBD_multiError$maxL_IBD
#  
#  # Re-impute marker dosages
#  new_dosages <- impute_dosages(IBD_list=IBD_4x.err,dosage_matrix=error_dosages,
#                                min_error_prior = 0.01,
#                                rounding_error = 0.2)
#  

## ---- echo = FALSE------------------------------------------------------------
cat("____________________________________________________________________________

186 out of a total possible 186 markers were imputed
In the original dataset there were 0 % missing values among the 186 markers
In the imputed dataset there are now 1.4 % missing values among these, using a rounding threshold of 0.2

____________________________________________________________________________
The % of markers with changed dosage scores is as follows (0 = no change):
error.matrix
    0     1     2     3     4 
89.44  4.18  2.88  1.72  0.35 ")
new_dosages <- readRDS(file = "new_dosages.RDS")
error_dosages <- readRDS(file = "error_dosages.RDS")
F1cols <- 3:ncol(error_dosages); ndose <- length(error_dosages[,F1cols])

## -----------------------------------------------------------------------------
# Check the results:
round(length(which(error_dosages[,F1cols] != SNP_dosages.4x[,F1cols])) / ndose,2)

round(length(which(new_dosages[,F1cols] != SNP_dosages.4x[,F1cols])) / ndose,2)

