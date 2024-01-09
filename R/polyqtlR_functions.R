#' Calculate Best Linear Unbiased Estimates using linear mixed model from \code{nlme} package
#' @description Calculation of BLUEs from data frame of genotype names and phenotypes (assuming repeated measurements)
#' @param data Data frame of genotype codes and corresponding phenotypes
#' @param model The model specification of fixed terms, eg. Yield ~ Clones
#' @param random The random component of the model (repeat structure, can be nested), eg. ~1 | Blocks if only Blocks are used
#' @param genotype.ID The colname used to describe genotypes, e.g. "Clones"
#' @return A data-frame with columns "geno" for the genotype names, and "blue" for the BLUEs.
#' @examples
#' data("Phenotypes_4x")
#' blue <- BLUE(data = Phenotypes_4x,model = pheno~geno,random = ~1|year,genotype.ID = "geno")
#' @export
BLUE <- function(data,
                 model,
                 random,
                 genotype.ID) {
  
  if (!requireNamespace("nlme", quietly = TRUE)) {
    stop("Package \"nlme\" needed for this function to work. Please install it first.",
         call. = FALSE)
  }
  
  lmeRES <- nlme::lme(fixed = model, data = data, random = random, na.action = na.omit) #At the moment I omit missing values.. later I could add imputation here.
  
  blue <- lmeRES$coef$fixed
  blue <- c(blue[1], blue[-1] + blue[1]) #Add the intercept term, then remove it
  
  ## This causes an error if e.g. genotype.ID = "F1" and all offspring are F1_001 etc:
  # names(blue)[2:length(blue)] <- gsub(genotype.ID, "", names(blue)[2:length(blue)])
  
  n1 <- nchar(genotype.ID)
  n2 <- nchar(names(blue)[2:length(blue)])
  
  names(blue)[2:length(blue)] <- substr(x = names(blue)[2:length(blue)],
                                        start = n1 + 1,
                                        stop = n2)
  names(blue)[1] <- setdiff(unique(data[,genotype.ID]),names(blue)[2:length(blue)])
  
  return(data.frame("geno" = names(blue),
                    "blue" = blue))
} #BLUE


#' Build a multi-QTL model using step-wise procedure of checking genetic co-factors.
#' @description The function \code{check_cofactors} initially fits all significant QTL positions as co-factors, both individually and in combination. Significance thresholds
#' are re-estimated each time, yielding threshold-corrected LOD scores. If this leads to a change in the estimated position of QTL, or detection of subsequent peaks, a second 
#' round of co-factor inclusion is performed for all new QTL or novel QTL combinations. Finally, the multi-QTL model that maximises the individual significance of each
#' QTL is returned as a data.frame. This can be directly passed to the function \code{\link{PVE}} to estimate the percentage variance explained by the full
#' multi-QTL model and all possible sub-models. 
#' Note: this function estimates the most likely QTL positions by maximising the threshold-corrected LOD at QTL peaks.
#' Non-additive interactions between QTL may be missed as a result. It is recommended to run a manual co-factor analysis as well,
#' as described in the package vignette.
#' @param IBD_list List of IBD_probabilities as estimated using one of the various methods available (e.g. \code{\link{estimate_IBD}}).
#' @param Phenotype.df A data.frame containing phenotypic values
#' @param genotype.ID The colname of \code{Phenotype.df} that contains the population identifiers (F1 names) (must be a colname of \code{Phenotype.df})
#' @param trait.ID The colname of \code{Phenotype.df} that contains the response variable to use in the model (must be a colname of \code{Phenotype.df})
#' @param LOD_data Output of \code{\link{QTLscan}} function. Since v.0.1.0 this argument is optional - function will re-run a \code{QTLscan} if not provided. Indeed, it may be desirable
#' to not specify \code{LOD_data} if argument \code{test_full_model} is \code{TRUE}, as this will first combine the best results using additive-effects or allelic interaction-effects models before
#' searching for additional QTL.
#' @param min_res The minimum genetic distance (resolution) assumed possible to consider 2 linked QTL (on the same linkage group) as independent. By default a value of 20 cM is used.
#' This is not to suggest that 20 cM is a realistic resolution in a practical mapping study, but it provides the function with a criterion to consider 2 significant QTL within this distance as one and the same. 
#' For this purpose, 20 cM seems a reasonable value to use. In practice, closely linked QTL will generally "explain" all the variation at nearby positions, making it unlikely to
#' be able to disentangle their effects. QTL positions will vary slightly when co-factors are introduced, but again this variation is presumed not to exceed 20 cM either side.
#' @param test_full_model By default \code{FALSE}, in which case the normal additive-effects model is used in \code{\link{QTLscan}}. If set to \code{TRUE},
#' then both the additive and full models are run for each genome-wide scan. 
#' @param verbose Logical, by default \code{TRUE} - should progress messages be printed to the console?
#' @param \dots Option to pass extra arguments to \code{QTLscan}, for example specifying \code{ncores} for parallel processing, or changing
#' the default settings of the permutation test (by default the number of permutations to perform = 1000 and alpha = 0.05). For a full list of options see the documentation of \code{\link{QTLscan}}.
#' @return Data frame with the following columns:
#' \describe{
#' \item{LG}{Linkage group identifier}
#' \item{cM}{CentiMorgan position}
#' \item{deltaLOD}{The difference between the LOD score at the peak and the significance threshold (always positive, otherwise the QTL would not be significant)}
#' \item{CofactorID}{An identifier giving the co-factor model used in detecting the QTL (if no co-factors were included then \code{NA}). The co-factor model is described
#' by concatenating all co-factor positions with a '+', so for example 1_10+4_20 would mean a co-factor model with 2 positions included as co-factors, namely 10 cM on linkage
#' group 1 and 20 cM on linkage group 4.}
#' }
#' @examples 
#' data("IBD_4x","BLUEs.pheno","qtl_LODs.4x")
#' check_cofactors(IBD_list=IBD_4x,Phenotype.df=BLUEs.pheno,
#' genotype.ID="Geno",trait.ID="BLUE",LOD_data=qtl_LODs.4x)
#' @export
check_cofactors <- function(IBD_list,
                            Phenotype.df,
                            genotype.ID,
                            trait.ID,
                            LOD_data = NULL,
                            min_res = 20,
                            test_full_model = FALSE,
                            verbose = TRUE,
                            ...){
  
  if(is.null(LOD_data)){ #The user has not provided output from a QTL scan

    if(test_full_model){
      
      if(verbose) message(paste0("No LOD_data detected. Running an initial QTL scan (with argument allelic_interaction = TRUE and FALSE)..."))
      
      LOD_data1 <- QTLscan(IBD_list = IBD_list,
                          Phenotype.df = Phenotype.df,
                          genotype.ID = genotype.ID,
                          trait.ID = trait.ID,
                          perm_test = TRUE,
                          allelic_interaction = FALSE,
                          verbose = verbose,
                          ...)
      
      LOD_data2 <- QTLscan(IBD_list = IBD_list,
                           Phenotype.df = Phenotype.df,
                           genotype.ID = genotype.ID,
                           trait.ID = trait.ID,
                           perm_test = TRUE,
                           allelic_interaction = TRUE,
                           verbose = verbose,
                           ...)
      
      
      ## Combine the best results of additive and non-additive scans
      qtl.df1 <- findQTLpeaks(LOD_data1)
      qtl.df2 <- findQTLpeaks(LOD_data2)
      combinedQTL <- rbind(qtl.df1,qtl.df2)
      
      qtl.df <- do.call(rbind,lapply(sort(unique(combinedQTL$LG)), function(lg){
        temp <- combinedQTL[combinedQTL$LG == lg,]
        
        out <- temp[which.max(temp$deltaLOD),]
        
        if(max(abs(out$cM - temp$cM)) >= min_res){
          rest <- temp[-which.max(temp$deltaLOD),]
          rest <- rest[abs(rest$cM - out$cM) >= min_res,]
          if(nrow(rest) >= 1){
            out <- rbind(out,rest[which.max(rest$deltaLOD),])
          }
        }
        
        return(out)
      }))
      
    } else{
      ## Only run the additive-effects scan
      if(verbose) message(paste0("No LOD_data detected. Running an initial QTL scan (with argument allelic_interaction = ",test_full_model,")..."))
      
      LOD_data <- QTLscan(IBD_list = IBD_list,
                          Phenotype.df = Phenotype.df,
                          genotype.ID = genotype.ID,
                          trait.ID = trait.ID,
                          perm_test = TRUE,
                          allelic_interaction = test_full_model,
                          verbose = verbose,
                          ...)
      
      qtl.df <- findQTLpeaks(LOD_data)
    }

  } else{
    if(is.null(LOD_data$Perm.res)) stop("This function only fits significant QTL positions as co-factors. Please re-run QTLscan with perm_test = TRUE to check significance.")
    qtl.df <- findQTLpeaks(LOD_data) #qtl.df <- polyqtlR:::findQTLpeaks(LOD_data) #debugging
  }
  
  ## Check that there are significant QTL to proceed:
  if(nrow(qtl.df) == 0) stop("This function only fits significant QTL positions as co-factors. If you still want to add non-significant positions as co-factors, use function QTLscan instead.")
  
  if(verbose) message(paste0("There ",ifelse(nrow(qtl.df) == 1,"was ","were "),nrow(qtl.df)," significant QTL detected in initial scan, proceeding with these..."))
  
  ## Generate all possible combinations of QTL:
  QTLmodel.list <- lapply(1:nrow(qtl.df), function(n) combn(1:nrow(qtl.df),n))
  
  if(verbose) message("Searching for extra QTL given detected co-factors....")
  
  ## Define an internal function for fitting cofactor, used 4 times subsequently:
  fit_cofactor <- function(ibds,
                           phenos,
                           genoID,
                           traitID,
                           cofdf,
                           ai){
    
    temp.QTL <- QTLscan(IBD_list = ibds,
                        Phenotype.df = phenos,
                        genotype.ID = genoID,
                        trait.ID = traitID,
                        cofactor_df = cofdf,
                        allelic_interaction = ai,
                        verbose = FALSE,
                        ...)
    
    resids <- temp.QTL$Residuals
    
    temp.QTL <- QTLscan(IBD_list = ibds,
                        Phenotype.df = resids,
                        genotype.ID = genoID,
                        trait.ID = "Pheno",
                        perm_test = TRUE,
                        allelic_interaction = ai,
                        verbose = FALSE,
                        ...)
    
    temp_qtl.df <- findQTLpeaks(temp.QTL) #temp_qtl.df <- polyqtlR:::findQTLpeaks(temp.QTL) #debugging
    
    return(temp_qtl.df)
  }
  
  fullQTL <- do.call(rbind,lapply(QTLmodel.list, function(q_mat = QTLmodel.list[[1]]){
    do.call(rbind,apply(q_mat,2,function(sub_q = 1){
      
      temp_qtl.df <- fit_cofactor(ibds = IBD_list,
                                  phenos = Phenotype.df,
                                  genoID = genotype.ID,
                                  traitID = trait.ID,
                                  cofdf = qtl.df[sub_q,c("LG","cM","model")],
                                  ai = FALSE) #run first with allelic_interaction = FALSE
      
      # if(nrow(temp_qtl.df) > 0) temp_qtl.df$CofactorID <- paste(paste(qtl.df$LG[sub_q],round(qtl.df$cM[sub_q],2),"a",sep="_"),collapse="+")
      if(nrow(temp_qtl.df) > 0) temp_qtl.df$CofactorID <- paste(paste(qtl.df$LG[sub_q],round(qtl.df$cM[sub_q],2),qtl.df$model[sub_q],sep="_"),collapse="+")
      
      ## Runs if user would like QTL model with allelic interactions to also be checked (v.0.1.0):
      if(test_full_model){
        temp_qtl.df2 <- fit_cofactor(ibds = IBD_list,
                                     phenos = Phenotype.df,
                                     genoID = genotype.ID,
                                     traitID = trait.ID,
                                     cofdf = qtl.df[sub_q,c("LG","cM","model")],
                                     ai = TRUE) #re-run with allelic_interaction = TRUE
        
        
        # if(nrow(temp_qtl.df2) > 0) temp_qtl.df2$CofactorID <- paste(paste(qtl.df$LG[sub_q],round(qtl.df$cM[sub_q],2),"f",sep="_"),collapse="+")
        if(nrow(temp_qtl.df2) > 0) temp_qtl.df2$CofactorID <- paste(paste(qtl.df$LG[sub_q],round(qtl.df$cM[sub_q],2),qtl.df$model[sub_q],sep="_"),collapse="+")
        
        temp_qtl.df <- rbind(temp_qtl.df,temp_qtl.df2)
      }
      
      return(temp_qtl.df)
    }))
  })
  )
  
  ## Add the original QTL positions as well:
  qtl.df$CofactorID <- NA 
  fullQTL <- rbind(fullQTL,qtl.df) 
  
  ## Need to reduce this down to best positions
  refinedQTL <- do.call(rbind,lapply(sort(unique(fullQTL$LG)), function(lg){
    temp <- fullQTL[fullQTL$LG == lg,]
    
    out <- temp[which.max(temp$deltaLOD),]
    
    if(max(abs(out$cM - temp$cM)) >= min_res){
      rest <- temp[-which.max(temp$deltaLOD),]
      rest <- rest[abs(rest$cM - out$cM) >= min_res,]
      if(nrow(rest) >= 1){
        out <- rbind(out,rest[which.max(rest$deltaLOD),])
      }
    }
    
    return(out)
  }))
  
  if(nrow(refinedQTL) > nrow(qtl.df)){
    ## There were new QTL discovered in the process, run the procedure again 
    if(verbose) message("New QTL detected: attempting to refine QTL positions further....")
    
    overlap <- NULL
    for(r in 1:nrow(refinedQTL)){
      hit <- refinedQTL[r,"LG"] == qtl.df$LG & refinedQTL[r,"cM"] == qtl.df$cM
      if(any(hit)) overlap <- c(overlap,r)
    }
    
    QTLmodel.list2 <- lapply(1:nrow(refinedQTL), function(n) {
      out <- combn(1:nrow(refinedQTL),n)
      if(length(overlap) > 0){
        rem <- apply(out,2,function(x) all(x%in%overlap))
        if(any(rem)) out <- out[,-which(rem),drop = FALSE]  
      }
      return(out)
    })
    
    refinedQTL2 <- do.call(rbind,lapply(QTLmodel.list2, function(q_mat){
      do.call(rbind,apply(q_mat,2,function(sub_q){
        
        temp_qtl.df <- fit_cofactor(ibds = IBD_list,
                                    phenos = Phenotype.df,
                                    genoID = genotype.ID,
                                    traitID = trait.ID,
                                    cofdf = refinedQTL[sub_q,],
                                    ai = FALSE)
        
        
        if(nrow(temp_qtl.df) > 0) temp_qtl.df$CofactorID <- paste(paste(refinedQTL$LG[sub_q],round(refinedQTL$cM[sub_q],2),"a",sep="_"),collapse="+")
        
        ## Runs if user would like QTL model with allelic interactions to also be checked (v.0.1.0):
        if(test_full_model){
          temp_qtl.df2 <- fit_cofactor(ibds = IBD_list,
                                       phenos = Phenotype.df,
                                       genoID = genotype.ID,
                                       traitID = trait.ID,
                                       cofdf = refinedQTL[sub_q,],
                                       ai = TRUE)
          
          
          if(nrow(temp_qtl.df2) > 0) temp_qtl.df2$CofactorID <- paste(paste(refinedQTL$LG[sub_q],round(refinedQTL$cM[sub_q],2),"f",sep="_"),collapse="+")
          
          temp_qtl.df <- rbind(temp_qtl.df,temp_qtl.df2)
        }
        
        return(temp_qtl.df)
      }))
    }))
    
    refinedQTL2 <- rbind(refinedQTL2,refinedQTL)
    
    outputQTL <- do.call(rbind,lapply(sort(unique(refinedQTL2$LG)), function(lg){
      temp <- refinedQTL2[refinedQTL2$LG == lg,]
      
      out <- temp[which.max(temp$deltaLOD),]
      
      if(max(abs(out$cM - temp$cM)) >= min_res){
        rest <- temp[-which.max(temp$deltaLOD),]
        rest <- rest[abs(rest$cM - out$cM) >= min_res,]
        if(nrow(rest) >= 1){
          out <- rbind(out,rest[which.max(rest$deltaLOD),])
        }
      }
      
      return(out)
    }))
    
  } else{
    outputQTL <- refinedQTL
  }
  
  rownames(outputQTL) <-  paste0("Q",1:nrow(outputQTL))
  
  ## Add a column of QTL numbers (linking with PVE function numbering)
  # outputQTL$QTL <- paste0("Q",1:nrow(outputQTL))
  
  return(outputQTL)
  
} #check_cofactors


#' Predict recombination breakpoints using IBD probabilities
#' @description The function \code{count_recombinations} returns a list of all predicted recombination breakpoints. The output can be passed 
#' using the argument \code{recombination_data} to the function \code{\link{visualiseHaplo}}, where the predicted breakpoints overlay the haplotypes. 
#' Alternatively, a genome-wide visualisation of the recombination landscape both per linkage group and per individual can be generated using the function \code{\link{plotRecLS}}, 
#' which can be useful in identifying problematic areas of the linkage maps, or problematic individuals in the population. Currently, recombination break-points
#' are only estimated from bivalents in meiosis; any offspring resulting from a predicted multivalent is excluded from the analysis and will be returned with a \code{NA} value.
#' @param IBD_list List of IBD_probabilities as estimated using one of the various methods available (e.g. \code{\link{estimate_IBD}}).
#' @param plausible_pairing_prob The minimum probability of a pairing configuration needed to analyse an individual's IBD data.
#' The default setting of 0.3 accommodates scenarios where e.g. two competing plausible pairing scenarios are possible.
#' In such situations, both pairing configurations (also termed "valencies") would be expected to have a probability close to 0.5. Both are then considered,
#' and the output contains the probability of both situations. These can then be used to generate a probabilistic recombination landscape. In some cases,
#' it may not be possible to discern the pairing in one of the parents due to a lack of recombination (ie. full parental haplotypes were transmitted). In such cases, 
#' having a lower threshold here will allow more offspring to be analysed without affecting the quality of the predictions. If a more definite 
#' set of predictions is required, simply increase \code{plausible_pairing_prob} to eliminate such uncertainty. These individuals will then be 
#' returned with a \code{NA} value.
#' In any case, it is always helpful to visualise the output using the function \code{\link{visualiseHaplo}}.
#' @return A nested list corresponding to each linkage group. Within each LG, a list with 3 items is returned, specifying the \code{plausible_pairing_prob}, the \code{map} and 
#' the predicted \code{recombinations} in each individual in the mapping population. Per individual, all valencies with a probability greater than
#' \code{plausible_pairing_prob} are returned, specifying both the \code{Valent_probability} and the best estimate of the cM position of the
#' \code{recombination_breakpoints} involving pairs of homologues A, B, C etc. (in the order parent 1, parent 2). 
#' If no recombinations are predicted, a \code{NA} value is given instead.
#' @examples 
#' data("IBD_4x")
#' recom.ls <- count_recombinations(IBD_4x)
#' @export
count_recombinations <- function(IBD_list,
                                 plausible_pairing_prob = 0.3) {
  
  IBD_list <- test_IBD_list(IBD_list)
  
  IBDarray <- IBD_list[[1]]$IBDarray
  IBDtype <- IBD_list[[1]]$IBDtype
  bivalent_decoding <- IBD_list[[1]]$biv_dec
  
  if(!bivalent_decoding) warning("Currently recombinations are predicted from bivalents only.")
  
  gap <- IBD_list[[1]]$gap
  
  ploidy <- IBD_list[[1]]$ploidy
  ploidy2 <- IBD_list[[1]]$ploidy2
  
  if(IBDtype == "genotypeIBD"){
    # if(IBD_list[[1]]$method == "hmm_TO"){
    #   ## We are using the output of TetraOrigin
    #   GenCodes <- t(sapply(IBD_list[[1]]$genocodes,function(x) as.numeric(strsplit(as.character(x),"")[[1]])))
    #   Nstates <- nrow(GenCodes)
    #   indicatorMatrix <- matrix(0,nrow=Nstates, ncol=ploidy + ploidy2)
    #   for(r in 1:nrow(indicatorMatrix)){
    #     for(h in 1:ncol(GenCodes)){
    #       indicatorMatrix[r,GenCodes[r,h]] <- indicatorMatrix[r,GenCodes[r,h]] + 1
    #     }
    #   }
    #   
    # } else{
    ## We are using the output of estimate_IBD
    GenCodes <- t(sapply(sapply(IBD_list[[1]]$genocodes,function(x) strsplit(x,"")), function(y) sapply(y,function(z) which(letters == z))))
    Nstates <- nrow(GenCodes)
    indicatorMatrix <- matrix(0,nrow=Nstates, ncol=ploidy + ploidy2)
    for(r in 1:nrow(indicatorMatrix)){
      for(h in 1:ncol(GenCodes)){
        indicatorMatrix[r,GenCodes[r,h]] <- indicatorMatrix[r,GenCodes[r,h]] + 1
      }
    }
    
    # }
    
  } else if(IBDtype == "haplotypeIBD"){
    Nstates <- ploidy + ploidy2 #this is not actually needed..
    indicatorMatrix <- NULL
  } else{
    stop("Unable to determine IBD type, please check IBD_list elements have a valid $IBDtype tag.")
  }
  
  out.ls <- setNames(lapply(seq(length(IBD_list)), function(ch){
    nind <- dim(IBD_list[[ch]]$IBDarray)[3]
    
    if(!bivalent_decoding){
      # ML <- IBD_list[[ch]]$marginal.likelihood
      # ML <- sapply(unname(IBD_list[[ch]]$marginal.likelihoods),function(x) setNames(max(x),rownames(x)[which.max(x)]))
      
      # f1ok <- names(ML) == "BB" #simplest for now to take only BB offspring rather than messing with half-bits of offspring
      f1ok <- IBD_list[[ch]]$valency == "BB"
      
      if(length(which(f1ok)) == 0) stop(paste("Cannot proceed on LG",ch,"- no offspring predicted from bivalent-only meioses!"))
    } else{
      f1ok <- rep(TRUE,nind)
    }
    
    pair.ls <- IBD_list[[ch]]$pairing
    
    if(is.null(gap)){
      cMpositions <- IBD_list[[ch]]$map$position
    } else{
      cMpositions <- as.numeric(sapply(strsplit(dimnames(IBD_list[[ch]]$IBDarray)[[1]],"cM"),
                                       "[",2))
    }
    
    F1_recombinants <- lapply(1:nind, function(f1 = 9){
      
      if(f1 %in% which(f1ok)){
        ## Need to weight by the probabilities of the pairing; choose only plausible pairings:
        plausible <- pair.ls[[f1]][pair.ls[[f1]] >= plausible_pairing_prob]
        
        ## Generate haplotype probabilities
        haploprobs <- IBD_list[[ch]]$IBDarray[,,f1] %*% indicatorMatrix
        colnames(haploprobs) <- LETTERS[1:ncol(haploprobs)]
        
        f1.recoms <- setNames(lapply(names(plausible), function(pair = names(plausible)) {
          prob.pair <- strsplit(pair,"-")[[1]]
          
          recoms <- setNames(lapply(prob.pair, function(valent = prob.pair[2]){
            val.homs <- strsplit(valent,"")[[1]]
            
            ## For now this should be in a bivalents-context
            if(length(val.homs) != 2) stop("Currently only bivalents considered!")
            # round(haploprobs[,3:4],2)
            v1pred <- haploprobs[,val.homs[1]] > haploprobs[,val.homs[2]] #assumes only 2 columns
            
            rec.vec <- v1pred[-length(v1pred)] != v1pred[-1]
            rec.pos <- rowMeans(cbind(cMpositions[which(rec.vec)],
                                      cMpositions[which(rec.vec) + 1]))
            if(length(rec.pos) == 0) rec.pos <- NA
            
            return(rec.pos)
          }),prob.pair)
          
          return(list("Valent_probability" = round(plausible[[pair]],3),
                      "recombination_breakpoints" = recoms))
          
        }),names(plausible))
        
        if(length(f1.recoms) == 0) f1.recoms <- NA
        
      } else{
        f1.recoms <- NA
      }
      
      return(f1.recoms)
    })
    
    return(list("recombinations" = setNames(F1_recombinants, dimnames(IBD_list[[ch]]$IBDarray)[[3]]),
                "map" = IBD_list[[ch]]$map,
                "plausible_pairing_prob" = plausible_pairing_prob))
  }), names(IBD_list))
  
  return(out.ls)
} #count_recombinations


#' Function to extract the phased map from a mappoly.map object
#' @description Convert MAPpoly.map object into a phased maplist, needed for IBD estimation
#' @param mappoly_object An object of class 'mappoly.map', for example output of the function \code{mappoly::est_rf_hmm_sequential}
#' @return A phased.maplist, with linkage group names LG1 etc. Each list item is a data.frame with columns marker, position followed
#' by the phased map, coded in 1 and 0 for presence/absence of SNP (alternative) allele on parental homologues (h) numbered 1:ploidy for parent 1
#' and ploidy + 1 : 2*ploidy for parent 2.
#' @examples
#' \dontrun{
#' library("mappoly")
#' phased.maplist <- convert_mappoly_to_phased.maplist(maps.hexafake)
#' }
#' @export
convert_mappoly_to_phased.maplist <- function(mappoly_object){
  
  if(!inherits(mappoly_object, what = "list")) stop("list input expected!")
  
  ploidy <- mappoly_object[[1]]$info$ploidy
  
  ## The mappoly update 0.2.3 -> 0.3.0 changed this object's name from $m to $ploidy. 
  ## Try to add backward compatibility with a test on ploidy for people still running mappoly v0.2.3 or earlier: 
  if(length(ploidy) != 1 || !is.numeric(ploidy)){
    ploidy <- mappoly_object[[1]]$info$m #this test should be redundant but failed CRAN testing r-oldrel-macos-x86_64
  }
  
  nLG <- length(mappoly_object)
  
  out.ls <- lapply(mappoly_object, function(x){
    
    p1.phase <- do.call(rbind, lapply(x$maps[[1]]$seq.ph$P, function(x){
      temp <- rep(0,ploidy)
      temp[x[x>0]]<-1
      return(temp)
    }))
    
    p2.phase <- do.call(rbind, lapply(x$maps[[1]]$seq.ph$Q, function(x){
      temp <- rep(0,ploidy)
      temp[x[x>0]]<-1
      return(temp)
    }))
    
    phase <- cbind(p1.phase, p2.phase)
    colnames(phase) <- paste0("h",1:(2*ploidy))
    
    r <- x$maps[[1]]$seq.rf
    
    # Use Haldane's mapping function (cM)
    d <- -50*log(1-2*r) 
    
    return(cbind(data.frame("marker" = x$info$mrk.names,
                            "position" = round(cumsum(c(0, d)),3)),phase))
  })
  
  names(out.ls) <- paste0("LG",1:nLG)
  
  return(out.ls)
} #convert_mappoly_to_phased.maplist


#' Estimate the Genotypic Information Coefficient (GIC)
#' @description Function to estimate the GIC per homologue using IBD probabilities
#' @param IBD_list List of IBD probabilities
#' @return A nested list; each list element (per linkage group) contains the following items:
#' \describe{
#' \item{GIC : }{Matrix of GIC values estimated from the IBD probabilities}
#' \item{map : }{Integrated linkage map positions of markers used in IBD calculation}
#' \item{parental_phase : }{The parental marker phasing, coded in 1 and 0's}
#' }
#' @examples
#' data("IBD_4x")
#' GIC_4x <- estimate_GIC(IBD_list = IBD_4x)
#' @export
estimate_GIC <- function(IBD_list){
  
  IBD_list <- test_IBD_list(IBD_list)
  
  ploidy <- IBD_list[[1]]$ploidy
  ploidy2 <- IBD_list[[1]]$ploidy2
  
  ## Check the list depth:
  # listdepth <- list.depth(IBD_list) #polyqtlR:::list.depth(IBD_list)
  # if(listdepth != 3) stop("Unexpected input - IBD_list is expected to be a nested list representing 1 or more chromosomes! Please check input.")
  
  IBDarray <- IBD_list[[1]]$IBDarray
  IBDtype <- IBD_list[[1]]$IBDtype
  bivalent_decoding <- IBD_list[[1]]$biv_dec
  gap <- IBD_list[[1]]$gap
  
  popSize <- dim(IBDarray)[3]
  
  if(IBDtype == "genotypeIBD"){
    
    if(!bivalent_decoding) warning("Cannot estimate the GIC from multivalent-derived offspring due to possible double reduction.
                                   \nScreening out such individuals and proceeding...")
    
    # if(IBD_list[[1]]$method == "hmm_TO"){
    #   ## We are using the output of TetraOrigin
    #   # Nstates <- Nstates.fun(biv_dec = bivalent_decoding, pl = ploidy, pl2 = ploidy2)
    #   # mname <- paste0("GenotypeMat",Nstates)
    #   # indicatorMatrix <- as.matrix(get(mname, envir = getNamespace("polyqtlR")))
    #   
    #   GenCodes <- t(sapply(IBD_list[[1]]$genocodes,function(x) as.numeric(strsplit(as.character(x),"")[[1]])))
    #   Nstates <- nrow(GenCodes)
    #   indicatorMatrix <- matrix(0,nrow=Nstates, ncol=ploidy + ploidy2)
    #   for(r in 1:nrow(indicatorMatrix)){
    #     for(h in 1:ncol(GenCodes)){
    #       indicatorMatrix[r,GenCodes[r,h]] <- indicatorMatrix[r,GenCodes[r,h]] + 1
    #     }
    #   }
    #   
    # } else{
    ## We are using the output of estimate_IBD
    GenCodes <- t(sapply(sapply(IBD_list[[1]]$genocodes,function(x) strsplit(x,"")), function(y) sapply(y,function(z) which(letters == z))))
    Nstates <- nrow(GenCodes)
    indicatorMatrix <- matrix(0,nrow=Nstates, ncol=ploidy + ploidy2)
    for(r in 1:nrow(indicatorMatrix)){
      for(h in 1:ncol(GenCodes)){
        indicatorMatrix[r,GenCodes[r,h]] <- indicatorMatrix[r,GenCodes[r,h]] + 1
      }
    }
    
    # }
    
    
  } else if(IBDtype == "haplotypeIBD"){
    
    Nstates <- ploidy + ploidy2 #this is not actually needed..
    indicatorMatrix <- NULL
    
  } else{
    stop("Unable to determine IBD type, please check IBD_list elements have a valid $IBDtype tag.")
  }
  
  GIC.ls <- lapply(seq(length(IBD_list)), function(l){
    
    IBDarray <- IBD_list[[l]]$IBDarray
    
    if(IBDtype == "genotypeIBD"){
      
      if(!bivalent_decoding){
        # ML <- IBD_list[[l]]$marginal.likelihood
        # f1ok <- names(ML) == "BB" #simplest for now to take only BB offspring rather than messing with half-bits of offspring
        f1ok <- IBD_list[[l]]$valency == "BB"
        nf1 <- length(which(f1ok))
        if(nf1 == 0) stop(paste("Cannot proceed on LG",l,"- no offspring predicted from bivalent-only meioses!"))
      } else{
        nf1 <- popSize
        f1ok <- rep(TRUE,nf1)
      }
      
      IBDhaps <- array(NA, dim=c(dim(IBDarray)[1], dim(indicatorMatrix)[2], nf1))
      IBDhaps[] <- apply(IBDarray[,,f1ok],3,function(x) x %*% indicatorMatrix)
      
      dimnames(IBDhaps)[[1]] <- dimnames(IBDarray)[[1]]
      dimnames(IBDhaps)[[2]] <- paste0("h", 1:(ploidy + ploidy2))
      dimnames(IBDhaps)[[3]] <- dimnames(IBDarray)[[3]][f1ok]
      
      GICs <- do.call(rbind, lapply(1:dim(IBDhaps)[1], function(p) {
        apply(IBDhaps[p,,], 1, function(x) 1 - 4*sum(x*(1-x))/nf1)
      }))
      
    } else{
      
      GICs <- do.call(rbind, lapply(1:dim(IBDarray)[1], function(p) {
        apply(IBDarray[p,,], 1, function(x) 1 - 4*sum(x*(1-x))/popSize)
      }))
      
    }
    
    if(!is.null(gap)){
      cMpositions <- as.vector(sapply(dimnames(IBDarray)[[1]], function(name) as.numeric(unlist(strsplit(name,"cM"))[2])))
    } else{
      cMpositions <- round(IBD_list[[l]]$map$position,2)
    }
    
    rownames(GICs) <- cMpositions
    # colnames(GICs) <- paste0("h",1:(ploidy + ploidy2))
    
    return(list(GIC = GICs, 
                map = IBD_list[[l]]$map, 
                parental_phase = IBD_list[[l]]$parental_phase,
                ploidy = IBD_list[[l]]$ploidy,
                ploidy2 = IBD_list[[l]]$ploidy2
    ))
  }
  )
  
  names(GIC.ls) <- names(IBD_list)
  
  return(GIC.ls)
} #estimate_GIC


#' Generate IBD probabilities from marker genotypes and a phased linkage map
#' @description \code{estimate_IBD} is a function for creating identity-by-descent (IBD) probabilities. Two computational methods are offered:
#' by default IBD probabilites are estimated using hidden Markov models, but a heuristic method based on Bourke et al. (2014) is also included.
#' Basic input data for this function are marker genotypes (either discrete marker dosages (ie scores 0, 1, ..., ploidy representing the number of copies of the marker allele),
#' or the probabilities of these dosages) and a phased linkage map. Details on each of the methods are included under \code{method}
#' @param input_type Can be either one of 'discrete' or 'probabilistic'. For the former (default), \code{dosage_matrix} must be supplied,
#' while for the latter \code{probgeno_df} must be supplied. Note that probabilistic genotypes can only be accepted if the \code{method} is default ('hmm').
#' @param genotypes Marker genotypes, either a 2d matrix of integer marker scores or a data.frame of dosage probabilities. 
#' Details are as follows:
#' \describe{
#' \item{discrete : }{
#'  If \code{input_type} is 'discrete', \code{genotypes} is a matrix of marker dosage scores with markers in rows and individuals in columns.
#'  Both (marker) rownames and (individual or sample) colnames are needed.
#' }
#' \item{probabilistic : }{
#' If \code{input_type} is 'probabilistic', \code{genotypes} is a data frame as read from the scores file produced by function \code{saveMarkerModels} of R package 
#' \code{fitPoly}, or alternatively, a data frame containing at least the following columns:
#' \describe{
#' \item{SampleName : }{
#' Name of the sample (individual)
#' }
#' \item{MarkerName : }{
#' Name of the marker
#' }
#' \item{P0 : }{
#' Probabilities of dosage score '0'
#' }
#' \item{P1, P2... etc. : }{
#' Probabilities of dosage score '1' etc. (up to max offspring dosage, e.g. P4 for tetraploid population)
#' }
#' }
#' }
#' }
#' @param phased_maplist A list of phased linkage maps, the output of \code{polymapR::create_phased_maplist}
#' @param method The method used to estimate IBD probabilities, either \code{"hmm"} or \code{"heur"}. By default, the Hidden Markov Model (hmm) method is used.
#' This uses an approach developed by Zheng et al (2016), and implemented in the 'TetraOrigin' package. However, unlike the original TetraOrigin software, it does not 
#' re-estimate parental linkage phase, as this is assumed to have been generated during map construction. Alternatively, a heuristic algorithm can be employed (\code{method = "heur"}), providing
#' computational efficiency at higher ploidy levels (hexaploid, octoploid etc.), but at the cost of some accuracy. 
#' If \code{method = "hmm"} is specified, only diploid, triploid, autotetraploid and autohexaploid populations are currently allowed, while \code{method = "heur"} 
#' caters for all possible ploidy levels. Furthermore, the argument \code{bivalent_decoding} can only be set to \code{FALSE} in the case of 
#' the 'hmm' method (i.e. allowing for the possibility of multivalent formation and double reduction).
#' @param remove_markers Optional vector of marker names to remove from the maps. Default is \code{NULL}.
#' @param ploidy Integer. Ploidy of the organism.
#' @param ploidy2 Optional integer, by default \code{NULL}. Ploidy of parent 2, if different from parent 1.
#' @param parent1 Identifier of parent 1, by default assumed to be \code{"P1"}
#' @param parent2 Identifier of parent 2, by default assumed to be \code{"P2"}
#' @param individuals By default "all" offspring are included, but otherwise a subset can be selected, using a vector of offspring indexing numbers (1,2, etc.)
#' according to their order in \code{dosage_matrix}
#' @param log Character string specifying the log filename to which standard output should be written. If \code{NULL} log is send to stdout.
#' @param map_function Mapping function to use when converting map distances to recombination frequencies.
#' Currently only \code{"haldane"} or \code{"kosambi"} are allowed.
#' @param bivalent_decoding Option to consider only bivalent pairing during formation of gametes (ignored for diploid populations, as only bivalents possible there), by default \code{TRUE}
#' @param error The (prior) probability of errors in the offspring dosages, usually assumed to be small but non-zero
#' @param full_multivalent_hexa Option to allow multivalent pairing in both parents at the hexaploid level, by default \code{FALSE}. Note that if \code{TRUE},
#' a very large available RAM may be required (>= 32Gb) to process the data. 
#' @param verbose Logical, by default \code{TRUE}. Should progress messages be written?
#' @param ncores How many CPU cores should be used in the evaluation? By default 1 core is used.
#' @param fix_threshold If \code{method = "heur"}, the threshold to fix the IBD probabilities while correcting for the sum of probabilities.
#' @param factor_dist If \code{method = "heur"}, the factor by which to increase or decrease the recombination frequencies as calculated from the map distances.
#' @return A list of IBD probabilities, organised by linkage group (as given in the input \code{phased_maplist}). Each
#' list item is itself a list containing the following:
#' \describe{
#' \item{IBDtype}{The type of IBD; for this function only "genotypeIBD" are calculated.}
#' \item{IBDarray}{A 3d array of IBD probabilities, with dimensions marker, genotype-class and F1 individual.}
#' \item{map}{A 3-column data-frame specifying chromosome, marker and position (in cM)}
#' \item{parental_phase}{Phasing of the markers in the parents, as given in the input \code{phased_maplist}}
#' \item{marginal.likelihoods}{A list of marginal likelihoods of different valencies if method "hmm" was used, otherwise \code{NULL}}
#' \item{valency}{The predicted valency that maximised the marginal likelihood, per offspring. For method "heur", \code{NULL}}
#' \item{offspring}{Offspring names}
#' \item{biv_dec}{Logical, whether bivalent decoding was used in the estimation of the F1 IBD probabilities.}
#' \item{gap}{The size of the gap (in cM) used when interpolating the IBD probabilities. See function \code{\link{spline_IBD}} for details.}
#' \item{genocodes}{Ordered list of genotype codes used to represent different genotype classes.}
#' \item{pairing}{log likelihoods of each of the different pairing scenarios considered (can be used e.g. for post-mapping check of preferential pairing)}
#' \item{ploidy}{ploidy of parent 1}
#' \item{ploidy2}{ploidy of parent 2}
#' \item{method}{The method used, either "hmm" (default) or "heur". See argument \code{method}}
#' \item{error}{The error prior used, if method "hmm" was used, otherwise \code{NULL}}
#' }
#' @references
#' \itemize{
#' \item{Durbin R, Eddy S, Krogh A, Mitchison G (1998) Biological sequence analysis: Probabilistic models of proteins and nucleic acids. Cambridge: Cambridge University Press.}
#' \item{Hackett et al. (2013) Linkage analysis and QTL mapping using SNP dosage data in a tetraploid potato mapping population. PLoS One 8(5): e63939}
#' \item{Zheng et al. (2016) Probabilistic multilocus haplotype reconstruction in outcrossing tetraploids. Genetics 203: 119-131}
#' \item{Bourke P.M. (2014) QTL analysis in polyploids: Model testing and power calculations. Wageningen University (MSc thesis)}
#' }
#' @examples
#' data("phased_maplist.4x", "SNP_dosages.4x")
#' estimate_IBD(phased_maplist=phased_maplist.4x,genotypes=SNP_dosages.4x,ploidy=4)
#' @export
estimate_IBD <- function(input_type = "discrete",
                         genotypes,
                         phased_maplist,
                         method = "hmm",
                         remove_markers = NULL,
                         ploidy,
                         ploidy2 = NULL,
                         parent1 = "P1",
                         parent2 = "P2",
                         individuals = "all",
                         log = NULL,
                         map_function = "haldane",
                         bivalent_decoding = TRUE,
                         error = 0.01,
                         full_multivalent_hexa = FALSE,
                         verbose = FALSE,
                         ncores = 1,
                         fix_threshold = 0.1,
                         factor_dist = 1
){
  method <- match.arg(method, c("hmm","heur"))
  
  if(method == "hmm"){
    out <- hmm_IBD(input_type = input_type,
                   genotypes = genotypes,
                   phased_maplist = phased_maplist,
                   remove_markers = remove_markers,
                   ploidy = ploidy,
                   ploidy2 = ploidy2,
                   parent1 = parent1,
                   parent2 = parent2,
                   individuals = individuals,
                   log = log,
                   map_function = map_function,
                   bivalent_decoding = bivalent_decoding,
                   error = error,
                   full_multivalent_hexa = full_multivalent_hexa,
                   verbose = verbose,
                   ncores = ncores)
  } else{
    if(input_type != "discrete") stop("Probabilistic genotypes are only accepted when method = hmm")
    
    out <- fast_IBD(phased_maplist = phased_maplist,
                    dosage_matrix = genotypes,
                    map_function = map_function,
                    ploidy = ploidy,
                    ploidy2 = ploidy2,
                    fix_threshold = fix_threshold,
                    factor_dist = factor_dist,
                    ncores = ncores)
  }
  
  return(out)
} #estimate_IBD


#' Explore the possible segregation type of a QTL peak using Schwarz Information Criterion
#' @description Function to explore the possible segregation type at a QTL position using the Schwarz Information Criterion
#' @param IBD_list List of IBD probabilities
#' @param Phenotype.df A data.frame containing phenotypic values
#' @param genotype.ID The colname of \code{Phenotype.df} that contains the population identifiers (F1 names) (must be a colname of \code{Phenotype.df})
#' @param trait.ID The colname of \code{Phenotype.df} that contains the response variable to use in the model (must be a colname of \code{Phenotype.df})
#' @param linkage_group Numeric identifier of the linkage group being tested, based on the order of \code{IBD_list}.
#' Only a single linkage group is allowed.
#' @param LOD_data  Output of \code{\link{QTLscan}} function
#' @param cM By default \code{NULL}, in which case the position of maximum LOD score is taken as the position of interest. 
#' Otherwise, the cM position to be explored.
#' @param QTLconfig Nested list of homologue configurations and modes of action of QTL to be explored and compared, the output of \code{\link{segMaker}}.
#' Note that a default List is available of all possible bi-allelic QTL if none is provided.
#' Each list element is itself a list with components
#' \describe{
#'  \item{homs : }{ a vector of length at least 1, describing the proposed homologues the functional allele Q is on}
#'  \item{mode : }{ Vector of same length as \code{homs} with codes "a" for additive and "d" for dominant.}
#' }
#' @param plotBIC Logical, with default \code{TRUE} - should the calculated BIC values be plotted?
#' @param deltaBIC Numeric, by default 6. Configurations within this distance of the minimum BIC are considered plausible.
#' @param testAllele_Effects  Logical, with default \code{TRUE} - should the effects of the different alleles be tested
#' using the most likely QTL configuration?
#' @param log Character string specifying the log filename to which standard output should be written. If \code{NULL} log is send to stdout.
#' @return List with the following items:
#' \describe{
#' \item{linkage_group}{Linkage group of the QTL peak being explored}
#' \item{cM}{CentiMorgan position of the locus being explored}
#' \item{BIC}{Vector of BIC values corresponding to elements of \code{QTLconfig} provided for testing}
#' \item{Allele.effects}{Summary of the means and standard errors of groups with (+)
#' and without(-) the specified allele combinations for the most likely QTLconfig
#' if \code{testAllele_Effects} = \code{TRUE} (\code{NULL} otherwise).}
#' \item{genotype.means}{A one-column matrix of mean phenotype values of offspring classes, with rownames
#' corresponding to the genotype class. If the probability of certain genotype classes is 0 (e.g. double reduction
#' classes where no double reduction occurred), then the genotype mean for that class will be \code{NA}}
#' }
#' @examples
#' data("IBD_4x","BLUEs.pheno","qtl_LODs.4x")
#' exploreQTL(IBD_list = IBD_4x,
#'            Phenotype.df = BLUEs.pheno,
#'            genotype.ID = "Geno",
#'            trait.ID = "BLUE",
#'            linkage_group = 1,
#'            LOD_data = qtl_LODs.4x)
#' @export
exploreQTL <- function(IBD_list,
                       Phenotype.df,
                       genotype.ID,
                       trait.ID,
                       linkage_group,
                       LOD_data,
                       cM = NULL,
                       QTLconfig = NULL,
                       plotBIC = TRUE,
                       deltaBIC = 6,
                       testAllele_Effects = TRUE,
                       log = NULL){
  # Test input:
  IBD_list <- test_IBD_list(IBD_list)
  
  if(length(linkage_group)!=1 | !is.numeric(linkage_group))
    stop("linkage_group should be a single numeric identifier.")
  
  if(!is.numeric(deltaBIC)) stop("deltaBIC should be a numeric value (typically 6 or more)")
  
  ploidy <- IBD_list[[1]]$ploidy
  ploidy2 <- IBD_list[[1]]$ploidy2
  
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  ploidyF1 <- (ploidy + ploidy2)/2
  
  if(is.null(QTLconfig)){
    
    if(ploidyF1 == 3){ #need a work-around for triploid populations
      if(ploidy > ploidy2){
        sname <- "segList_3x" #this was the original version
      } else{
        sname <- "segList_3x_24" #update since version 0.1.0
      }
    } else{
      sname <- paste0("segList_",ploidyF1,"x")
    }
    
    QTLconfig <- get(sname)
  }
  
  QTLmode <- sapply(seq(length(QTLconfig)),function(n) {
    match.arg(QTLconfig[[n]]$mode,c("a","d"))
  })
  
  QTLseg <- lapply(seq(length(QTLconfig)),function(n) {
    unique(match.arg(arg = as.character(QTLconfig[[n]]$homs),
                     choices = as.character(1:(ploidy + ploidy2)),
                     several.ok = TRUE))
  })
  
  add_vect <- order(QTLmode)[QTLmode[order(QTLmode)] == "a"] 
  dom_vect <- order(QTLmode)[QTLmode[order(QTLmode)] == "d"]
  
  ## If the analysis is being done on a subset of linkage groups, this can cause issues:
  IBDarray <- IBD_list[[linkage_group]]$IBDarray
  IBDtype <- IBD_list[[linkage_group]]$IBDtype
  bivalent_decoding <- IBD_list[[1]]$biv_dec
  gap <- IBD_list[[1]]$gap
  
  if(IBDtype == "genotypeIBD"){
    
    # if(IBD_list[[1]]$method == "hmm_TO"){
    #   ## We are using the output of TetraOrigin
    #   # Nstates <- Nstates.fun(biv_dec = bivalent_decoding, pl = ploidy, pl2 = ploidy2)
    #   # mname <- paste0("GenotypeMat",Nstates)
    #   # indicatorMatrix <- as.matrix(get(mname, envir = getNamespace("polyqtlR")))
    #   # GenotypeCodes <- get(paste0("GenotypeCodes",Nstates))
    #   
    #   GenotypeCodes <- t(sapply(IBD_list[[1]]$genocodes,function(x) as.numeric(strsplit(as.character(x),"")[[1]])))
    #   Nstates <- nrow(GenotypeCodes)
    #   indicatorMatrix <- matrix(0,nrow=Nstates, ncol=ploidy + ploidy2)
    #   for(r in 1:nrow(indicatorMatrix)){
    #     for(h in 1:ncol(GenotypeCodes)){
    #       indicatorMatrix[r,GenotypeCodes[r,h]] <- indicatorMatrix[r,GenotypeCodes[r,h]] + 1
    #     }
    #   }
    #   
    #   
    # } else{
    ## We are using the output of estimate_IBD
    GenotypeCodes <- t(sapply(sapply(IBD_list[[1]]$genocodes,function(x) strsplit(x,"")), function(y) sapply(y,function(z) which(letters == z))))
    Nstates <- nrow(GenotypeCodes)
    indicatorMatrix <- matrix(0,nrow=Nstates, ncol=ploidy + ploidy2)
    for(r in 1:nrow(indicatorMatrix)){
      for(h in 1:ncol(GenotypeCodes)){
        indicatorMatrix[r,GenotypeCodes[r,h]] <- indicatorMatrix[r,GenotypeCodes[r,h]] + 1
      }
    }
    # }
    
  } else if(IBDtype == "haplotypeIBD"){
    Nstates <- ploidy + ploidy2 #this is not actually needed..
    indicatorMatrix <- NULL
  } else{
    stop("Unable to determine IBD type, please check IBD_list elements have a valid $IBDtype tag.")
  }
  
  # listdepth <- list.depth(IBD_list)#polyqtlR:::list.depth(IBD_list)
  # if(listdepth != 3) stop("Unexpected input - IBD_list is expected to be a nested list representing 1 or more chromosomes! Please check input.")
  
  if(dim(IBDarray)[2] %% Nstates != 0) stop("Incompatible input detected..")
  
  if(!is.null(gap)) {
    cMpositions <- as.numeric(unlist(lapply(strsplit(dimnames(IBDarray)[[1]],"cM"),function(x) x[2])))
  } else{
    cMpositions <- IBD_list[[linkage_group]]$map$position
  }
  
  LG.data <- LOD_data$QTL.res[LOD_data$QTL.res$chromosome == linkage_group,]
  
  if(is.null(cM)){
    cM <- LG.data[which.max(LG.data$LOD),"position"]
    
    if(!is.null(LOD_data$Perm.res)){
      if(LG.data[which.max(LG.data$LOD),"LOD"] < LOD_data$Perm.res$threshold)
        stop(paste("All positions on linkage group", linkage_group,"fall below the significance threhold of",
                   round(LOD_data$Perm.res$threshold,2)))
    }
  } else{
    if(!cM %in% cMpositions) {
      message("Allowed positions:")
      print(cMpositions)
      warning(paste("Postion",cM,"not exactly matched. Proceeding with position",cMpositions[which.min(abs(cMpositions - cM))],"instead"))
      cM <- cMpositions[which.min(abs(cMpositions - cM))]
    }
    
    if(!is.null(LOD_data$Perm.res)){
      if(LG.data[LG.data$position == cM,"LOD"] < LOD_data$Perm.res$threshold)
        warning(paste0("Position ",cM," on linkage group ", linkage_group,
                       " falls below the significance threhold of ",
                       round(LOD_data$Perm.res$threshold,2), ". Proceeding regardless..."))
    }
  }
  
  write(paste("Fitting specified QTL models at position",cM,"cM on linkage group",linkage_group),
        log.conn)
  
  ## First get the starting set - the ones that were phenotyped and genotyped:
  
  if(any(is.na(Phenotype.df[,genotype.ID]))){
    warning("Missing genotype.ID values detected. Removing and proceeding without...")
    Phenotype.df <- Phenotype.df[!is.na(Phenotype.df[,genotype.ID]),]
  }
  
  phenoGeno <- intersect(
    dimnames(IBDarray)[[3]],
    unique(Phenotype.df[,genotype.ID]))
  
  write(paste("\nThere are",length(phenoGeno),"individuals with matching identifiers between phenotypic and genotypic datasets.\n"),log.conn)
  
  pheno <- Phenotype.df[Phenotype.df[,genotype.ID] %in% phenoGeno &
                          !is.na(Phenotype.df[,trait.ID]),c(genotype.ID,trait.ID)]
  
  pheno <- pheno[order(pheno[,1]),]
  phenoGeno <- unique(pheno[,1])
  popSize <- length(phenoGeno)
  
  if(popSize != nrow(pheno)) stop("Check input - unequal dimensions of phenotypes and genotype data.\nSuggest to use polyqtlR:::BLUE() to generate consensus phenotypes first.")
  
  Phenotypes <- pheno[,2] #- mean(pheno[,2], na.rm=TRUE) #mean-centred phenotypes - is this really needed?
  
  peakProbs <- IBDarray[which.min(abs(cMpositions - cM)),,match(phenoGeno,dimnames(IBDarray)[[3]])]
  
  if(IBDtype == "genotypeIBD"){
    ALL.phenoSums <- peakProbs %*% Phenotypes
    ALL.genoSums <- rowSums(peakProbs)
    ALL.means <- ALL.phenoSums/ALL.genoSums
    ALL.means[is.nan(ALL.means[,1]),] <- NA #replace NaN with NA (as ALL.means is also in final output)
    nterms <- length(which(!is.na(ALL.means))) #use this instead of Nstates, though shouldn't make a big difference...
    
  } else{
    stop("BIC test has not yet been implemented for haplotype IBDs!")
  }
  
  ## Additive:
  if(length(add_vect) > 0){
    Bic.add <- sapply(QTLseg[add_vect], function(Qconfig){
      counts <- apply(GenotypeCodes, 1,function(x) length(which(x %in% Qconfig)))
      terms <- lapply(unique(counts), function(x) which(counts == x))
      nums <- sapply(terms,function(x) length(intersect(x,which(!is.na(ALL.means)))))
      
      aov.weights <- lm(ALL.means ~ counts, weights = ALL.genoSums, na.action = na.exclude)
      RSS1 <- anova(aov.weights)[2,2]
      
      return(round(2*log(nterms) + nterms*log(RSS1/nterms),4))
    })
  } else{
    Bic.add <- NULL
  }
  
  ## Dominant:
  if(length(dom_vect) > 0){
    Bic.dom <- sapply(QTLseg[dom_vect], function(Qconfig){
      # q.group <- which(apply(GenotypeCodes,1,function(x) length(intersect(x,Qconfig))) == 0) 
      # Q.group <- setdiff(1:Nstates,q.group)
      
      # Q.terms <- ALL.means[Q.group]
      # q.terms <- ALL.means[q.group]
      
      ## Redefine n1 and n2 based on the number of actual non-missing terms, as with empty IBD classes there is division by 0 in calculation of means.
      # n1 <- length(Q.terms[!is.na(Q.terms)])
      # n2 <- length(q.terms[!is.na(q.terms)])
      
      # sum.probs <- c(sum(ALL.genoSums[Q.group]), sum(ALL.genoSums[q.group]))
      
      # varQ.group <- var(Q.terms, na.rm = TRUE)
      
      # if(n2 > 1) {
      # varq.group <- var(q.terms, na.rm = TRUE)
      # logL <- (-Nstates/2)*(log(2*pi) + 1 + log((1/Nstates)*(varQ.group*(n1-1) + varq.group*(n2-1))))
      # 
      # } else{ #n2 = 1 so varq.group is undefined
      # logL <- (-Nstates/2)*(log(2*pi) + 1 + log((1/Nstates)*(varQ.group*(n1-1))))
      # }
      
      counts <- apply(GenotypeCodes, 1,function(x) length(which(x %in% Qconfig)))
      
      simplexcounts <- counts #only test for simplex dominance at the moment..
      
      simplexcounts[counts > 0] <- 1 #re-code as simplex dominant
      # duplexcounts[counts < 2] <- 0
      # duplexcounts[counts >= 2] <- 1 #re-code as duplex dominant
      
      aov.weights <- lm(ALL.means ~ simplexcounts, weights = ALL.genoSums, na.action = na.omit)
      RSS1 <- anova(aov.weights)[2,2]
      
      return(round(2*log(nterms) + nterms*log(RSS1/nterms),4))
    })
  } else{
    Bic.dom <- NULL
  }
  
  BIC <- c(Bic.add,Bic.dom)[order(c(add_vect,dom_vect))]
  
  plausible.configs <- order(BIC - min(BIC) - deltaBIC)[sort(BIC - min(BIC)) <= deltaBIC]
  
  if(plotBIC){
    plot(BIC,ylab="BIC",xlab="QTLconfig",pch="")
    text(1:length(BIC),BIC,1:length(BIC))
    
    # highlight BIC within deltaBIC of the minimum in green:
    text(plausible.configs,BIC[plausible.configs],plausible.configs,col="limegreen")
    
    # highlight min in red:
    text(which.min(BIC),min(BIC),which.min(BIC),col="red2",font=2)
    abline(h = min(BIC) + deltaBIC, col = "limegreen", lty = 3, lwd = 2)
  }
  
  ## Incorporate testAllele_Effects function here:
  if(testAllele_Effects){
    
    # Phenotypes <- pheno[,2] #Return to non mean-centred phenotypes
    outList <- lapply(seq(length(plausible.configs)), function(n) NULL) #Initialise
    counter <- 1
    
    if(IBDtype == "genotypeIBD"){
      haploProbs <- t(peakProbs) %*% indicatorMatrix
    } else{
      haploProbs <- t(peakProbs)
    }
    
    for(plaus.bic in plausible.configs){
      
      BICconfig <- QTLconfig[[plaus.bic]]
      temp.config <- unlist(BICconfig$homs) #This is the actual configuration being tested
      
      ## Begin with Additive case
      if(BICconfig$mode == "a"){
        
        DRtemp <- which(haploProbs[,temp.config] > 1)
        
        if(length(DRtemp) > 0){
          
          temp.haploProbs <- haploProbs
          temp.haploProbs[which(temp.haploProbs > 1)] <- 1
          
          if(length(temp.config) == 1){
            weights <- haploProbs[,temp.config] #If a single QTL allele, DR products are ok to use
          } else{
            weights <- apply(temp.haploProbs[,temp.config, drop = FALSE],1,prod) #If there are multi QTL alleles, cannot allow DR probabilities anymore
            ## as they would distort the product as well as the probability of absence (1-x operation below). So use replacement with 1 instead.
          }
          min.weights <- apply(apply(temp.haploProbs[,temp.config, drop = FALSE],2,function(x) 1-x),1,prod)
          
        } else{ #there are no double reduction probabilities...
          weights <- apply(haploProbs[,temp.config, drop = FALSE],1,prod)
          min.weights <- apply(apply(haploProbs[,temp.config, drop = FALSE],2,function(x) 1-x),1,prod)
        }
        
        ## Renormalise weights to sum to popSize:
        norm <- popSize / sum(c(weights,min.weights))
        min.weights <- min.weights*norm
        weights <- weights*norm
        
        plus.mean <- sum(weights*Phenotypes)/sum(weights)
        plus.sd <-  sqrt(weighted.var(Phenotypes,weights))
        
        minus.mean <- sum(min.weights*Phenotypes)/sum(min.weights)
        minus.sd <- sqrt(weighted.var(Phenotypes,min.weights))
        
      } else{ #Dominant case
        ## For the dominant case we have to handle things differently to generate a weighted mean:
        q.group <- unlist(sapply(1:Nstates, function(n) if(length(setdiff(GenotypeCodes[n,],temp.config)) == ploidyF1) n))
        Q.group <- setdiff(1:Nstates,q.group)
        
        n1 <- length(Q.group)
        n2 <- length(q.group)
        
        ## Need these for the output - to estimate the effective number of individuals in the 2 groups
        weights <- ALL.genoSums[Q.group]
        min.weights <- ALL.genoSums[q.group]
        
        if(min(n1,n2) == 0) stop("No phenotypic segregation expected for given QTL configuration.")
        if(n1+n2 != Nstates) stop("Unexpected number of groups encountered")
        
        Q.terms <- ALL.means[Q.group]
        q.terms <- ALL.means[q.group]
        
        plus.mean <- mean(Q.terms, na.rm = TRUE) #+  mean(pheno[,2], na.rm=TRUE) #Add back the overall mean to ALL.means terms (mean-centred)
        plus.sd <- sd(Q.terms,na.rm = TRUE)
        minus.mean <-  mean(q.terms, na.rm = TRUE) #+  mean(pheno[,2], na.rm=TRUE) #Add back the overall mean to ALL.means terms (mean-centred)
        minus.sd <- sd(q.terms, na.rm = TRUE)
      }
      
      ##Simplify the notation for the moment:
      noQTL <- "o"
      yesQTL <- "Q"
      
      seg.msg <- rep(noQTL, ploidy + ploidy2)
      seg.msg[unlist(BICconfig$homs)] <- yesQTL
      seg.msg <- paste(
        paste0(seg.msg[1:ploidy],collapse=""),
        paste0(seg.msg[(ploidy+1):(ploidy+ploidy2)],collapse=""),
        sep = " x ")
      mode.msg <- ifelse(BICconfig$mode == "a","additive","dominant")
      
      write(paste("\nThe QTL configuration for trait",trait.ID,
                  "at", cM,"cM on linkage group",linkage_group,
                  "with seg type:\n",seg.msg,"and",mode.msg,"gene action had BIC",BIC[plaus.bic],"\n"),log.conn)
      
      temp.out <- data.frame("mean" = round(c(plus.mean,minus.mean),3),
                             "sd" = round(c(plus.sd,minus.sd),3),
                             "N" = round(c(sum(weights),sum(min.weights))))
      
      rownames(temp.out) <- c("+","-")
      
      outList[[counter]] <- list("BIC" = BIC[plaus.bic],
                                 "Segtype" = seg.msg,
                                 "Action" = mode.msg,
                                 "Allele effects" = temp.out)
      names(outList)[counter] <- paste0(paste0(temp.config,collapse = "_"),substr(mode.msg,1,1))
      counter <- counter + 1
    }
    
  } else{
    outList = NULL
  }
  
  if (!is.null(log))
    close(log.conn)
  
  return(invisible(list(linkage_group = linkage_group,
                        cM = cM,
                        BIC = BIC,
                        Allele.effects = outList,
                        genotype.means = ALL.means)))
} #exploreQTL

#' Function to find the position of maximum LOD on a particular linkage group
#' @description Given QTL output, this function returns the position of maximum LOD for a specified linkage group.
#' @param LOD_data Output of \code{\link{QTLscan}} function.
#' @param linkage_group Numeric identifier of the linkage group being tested, based on the order of \code{IBD_list}.
#' Only a single linkage group is allowed.
#' @param verbose Should messages be written to standard output? By default \code{TRUE}.
#' @examples 
#' data("qtl_LODs.4x")
#' findPeak(LOD_data=qtl_LODs.4x,linkage_group=1)
#' @export
findPeak <- function(LOD_data,
                     linkage_group,
                     verbose = TRUE){
  if(is.null(LOD_data$Perm.res)){
    warning("No significance threshold available to check")
    thresh <- 0
  } else{
    thresh <- LOD_data$Perm.res$threshold
  }
  
  lgdata <- LOD_data$QTL.res[LOD_data$QTL.res$chromosome == linkage_group,]
  maxdata <- lgdata[which.max(lgdata$LOD),]
  
  if(maxdata$LOD < thresh) warning(paste("Peak on LG",linkage_group,
                                         "does not exceed the significance threshold"))
  
  if(verbose) print(cbind(maxdata,thresh = thresh))
  
  return(maxdata$position)
} #findPeak


#' Function to find a LOD - x support interval around a QTL position
#' @description Given QTL output, this function returns the LOD - x support for a specified linkage group, taking the 
#' maximum LOD position as the desired QTL peak.
#' @param LOD_data Output of \code{\link{QTLscan}} function.
#' @param linkage_group Numeric identifier of the linkage group being tested, based on the order of \code{IBD_list}.
#' Only a single linkage group is allowed.
#' @param LOD_support The level of support around a QTL peak, by default 2 (giving a LOD - 2 support interval, the 
#' range of positions with a LOD score within 2 LOD units of the maximum LOD on that linkage group). 
#' @examples 
#' data("qtl_LODs.4x")
#' findSupport(LOD_data=qtl_LODs.4x,linkage_group=1)
#' @export
findSupport <- function(LOD_data,
                        linkage_group,
                        LOD_support = 2){
  
  if(!LOD_support > 0) stop("LOD_support should be a positive number!")
  
  lgdata <- LOD_data$QTL.res[LOD_data$QTL.res$chromosome == linkage_group,]
  LODrange <- lgdata$LOD >= max(lgdata$LOD) - LOD_support
  
  return(range(lgdata$position[LODrange]))
} #findSupport


#' Import IBD probabilities as estimated by TetraOrigin or PolyOrigin
#' @description Imports the IBD probability output of TetraOrigin (Mathematica software) or PolyOrigin (julia software) into the same
#' format as natively-estimated IBD probabilities from the polyqtlR package.  
#' @param method The method used for IBD estimation, either "TO" for TetraOrigin or "PO" for PolyOrigin
#' @param folder The path to the folder in which the Tetra/PolyOrigin (ie. TetraOrigin or PolyOrigin) output is contained, 
#' default is \code{NULL} if files are in working directory.
#' @param filename If method = "TO", the (vector of) character filename stem(s) of the \code{.csv} file(s) containing the output of TetraOrigin (stem = without ".csv"). 
#' Should be in order according to LG/chromosome numbering. If method = "PO", then simply specify the PolyOrigin filename stem here (as the output is not split into separate linkage
#' groups in PolyOrigin). A PolyOrigin file with name <filename>_polyancestry.csv and its corresponding log file <filename>.log will then be searched for.
#' @param bivalent_decoding Logical, if method = "TO" you must specify \code{TRUE} if only bivalent pairing was allowed in TetraOrigin (in offspring deciding step),
#'  otherwise specify \code{FALSE} if multivalent pairing was also allowed. If method = "PO", this will be automatically detected, so no need to specify (will be ignored).
#' @param error If method = "TO", the offspring error prior used in the offspring decoding step of TetraOrigin, by default assumed to be 0.01. For method = "PO", this is automatically read in.
#' @param log Character string specifying the log filename to which standard output should be written. If \code{NULL} log is send to stdout.
#' @return Returns a list with the following items:
#' \item{IBDtype : }{Always "genotypeIBD" for the output of TetraOrigin}
#' \item{IBDarray : }{An array of IBD probabilities. The dimensions of the array are: markers, genotype classes and individuals.}
#' \item{map : }{Integrated linkage map positions of markers used in IBD calculation}
#' \item{parental_phase : }{The parental marker phasing as used by TetraOrigin, recoded in 1 and 0's}
#' \item{marginal.likelihoods : }{A list of marginal likelihoods of different valencies, currently \code{NULL}}
#' \item{valency : }{The predicted valency that maximised the marginal likelihood, per offspring. Currently \code{NULL}}
#' \item{offspring : }{Offspring names}
#' \item{biv_dec : }{Logical, the bivalent_decoding parameter specified.}
#' \item{gap : }{The gap size used in IBD interpolation if performed by \code{\link{spline_IBD}}. At this stage, \code{NULL}}
#' \item{genocodes : }{Ordered list of genotype codes used to represent different genotype classes.}
#' \item{pairing : }{log likelihoods of each of the different pairing scenarios considered 
#' (can be used e.g. for post-mapping check of preferential pairing)}
#' \item{ploidy : }{The ploidy of parent 1, by default assumed to be 4}
#' \item{ploidy2 : }{The ploidy of parent 2, by default assumed to be 4}
#' \item{method : }{The method used, either "hmm_TO" (TetraOrigin) or "hmm_PO" (PolyOrigin)}
#' \item{error : }{The error prior used in the calculation in TetraOrigin, assumed to be 0.01}
#' @examples
#' \dontrun{
#' ## These examples demonstrate the function call for both methods, but won't run without input files
#' ## from either package, hence this call will normally result in an Error:
#' IBD_TO <- import_IBD(method = "TO", filename = paste0("test_LinkageGroup",1:5,"_Summary"),
#' bivalent_decoding = FALSE, error = 0.05)
#' ## Equivalent call for PolyOrigin output:
#' IBD_PO <- import_IBD(method = "PO",filename = "test")
#' }
#' @export
import_IBD <- function(method,
                       folder=NULL,
                       filename,
                       bivalent_decoding = TRUE,
                       error = 0.01,
                       log=NULL){
  
  if(is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  method <- match.arg(method, choices = c("TO","PO"))
  
  if(method == "TO"){
    outlist <- lapply(filename, function(fname){
      if(!is.null(folder)) {
        path.to.file <- file.path(folder,paste0(fname,".csv"))
      } else {
        path.to.file <- paste0(fname,".csv")
      }
      
      suppressWarnings(input_lines <- readLines(path.to.file))
      sink(file = log.conn, append=TRUE)
      placeholders <- grep("inferTetraOrigin-Summary",input_lines)
      
      write(paste("Importing map data under description",input_lines[placeholders[1]]),log.conn)
      mapData <- read.csv(path.to.file, skip = 1,
                          nrows = (placeholders[2] - placeholders[1] - 2),
                          stringsAsFactor=FALSE,
                          header = TRUE)
      colnames(mapData) <- c("marker","chromosome","position") #make compatible with polyqtlR
      
      
      write(paste("Importing parental phasing under description",input_lines[placeholders[2]]),log.conn)
      parentalPhase <- read.csv(path.to.file,
                                skip = placeholders[2],
                                nrows = (placeholders[3] - placeholders[2] - 2),
                                stringsAsFactor=FALSE,header = TRUE)
      rownames(parentalPhase) <- parentalPhase[,1]
      parentalPhase <- t(parentalPhase[-1,-1])
      parentalPhase[parentalPhase==1] <- 0
      parentalPhase[parentalPhase==2] <- 1
      
      write(paste("Importing pairing likelihoods under description",input_lines[placeholders[4]]),log.conn)
      ln_valent <- read.csv(path.to.file,
                            skip = placeholders[4],
                            nrows = (placeholders[5] - placeholders[4] - 2),
                            stringsAsFactor=FALSE,header = TRUE)
      newcolnames <- sapply(sapply(ln_valent[1,2:ncol(ln_valent)], strsplit,"-"),
                            function(x) paste0(sapply(x, function(y) paste0(LETTERS[as.numeric(unlist(strsplit(y,"")))], collapse="")),
                                               collapse="-")
      )
      ln_valent <- ln_valent[-1,]
      
      ## Transform to a list to conform to correct format:
      ln_valent.ls <- setNames(lapply(1:nrow(ln_valent), function(r) ln_valent[r,2:ncol(ln_valent)]),ln_valent[,1])
      
      ## Express these as relative likelihoods, rather than log likelihoods..
      ln_valent.ls <- lapply(ln_valent.ls, function(x) {
        out <- exp(as.numeric(x) - max(as.numeric(x)))
        ## Re-normalise:
        out <- out / sum(out)
        return(setNames(out, newcolnames))
      })
      
      
      write(paste("Importing genocodes under description",input_lines[placeholders[6]]),log.conn)
      genocodes <- read.csv(path.to.file,
                            skip = placeholders[6],
                            nrows = (placeholders[7] - placeholders[6] - 2),
                            stringsAsFactor=FALSE,header = TRUE)
      genocodes <- genocodes$Code
      ## polyqtlR uses letter format for the genotypes, convert:
      genocodes <- sapply(genocodes, function(x=genocodes[1]) paste0(letters[as.numeric(strsplit(as.character(x),"")[[1]])],collapse = ""))
      
      write(paste("Importing IBD data under description",input_lines[placeholders[7]]),log.conn)
      
      if(length(placeholders) != 8){
        if(length(placeholders) != 7) stop("Insufficient data to proceed. Please check TetraOrigin summary file and/or re-run TetraOrigin.")
        
        warning("Not all expected sections encountered in TetraOrigin summary file. Attempting to proceed with current data file...")
        genotypeProbs <- read.csv(path.to.file,
                                  skip = placeholders[7],
                                  stringsAsFactor=FALSE,header = TRUE)
      } else{
        genotypeProbs <- read.csv(path.to.file,
                                  skip = placeholders[7],
                                  nrows = (placeholders[8] - placeholders[7] - 2),
                                  stringsAsFactor=FALSE,header = TRUE)
      }
      
      outRowNames <- genotypeProbs[2:nrow(genotypeProbs),1]
      F1names <- unique(sapply(outRowNames, function(x) strsplit(x,"_genotype")[[1]][1]))
      genotypeProbs <- as.matrix(genotypeProbs[-1,-1])
      genotypeProbs <- apply(genotypeProbs,2,as.numeric)
      
      Nstates <- Nstates.fun(biv_dec = bivalent_decoding, pl = 4, pl2 = 4)
      
      if(nrow(genotypeProbs) %% Nstates != 0) stop("Incompatible input detected. Make sure bivalent_decoding is correctly specified.")
      
      F1size <- nrow(genotypeProbs)/Nstates
      
      IBDoutput <- genotypeProbs #No splining here
      
      # colnames(IBDoutput) <- paste0("cM",mapData[,3]) #why did I do this? For compatibility with splined IBDs? Removed 16.06.23
      
      # rownames(IBDoutput) <- outRowNames
      
      if(!is.null(log)) {
        sink()
        close(log.conn)
      }
      
      ## Turn the IBDoutput into a 3d array (markers x genotype_class x F1_individual)
      # IBDarray <- array(t(IBDoutput),dim = c(ncol(IBDoutput),Nstates,F1size),
      #                   dimnames = list(colnames(IBDoutput),
      #                                   paste0("g",1:Nstates),
      #                                   F1names))
      
      
      IBDarray <- array(t(IBDoutput),dim = c(ncol(IBDoutput),Nstates,F1size),
                        dimnames = list(colnames(IBDoutput),
                                        genocodes, #updated 16.06.23
                                        F1names))
      
      return(list(IBDtype = "genotypeIBD",
                  IBDarray = IBDarray,
                  map = mapData,
                  parental_phase = parentalPhase,
                  marginal.likelihoods = NULL,
                  valency = NULL,
                  offspring = F1names,
                  biv_dec = bivalent_decoding,
                  gap = NULL,
                  genocodes = genocodes,
                  pairing = ln_valent.ls,
                  ploidy = 4,
                  ploidy2 = 4,
                  method = "hmm_TO",
                  error = error))
    })
    
    names(outlist) <- paste0("LG", seq(length(filename)))
    
  } else{
    ## We are dealing with PolyOrigin output
    
    ## In this case, filename is a single filename stem, not a vector of filenames (split up per linkage group) 
    if(!is.null(folder)) {
      path.to.file <- file.path(folder,paste0(filename,"_polyancestry.csv"))
      path.to.file2 <- file.path(folder,paste0(filename,".log"))
    } else {
      path.to.file <- paste0(filename,"_polyancestry.csv")
      path.to.file2 <- paste0(filename,".log")
    }
    
    suppressWarnings(input_lines <- readLines(path.to.file))
    # sink(file = log.conn, append=TRUE)
    placeholders <- grep("PolyOrigin-PolyAncestry",input_lines)
    
    ## Import the parental phase and map positions (parentgeno)
    parentalPhase <- read.csv(path.to.file,
                              skip = placeholders[8],
                              nrows = (placeholders[9] - placeholders[8] - 2),
                              stringsAsFactor=FALSE,header = TRUE)
    
    p1phase <- parentalPhase$P1
    p2phase <- parentalPhase$P2
    
    P1 <- do.call(rbind,strsplit(p1phase,"\\|"))
    P1 <- sapply(1:ncol(P1), function(x) as.numeric(P1[,x])) - 1 #code in 0 and 1 instead of 1 and 2
    ploidy <- ncol(P1)
    
    P2 <- do.call(rbind,strsplit(p2phase,"\\|"))
    P2 <- sapply(1:ncol(P2), function(x) as.numeric(P2[,x])) - 1
    ploidy2 <- ncol(P2)
    
    mapData <- parentalPhase[,c("chromosome","marker","position")]
    
    parphase <- cbind(P1,P2)
    colnames(parphase) <- paste0("h",1:ncol(parphase))
    rownames(parphase) <- parentalPhase[,"marker"]
    
    ## Extract offspring names
    offspringinfo <- read.csv(path.to.file,
                              skip = placeholders[3],
                              nrows = (placeholders[4] - placeholders[3] - 2),
                              stringsAsFactor=FALSE,header = TRUE)
    pop <- offspringinfo$individual
    
    ## Import pairing likelihoods (valentprob)
    valentprob <- read.csv(path.to.file,
                           skip = placeholders[7],
                           nrows = (placeholders[8] - placeholders[7] - 2),
                           stringsAsFactor=FALSE,header = TRUE)
    
    if(ploidy == 4){
      valent.loglike <- lapply(sapply(valentprob[,"loglike"],strsplit,"\\|"),function(x) as.numeric(x))
      valent.prob <- lapply(sapply(valentprob[,"valentprob"],strsplit,"\\|"),function(x) as.numeric(x))
      valent.index <- lapply(sapply(valentprob[,"valentindex"],strsplit,"\\|"),function(x) as.numeric(x))
    } else if(ploidy == 2){
      valent.loglike <- valentprob[,"loglike"]
      valent.prob <- valentprob[,"valentprob"]
      valent.index <- valentprob[,"valentindex"]
    } else{
      stop("Function not yet checked on data other than 2x and 4x. Please contact the developer!")
    }
    
    suppressWarnings(logfile_lines <- readLines(path.to.file2))
    
    ## To generate $bivalent_decoding field, need to check the log file (note: by default chrpairing == 44)
    if(ploidy == 4){
      
      bd <- grep(pattern = "chrpairing =", logfile_lines)
      
      if(length(bd) == 0){
        warning("Could not match chrpairing in PolyOrigin logfile. Assuming multivalent pairing was allowed in offspring phasing!")
        bivdec <- FALSE #Assume chrpairing = 44
      } else{
        cp <- strsplit(logfile_lines[bd][1],"chrpairing = ")[[1]][2]
        if(cp == "44"){
          bivdec <- FALSE
        } else{
          bivdec <- TRUE 
        }
      }
      
    } else if(ploidy == 2){
      bivdec <- TRUE #by definition a diploid pairs in bivalents, so chrpairing == 44 has no meaning
    } else{
      stop("Function not yet checked on data other than 2x and 4x. Please contact the developer!")
    }
    
    
    ## To generate $error field, also need to check the log file 
    polyreconstruct <- grep("PolyOrigin, polyReconstruct", logfile_lines)
    er <- grep(pattern = "epsilon = ", logfile_lines)
    er.hit <-er[er > max(polyreconstruct)]
    
    if(length(er.hit) == 0){
      warning("Could not match epsilon in PolyOrigin logfile. Assuming epsilon = 0.01")
      eps <- 0.01
    } else{
      eps <- as.numeric(strsplit(logfile_lines[er.hit],"epsilon = ")[[1]][2])
    }
    
    ## Import the coding for the valenctindex:
    valentlist <- read.csv(path.to.file,
                           skip = placeholders[6],
                           nrows = (placeholders[7] - placeholders[6] - 2),
                           stringsAsFactor=FALSE,header = TRUE)
    if(ploidy == 4){
      
      P1Q <- grep(pattern = paste0(1:ploidy,collapse=":"),x = valentlist$valent)
      P2Q <- grep(pattern = paste0((ploidy+1):(ploidy+ploidy2),collapse=":"),x = valentlist$valent)
      
      QQ <- intersect(P1Q, P2Q)
      QB <- setdiff(P1Q,QQ)
      BQ <- setdiff(P2Q,QQ)  
      BB <- setdiff(1:nrow(valentlist),c(QQ,QB,BQ))
      
      vl <- 1:nrow(valentlist)
      for(r in 1:length(vl)){
        if(vl[r] %in% BB) names(vl)[r] <- "BB"
        if(vl[r] %in% QB) names(vl)[r] <- "QB"
        if(vl[r] %in% BQ) names(vl)[r] <- "BQ"
        if(vl[r] %in% QQ) names(vl)[r] <- "QQ"
      }
      
      ## This will be the $valency output (workaround to maintain compatibility)
      maxL.valent <- sapply(1:length(valent.prob), function(n=1){
        return(names(vl[valent.index[[n]][which.max(valent.prob[[n]])]]))
      } )
      
      
    } else if(ploidy == 2){
      # vl <- "BB" #Assumes length of valentlist is 1, which is true for a diploid
      maxL.valent <- rep("BB",length(valent.prob))
      
    } else{
      ## This Error would already have been given, but just for completeness:
      stop("Function not yet checked on data other than 2x and 4x. Please contact the developer!")
    }
    
    ## I will not try to replicate the output $marginal likelihoods for now
    ## as these are just logL of BB, BQ, QB and QQ, while PO gives logL of 16 valencies.
    ## In time, update estimate_IBD to be compatible (I guess it unlikely that users will continue
    ## to use TetraOrigin given the availability of PolyOrigin now)
    
    ## To generate $pairing, need to transform the output of PO into polyqtlR format.
    ## Convert the names of valentlist into polyqtlR coding
    valent.recoding <- sapply(valentlist$valent, function(x = valentlist$valent[16]){
      psplit <- strsplit(x,"\\|")[[1]]
      bivs <- grep("-",psplit) #the dash divides bivalents (in 4x, still need to check this works in 2x)
      bivpairs <- quadpairs <- NULL
      
      if(length(bivs) > 0){
        split2 <- do.call(c,lapply(sapply(psplit[bivs], strsplit,"\\-"),strsplit,"\\:"))
        bivpairs <- sapply(split2,function(x) paste0(LETTERS[as.numeric(x)], collapse = ""))
        if(length(bivs) != ploidy/2){ #there are quads
          quadpairs <- sapply(sapply(psplit[setdiff(1:(ploidy/2),bivs)], strsplit,"\\:"), function(x) paste0(LETTERS[as.numeric(x)], collapse = ""))
        }
      } else{
        quadpairs <- sapply(sapply(psplit, strsplit,"\\:"), function(x) paste0(LETTERS[as.numeric(x)], collapse = ""))
      }
      return(paste0(sort(c(bivpairs,quadpairs)),collapse="-"))
    })
    
    valent.loglike.recoded <- sapply(1:length(valent.loglike), function(n=1){
      setNames(exp(valent.loglike[[n]])/sum(exp(valent.loglike[[n]])),valent.recoding[valent.index[[n]]])
    })
    
    ## In estimate_IBD we took the approach of classifying an offspring as coming from bivalent or 
    ## multivalent pairing, as proposed in TetraOrigin. The approach of PolyOrigin appears to be
    ## more flexible in this regard - it returns the likelihoods of all valencies with non-zero probability
    ## This could cause potentially unwanted behaviour in some polyqtlR functions e.g. estimate_GIC. 
    
    ## For now, I will classify the PolyOrigin output as Q if the most likely configuration includes a 
    ## quadrivalent, and B otherwise. May need to revisit this, or update estimate_IBD to make compatible 
    ## with the PolyOrigin approach (which is broadly similar)
    
    ## Also import the ancestralgenotype:
    ancestralgenotype <- read.csv(path.to.file,
                                  skip = placeholders[9],
                                  nrows = (placeholders[10] - placeholders[9] - 2),
                                  stringsAsFactor=FALSE,header = TRUE)
    
    ancestrycodes <- do.call(rbind,sapply(ancestralgenotype$state, function(x) strsplit(x,"-")))
    genocodes <- apply(sapply(1:((ploidy+ploidy2)/2), function(cl) letters[as.numeric(ancestrycodes[,cl])]),1,paste0,collapse="")
    
    LG.names <- LGs <- unique(parentalPhase$chromosome) 
    # LGs = c("1","2")
    # LG.names will be used as the output list names also. Check that they are not purely numeric; if so, paste LG onto them:
    numeric.test <- suppressWarnings(as.numeric(LG.names))
    
    if(any(!is.na(numeric.test))) {
      warning(paste("Purely numeric LG names detected for:",paste0(LG.names[!is.na(numeric.test)],collapse=", "),"\nPasting LG as suffix for these."))
      LG.names[!is.na(numeric.test)] <- paste0("LG",LG.names[!is.na(numeric.test)])
    }
    
    
    ## Finally, extract the IBD array:
    genoprobs <- read.csv(path.to.file,
                          skip = placeholders[10],
                          stringsAsFactor=FALSE,header = TRUE)
    
    IBD.ls <- setNames(lapply(LGs, function(lg = "LG12"){
      
      lghits <- which(parentalPhase$chromosome == lg)  #these are the rows for this LG
      
      genoProbs <- array(0, dim=c(length(lghits), length(genocodes), length(pop)),
                         dimnames = list(parentalPhase[lghits,"marker"],
                                         genocodes,pop))
      
      counter <- 1
      
      for(r in lghits){
        index.probs <- sapply(genoprobs[r,pop],strsplit,"=>")
        ## Need to assign names from 1st elements to second, and handle in a list:
        named.probs <- lapply(index.probs, function(x=index.probs[[4]]) {
          y = sapply(x,strsplit,"\\|")
          return(setNames(as.numeric(y[[2]]), genocodes[as.numeric(y[[1]])]))
        })
        
        for(i in 1:length(named.probs)){
          genoProbs[counter,names(named.probs[[i]]),pop[i]] <- named.probs[[i]]
        }
        counter <- counter + 1
      }
      
      return(genoProbs)
    }),LG.names) #this could be causing a bug..
    
    Npop <- length(pop)
    
    ## Put it all together, per linkage group:
    counter <- 1
    outlist <- setNames(lapply(1:length(LGs), function(i){
      
      lghits <- which(parentalPhase$chromosome == LGs[i])  #these are the rows for this LG
      
      # F1hits <- ((lg-1)*Npop + 1):(Npop*lg) #bug here if lg is non-numeric. Use counter instead
      F1hits <- ((counter-1)*Npop + 1):(Npop*counter)
      counter <<- counter + 1 #workaround to increment counter..
      
      return(list(IBDtype = "genotypeIBD",
                  IBDarray = IBD.ls[[i]],
                  map = mapData[lghits,],
                  parental_phase = parphase[lghits,],
                  marginal.likelihoods = NULL,
                  valency = setNames(maxL.valent[F1hits],pop),
                  offspring = pop,
                  biv_dec = bivdec,
                  gap = NULL,
                  genocodes = genocodes,
                  pairing = setNames(valent.loglike.recoded[F1hits],pop),
                  ploidy = ploidy,
                  ploidy2 = ploidy2,
                  method = "hmm_PO",
                  error = eps))
      
    }), LG.names)
  }
  
  return(outlist)
} #import_IBD

#' Re-estimate marker dosages given IBD input estimated using a high error prior.
#' @description Function to correct marker dosage scores given a list of previously estimated IBD probabilities. This may
#' prove useful to correct genotyping errors. Running the \code{\link{estimate_IBD}} function with a high error prior will 
#' result in suppressed predictions of double recombination events, associated with genotyping errors. By forcing the HMM to penalise 
#' double recombinations heavily, a smoothed haplotype landscape is achieved in which individual genotype observations are down-weighted. 
#' This smoothed output is then used to re-estimate marker dosages, dependent on (correct) parental scores. 
#' An alternative strategy is to use the function \code{\link{maxL_IBD}} over a range of error priors first, and use the resulting \code{$maxL_IBD} output 
#' as input here (as the \code{IBD_list}). In this case, set the argument \code{min_error_prior} to a low value (0.005 say) to avoid issues. 
#' @param IBD_list List of IBD probabilities
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns. Note that probabilistic genotypes are not currently catered for here.
#' @param parent1 The identifier of parent 1, by default "P1"
#' @param parent2 The identifier of parent 2, by default "P2"
#' @param rounding_error The maximum deviation from an integer value that an inputed value can have, by default 0.05. For example, an imputed
#' score of 2.97 or 3.01 would both be rounded to a dosage of 3, while 2.87 would be deemed too far from an integer score, and would be made missing.
#' If you find the output contains too many missing values, a possibility would be to increase the \code{rounding_error}. However this may
#' also introduce more errors in the output! 
#' @param min_error_prior Suggestion for a suitably high error prior to be used in IBD calculations to ensure IBD smoothing is achieved. If IBD probabilities were estimated 
#' with a smaller error prior, the function aborts.
#' @param verbose Should messages be written to standard output?
#' @examples
#' \dontrun{
#' # Toy example only, as this will result in an Error: the original error prior was too low
#' data("IBD_4x","SNP_dosages.4x")
#' impute_dosages(IBD_list=IBD_4x,dosage_matrix=SNP_dosages.4x)
#' }
#' @export
impute_dosages <- function(IBD_list,
                           dosage_matrix,
                           parent1 = "P1",
                           parent2 = "P2",
                           rounding_error = 0.05,
                           min_error_prior = 0.1,
                           verbose = TRUE
){
  IBD_list <- test_IBD_list(IBD_list)
  
  ploidy <- IBD_list[[1]]$ploidy
  ploidy2 <- IBD_list[[1]]$ploidy2
  
  dosage_matrix <- test_dosage_matrix(dosage_matrix)
  
  if(!all(c(parent1,parent2) %in% colnames(dosage_matrix))) stop("Parental columns not found in dosage_matrix. Please check arguments parent1 and parent2!")
  
  if(IBD_list[[1]]$error < min_error_prior) stop(paste0("IBDs were estimated with an error prior of ", 
                                                        IBD_list[[1]]$error,
                                                        ", less than the set minimum of ",min_error_prior,
                                                        " provided. Either re-estimate IBDs using a higher prior, or reduce argument min_error_prior accordingly!"))
  
  IBDtype <- IBD_list[[1]]$IBDtype
  
  if(IBDtype == "genotypeIBD"){
    # if(IBD_list[[1]]$method == "hmm_TO"){
    #   ## We are using the output of TetraOrigin
    #   # Nstates <- Nstates.fun(biv_dec = bivalent_decoding, pl = ploidy, pl2 = ploidy2)
    #   # mname <- paste0("GenotypeMat",Nstates)
    #   # indicatorMatrix <- as.matrix(get(mname, envir = getNamespace("polyqtlR")))
    #   
    #   GenCodes <- t(sapply(IBD_list[[1]]$genocodes,function(x) as.numeric(strsplit(as.character(x),"")[[1]])))
    #   Nstates <- nrow(GenCodes)
    #   indicatorMatrix <- matrix(0,nrow=Nstates, ncol=ploidy + ploidy2)
    #   for(r in 1:nrow(indicatorMatrix)){
    #     for(h in 1:ncol(GenCodes)){
    #       indicatorMatrix[r,GenCodes[r,h]] <- indicatorMatrix[r,GenCodes[r,h]] + 1
    #     }
    #   }
    #   
    # } else{
    ## We are using the output of estimate_IBD
    GenCodes <- t(sapply(sapply(IBD_list[[1]]$genocodes,function(x) strsplit(x,"")), function(y) sapply(y,function(z) which(letters == z))))
    Nstates <- nrow(GenCodes)
    indicatorMatrix <- matrix(0,nrow=Nstates, ncol=ploidy + ploidy2)
    for(r in 1:nrow(indicatorMatrix)){
      for(h in 1:ncol(GenCodes)){
        indicatorMatrix[r,GenCodes[r,h]] <- indicatorMatrix[r,GenCodes[r,h]] + 1
      }
    }
    
    # }
    
  } else if(IBDtype == "haplotypeIBD"){
    Nstates <- ploidy + ploidy2 #this is not actually needed..
    indicatorMatrix <- NULL
  } else{
    stop("Unable to determine IBD type, please check IBD_list elements have a valid $IBDtype tag.")
  }
  
  ## Extract the names of the markers that can be imputed:
  ibd.markers <- do.call(c,lapply(seq(length(IBD_list)), function(lg) dimnames(IBD_list[[lg]]$IBDarray)[[1]]))
  
  common.markers <- intersect(rownames(dosage_matrix),ibd.markers)
  
  ## Check and correct marker compatibility between datasets if necessary:
  if(length(common.markers) == 0) stop("No common marker names between dosage_matrix and IBD_list. Perhaps a splined IBD list was provided? Please check input!")
  if(length(common.markers) < length(ibd.markers)) {
    warning(paste("It is only possible to identify dosages of",
                  length(common.markers),
                  "of the",length(ibd.markers),
                  "markers in IBD_list. Proceeding with these.."))
    for(lg in seq(length(IBD_list))){
      IBD_list[[lg]]$IBDarray <- IBD_list[[lg]]$IBDarray[intersect(dimnames(IBD_list[[lg]]$IBDarray)[[1]],
                                                                   common.markers),,]
      
      IBD_list[[lg]]$parental_phase <- IBD_list[[lg]]$parental_phase[intersect(rownames(IBD_list[[lg]]$parental_phase),
                                                                               common.markers),]
    }
  }
  
  ## Check and correct individual compatibility between datasets if necessary:
  common.inds <- intersect(setdiff(colnames(dosage_matrix),c(parent1,parent2)),dimnames(IBD_list[[1]]$IBDarray)[[3]])
  
  if(length(common.inds) == 0) stop("No common individual names between dosage_matrix and IBD_list. Please check both!")
  
  if(length(common.inds) < length(dimnames(IBD_list[[1]]$IBDarray)[[3]])) {
    warning(paste("It is only possible to identify dosages of",
                  length(common.inds),
                  "of the",length(dimnames(IBD_list[[1]]$IBDarray)[[3]]),
                  "individuals in IBD_list. Proceeding with these.."))
    for(lg in seq(length(IBD_list))){
      IBD_list[[lg]]$IBDarray <- IBD_list[[lg]]$IBDarray[,,intersect(dimnames(IBD_list[[lg]]$IBDarray)[[3]],
                                                                     common.inds)]
    }
  }
  
  imp <- function(lg=1){
    ibd.tmp <- IBD_list[[lg]]$IBDarray
    parphase.tmp <- IBD_list[[lg]]$parental_phase
    
    ## These can be combined to impute the most likely offspring scores at each marker position..
    dose.impute <- t(do.call(rbind, lapply(seq(dim(ibd.tmp)[3]), 
                                           function(ind) rowSums((ibd.tmp[,,ind] %*% indicatorMatrix) * parphase.tmp))))
    
    dose.impute[which(abs(round(dose.impute) - dose.impute) > rounding_error)] <- NA
    
    dose.impute <- round(dose.impute)
    
    dose.impute <- cbind(dosage_matrix[rownames(dose.impute),c(parent1, parent2)],
                         dose.impute)
    
    colnames(dose.impute)[3:ncol(dose.impute)] <- common.inds
    
    return(dose.impute)
  }
  
  new.dosages <- do.call(rbind, lapply(seq(length(IBD_list)), imp))
  
  ## Generate some comparison statistics
  if(verbose){
    old.dose <- dosage_matrix[rownames(new.dosages),]
    
    cat("\n____________________________________________________________________________\n\n")
    message(paste(nrow(new.dosages),"out of a total possible",nrow(dosage_matrix),"markers were imputed"))
    message(paste("In the original dataset there were",round(100*length(which(is.na(old.dose)))/length(old.dose),1),
                  "% missing values among the",nrow(new.dosages),"markers"))
    message(paste("In the imputed dataset there are now",round(100*length(which(is.na(new.dosages)))/length(new.dosages),1),
                  "% missing values among these, using a rounding threshold of",rounding_error))
    
    
    error.matrix <- abs(old.dose[,3:ncol(old.dose)] - new.dosages[,3:ncol(new.dosages)])
    cat("\n____________________________________________________________________________\n")
    message("The % of markers with changed dosage scores is as follows (0 = no change):")
    print(round(100*table(error.matrix)/length(error.matrix),2))
    
  }
  
  return(new.dosages)
} #impute_dosages

#' Wrapper function to run estimate_IBD function over multiple error priors
#' @description Function to run the \code{\link{estimate_IBD}} function over a range of possible error priors. The function returns
#' a merged set of results that maximise the marginal likelihood per individual, i.e. allowing a per-individual error rate within the options
#' provided in the errors argument.
#' @param errors Vector of offspring error priors to test (each between 0 and 1)
#' @param \dots Arguments passed to \code{\link{estimate_IBD}}. 
#' @return A list containing the following components:
#' \describe{
#' \item{maxL_IBD}{A nested list as would have been returned by the estimate_IBD function, but composite across error priors to maximise the 
#' marginal likelihoods. Note that the $error values per linkage group are now the average error prior across the population per linkage group}
#' \item{MML}{A 3d array of the maximal marginal likelihoods, per error prior. Dimensions are individuals, linkage groups, error priors.}
#' \item{error_per_ind}{A matrix of the most likely genotyping error rates per individual (in rows) for each linkage group (in columns)}
#' \item{errors}{The error priors used (i.e. the input vector is returned for later reference.)}
#' }
#' @examples
#' \dontrun{
#' data("phased_maplist.4x","SNP_dosages.4x")
#' maxL_IBD(phased_maplist=phased_maplist.4x,genotypes=SNP_dosages.4x,
#' ploidy=4,errors=c(0.01,0.02,0.05,0.1))
#' }
#' @export
maxL_IBD <- function(errors = c(0.01,0.05,0.1,0.2), 
                     ...){
  
  chkargs <- list(...)
  
  if("error" %in% names(chkargs)) stop("Cannot specify argument 'error' in this function, use
                                       argument 'errors' instead")
  if("method" %in% names(chkargs)){
    if(chkargs$method != "hmm") stop("This function only operates using method = 'hmm' (default)")
  }
  
  if(!all(sapply(errors,is.numeric))) stop("errors must be a vector of numeric values")
  if(length(errors) < 2) stop("At least two error priors are required for this function")
  if(any(errors <= 0) | any(errors >= 1)) stop("errors should be numeric in range (0,1)")
  
  IBD.list <- setNames(lapply(errors, function(er){
    estimate_IBD(error = er,...)
  }),paste0("e_",errors))
  
  # #debug:  
  # IBD.list <- setNames(lapply(errors, function(er){
  #   estimate_IBD(error = er,phased_maplist= phased_maplist.4x,
  #                genotypes=SNP_dosages.4x,
  #                ploidy=4)
  # }),paste0("e_",errors))
  
  # extract maximal marginal likelihoods per error prior:
  n.ind <- length(IBD.list[[1]][[1]]$offspring)
  mml.arr <- array(0,dim = c(n.ind,
                             length(IBD.list[[1]]),
                             length(errors)),
                   dimnames = list(IBD.list[[1]][[1]]$offspring,
                                   names(IBD.list[[1]]),
                                   paste0("e_",errors)
                   ))
  
  for(i in 1:length(errors)){
    mml.arr[,,i] <- sapply(IBD.list[[i]], function(x) sapply(x$marginal.likelihoods,max))
  }
  
  ## Generate a new IBD array composed of the values that maximise mml across error priors:
  IBDout <- setNames(lapply(1:length(IBD.list[[1]]), function(i){
    ## create a temporary template list, copied from the first error prior:
    temp <- IBD.list[[1]][[i]]
    
    ## get error priors that maximises max-mar-L for each individuals:
    max_mml <- apply(mml.arr[,i,],1,which.max)
    
    ## record the mean as approximate error prior for this LG:
    temp$error <- round(mean(errors[max_mml]),3)
    
    ## Record the marginal likelihoods of these:
    temp$marginal.likelihoods <- setNames(lapply(1:n.ind,function(n){
      IBD.list[[max_mml[n]]][[i]]$marginal.likelihoods[[n]]
    }),
    temp$offspring)
    
    ## Create an array with the correct entries, merging over max_mml as before:
    for(n in 1:n.ind){
      temp$IBDarray[,,n] <- IBD.list[[max_mml[n]]][[i]]$IBDarray[,,n]
    }
    
    temp$valency <- sapply(1:n.ind, function(n) IBD.list[[max_mml[n]]][[i]]$valency[n])
    
    temp$pairing <- setNames(lapply(1:n.ind, function(n) IBD.list[[max_mml[n]]][[i]]$pairing[[n]]),
                             temp$offspring)
    
    return(temp)
  }),
  names(IBD.list[[1]]))
  
  ## Extract the predicted errors per offspring, per LG:
  error_per_ind <- as.matrix(sapply(1:length(IBD.list[[1]]), function(i) errors[apply(mml.arr[,i,],1,which.max)]))
  colnames(error_per_ind) <- names(IBD.list[[1]])
  rownames(error_per_ind) <- IBD.list[[1]][[1]]$offspring
  
  return(list(maxL_IBD = IBDout,
              MML = mml.arr,
              error_per_ind = error_per_ind,
              errors = errors))
  
} #maxL_IBD

#' Generate a 'report' of predicted meiotic behaviour in an F1 population
#' @description Function to extract the chromosome pairing predictions as estimated by \code{\link{estimate_IBD}}. 
#' Apart from producing an overview of the pairing during parental meiosis (including counts of multivalents, per linkage group per parent),
#' the function also applies a simple chi-squared test to look for evidence of non-random pairing behaviour from the bivalent counts (deviations from a polysomic model)
#' @param IBD_list List of IBD probabilities as estimated by \code{\link{estimate_IBD}} using method 'hmm', or externally (e.g. using TetraOrigin)
#' @param visualise Logical, by default \code{FALSE}. If \code{TRUE}, a plot of the pairing results is produced per LG. In order to 
#' flag extreme deviations from the expected numbers (associated with polysomic inheritance, considered the Null hypothesis),
#' barplots are coloured according to the level of significance of the X2 test. Plots showing red bars indicate extreme deviations from
#' a polysomic pattern.
#' @param precision To how many decimal places should summed probabilities per bivalent pairing be rounded? By default 2. 
#' @return The function returns a nested list, with one element per linkage group in the same order as the input IBD list.
#' Per linkage group, a list is returned containing the following components:
#' \describe{
#' \item{P1_multivalents}{The count of multivalents in parent 1 (only relevant if \code{bivalent_decoding = FALSE} during IBD calculation)}
#' \item{P2_multivalents}{Similarly, the count of multivalents in parent 2}
#' \item{P1_pairing}{The counts of each bivalent pairing predicted in parent 1, with an extra column Pr(X2) which gives the p-value of 
#' the X2 test of the off-diagonal terms in the matrix. In the case of a tetraploid, pairing A with B automatically implies C with D pairing,
#' so the count table contains a lot of redundancy. The table should be read using both row and column names, so row A and column B corresponds
#' to the count of individuals with A and B pairing (and hence C and D pairing). In a hexaploid, A-B pairing does not imply a particular pairing
#' configuration in the remaining homologues. In this case, row A and column B is the count of individuals where A and B were predicted to have paired,
#' summed over all three bivalent configurations with A and B paired (AB-CD-EF, AB-CE-DF, AB-CF,DE).}
#' \item{P2_pairing}{Same as P1_pairing, except using parent 2}
#' \item{ploidy}{The ploidy of parent 1}
#' \item{ploidy2}{The ploidy of parent 2}
#' }
#' @examples 
#' data("IBD_4x")
#' mr.ls<-meiosis_report(IBD_list = IBD_4x)
#' @export
meiosis_report <- function(IBD_list,
                           visualise = FALSE,
                           precision = 2){
  
  IBD_list <- test_IBD_list(IBD_list)
  
  oldpar <- par(no.readonly = TRUE)    
  on.exit(par(oldpar))            # Thanks Gregor
  
  if(IBD_list[[1]]$IBDtype != "genotypeIBD") stop("Cannot produce meiosis report unless IBDs have been estimated using HMM method (or externally using TetraOrigin).")
  
  ploidy <- IBD_list[[1]]$ploidy
  ploidy2 <- IBD_list[[1]]$ploidy2
  
  biv_dec <- IBD_list[[1]]$biv_dec
  
  r1 <- 1:ploidy
  r2 <- (ploidy+1):(ploidy+ploidy2)
  
  P1grep <- sapply(LETTERS[r1], function(n) apply(sapply(LETTERS[r1], function(m) sort(c(n,m))),
                                                  2,paste0,collapse=""))
  
  P2grep <- sapply(LETTERS[r2], function(n) apply(sapply(LETTERS[r2], function(m) sort(c(n,m))),
                                                  2,paste0,collapse=""))
  
  pairing_list <- state_fun(ploidy,ploidy2,biv_dec)
  all.pairings <- unlist(sapply(pairing_list$state.ls, function(x) names(x)))
  
  if(!biv_dec){ #HMM approach was to define multivalent as association with all homologues only. 
    P1DR <- paste0(LETTERS[1:ploidy],collapse="")
    P2DR <- paste0(LETTERS[(ploidy+1):(ploidy+ploidy2)],collapse="")
  }
  
  out <- setNames(lapply(seq(length(IBD_list)), function(l){
    
    maxprobs <- sapply(IBD_list[[l]]$pairing,max)
    
    # plausible <- which(maxprobs >= plausible_pairing_prob)
    
    pairing.sums <- setNames(rep(0,length(all.pairings)),all.pairings)
    
    # for(y in IBD_list[[l]]$pairing[plausible]){
    #   pairing.sums[names(y)] <- pairing.sums[names(y)] + y
    # }
    
    for(y in IBD_list[[l]]$pairing){
      pairing.sums[names(y)] <- pairing.sums[names(y)] + y
    }
    
    P1countsDR <- P2countsDR <- 0
    
    if(!biv_dec){
      P1hitsDR <- grep(P1DR,all.pairings)
      P2hitsDR <- grep(P2DR,all.pairings)
      
      P1hits <- sapply(P1grep, function(b) grep(b,all.pairings[-P1hitsDR]))
      P2hits <- sapply(P2grep, function(b) grep(b,all.pairings[-P2hitsDR]))
      
      P1counts <- matrix(sapply(P1hits,function(i) sum(pairing.sums[-P1hitsDR][i])),ncol=ploidy,
                         dimnames = list(LETTERS[r1],LETTERS[r1]))
      P2counts <- matrix(sapply(P2hits,function(i) sum(pairing.sums[-P2hitsDR][i])),ncol=ploidy2,
                         dimnames = list(LETTERS[r2],LETTERS[r2]))
      
      P1countsDR <- sum(pairing.sums[P1hitsDR])
      P2countsDR <- sum(pairing.sums[P2hitsDR])
      
    } else{
      P1hits <- sapply(P1grep, function(b) grep(b,all.pairings))
      P2hits <- sapply(P2grep, function(b) grep(b,all.pairings))
      
      P1counts <- matrix(sapply(P1hits,function(i) sum(pairing.sums[i])),ncol=ploidy,
                         dimnames = list(LETTERS[r1],LETTERS[r1]))
      P2counts <- matrix(sapply(P2hits,function(i) sum(pairing.sums[i])),ncol=ploidy2,
                         dimnames = list(LETTERS[r2],LETTERS[r2]))
    }
    
    ## Apply a chi.sq test on off-diagonal terms
    
    P1.offdiag <- matrix(P1counts[upper.tri(P1counts) | lower.tri(P1counts)],
                         ncol = ploidy - 1,byrow = T)
    
    P2.offdiag <- matrix(P2counts[upper.tri(P2counts) | lower.tri(P2counts)],
                         ncol = ploidy2 - 1,byrow = T)
    
    ## Note than in tetraploid case there is redundancy here:
    P1.x2 <- suppressWarnings(apply(P1.offdiag,1,function(cnts) chisq.test(cnts)$p.val))
    P2.x2 <- suppressWarnings(apply(P2.offdiag,1,function(cnts) chisq.test(cnts)$p.val))
    
    P1out <- cbind(round(P1counts,precision),"Pr(X2)"=P1.x2)
    P2out <- cbind(round(P2counts,precision),"Pr(X2)"=P2.x2)
    
    return(list("P1_multivalents" = P1countsDR,
                "P2_multivalents" = P2countsDR,
                "P1_pairing" = P1out,
                "P2_pairing" = P2out,
                "ploidy" = ploidy,
                "ploidy2" = ploidy2
    ))
    
  }),names(IBD_list))
  
  ## Extra: add some visualisation options here also
  if(visualise){
    for(l in seq(length(IBD_list))){
      data <- out[[l]]
      
      P1expect <- sum(data$P1_pairing[1,1:ploidy])/(ploidy-1)
      P2expect <- sum(data$P2_pairing[1,1:ploidy2])/(ploidy2-1)
      
      P1plot <- data$P1_pairing[,1:ploidy] - P1expect
      diag(P1plot) <- 0
      P2plot <- data$P2_pairing[,1:ploidy2] - P2expect
      diag(P2plot) <- 0
      
      ## Make a colour ramp palette based on -log10(p) value
      p1lops <- -log10(data$P1_pairing[,ncol(data$P1_pairing)])
      p2lops <- -log10(data$P2_pairing[,ncol(data$P2_pairing)])
      
      
      #For now: introduce an arbitrary set of numbers to apply scale:
      seqpts <- c(seq(0,1.9,0.04),seq(2,7.5,0.5),seq(8,20,1),seq(21,120,4),500)
      colpal <- colorRampPalette(c("darkgreen","grey","red"))(length(seqpts))
      p1cols <- colpal[findInterval(p1lops,seqpts)]
      p2cols <- colpal[findInterval(p2lops,seqpts)]
      
      layout(
        rbind(c(1,1),
              c(2,2),
              apply(cbind(3:(ploidy+2),(ploidy + 3):(ploidy + ploidy2 + 2)),2,sort)
        ),
        heights = c(0.2,0.2,rep(1,ploidy))
      )
      
      par(mar = c(0,0,0,0))
      
      plot(NULL,xlim = c(0,1),ylim=c(0,1),axes=FALSE)
      xlen <- sum(abs(par("usr")[1:2]))
      text(x = 0.5*xlen,y = 0.5,names(IBD_list)[l],cex = 2.5, font = 2)
      
      
      plot(NULL,xlim = c(0,1),ylim=c(0,1),axes=FALSE)
      xlen <- sum(abs(par("usr")[1:2]))
      text(x = 0.25*xlen,y = 0.5,"P1",cex = 2, font = 2)
      text(x = 0.75*xlen,y = 0.5,"P2",cex = 2, font = 2)
      
      par(mar = c(4,4,0.5,0.5))
      sapply(1:nrow(P1plot), function(r) 
        barplot(P1plot[r,],
                names.arg = paste0(rownames(P1plot)[r],
                                   colnames(P1plot)),
                ylab = "Deviations",
                col = p1cols[r]))
      
      sapply(1:nrow(P2plot), function(r) 
        barplot(P2plot[r,],
                names.arg = paste0(rownames(P2plot)[r],
                                   colnames(P2plot)),
                ylab = "Deviations",
                col = p2cols[r]))
    }
  }
  
  return(out)
} #meiosis_report




#' Plot the results of QTL scan. 
#' @description Up to package v.0.0.9, there were three plotting functions for the output of \code{QTLscan}, namely \code{plotQTL}, \code{plotLinearQTL} and \code{plotLinearQTL_list}.
#' Since release 0.1.0, the functionality of all three functions has been combined into a single general plotting function, named \code{plotQTL}.
#' The plot layout is now specified by a new argument \code{layout}, allowing the user to plot results for single chromosomes separately, or together either adjacently or in a grid layout.
#' Results from multiple analyses can be overlaid. Previously, it was possible to call the function \code{plotQTL} multiple times and overlay subsequent plots using the argument \code{overlay = TRUE}.
#' This approach is no longer supported. Instead, if multiple results are to be overlaid, they can be provided as a list of \code{QTLscan} or \code{singleMarkerRegression} outputs. Note however that this
#' is only possible using the default layout. If significance thresholds are
#' present, the default behaviour is to rescale LOD values so that multiple plots can be combined with overlapping signficance thresholds. This rescaling behaviour can also be
#' disabled (by setting \code{rescale = FALSE}). Note that not all arguments may be appropriate for all layouts.
#' @param LOD_data Output of \code{\link{QTLscan}} function. If you wish to overlay multiple genome-wide \code{QTLscan} outputs, then first compile these into a single list and 
#' pass this to \code{LOD_data}, for example \code{LOD_data = list(qtl_res1, qtl_res2)}. If this is passed as a named list and \code{add_legend = TRUE}, these names will be used
#' in the legend as well. 
#' @param layout There are three possible plot layouts - single chromosome plots ("s"), genome-wide plots arranged adjacently in a linear fashion ("l") which is
#' also the default, and genome-wide plots arranged in a grid ("g"), i.e. a grid of single chromosome plots. In the latter case, a suitable grid dimension will be determined
#' based on the number of linkage groups detected in \code{LOD_data}. If you wish to overlay results from multiple multi-chromosome analyses, use the default layout.  
#' @param inter_chm_gap The gap size (in units of cM) between successive chromosomes when \code{layout = "l"}. By default a gap of 5 cM is used. Normally the user should not need to change this.
#' @param ylimits Use to specify ylimits of plot region, though by default \code{NULL} in which case a suitable plot region is automatically used.
#' @param sig.unit Label to use on the y-axis for significance units, by default assumed to be LOD score.
#' @param plot_type Plots can be either in line drawings ("lines", default) or scatter plot format ("points").
#' @param colour Vector of colours to be used in the plotting. A default set of 4 colours is provided, the first of which is used when results from a single QTL scan are to be plotted.
#' @param add_xaxis Should an x-axis be drawn? If multiple QTL analyses are performed on different traits, specifying this to be \code{FALSE}
#' and using \code{par(mar=c(0,4.1,4.1,2.1))} allows subsequent plots to be neatly stacked.
#' @param add_rug Logical, by default \code{TRUE} - should original marker points be added to plot?
#' @param add_thresh Logical, by default \code{TRUE} - should a significance threshold be added to plot?
#' @param override_thresh By default \code{NULL}. Can be used to specify a (numeric) value for the significance threshold, overriding any stored in \code{LOD_data}. If you wish to override
#' thresholds for multiple analyses (so, when \code{LOD_data} is a list of QTL outputs), can also provide a vector of numeric values here.
#' @param thresh.lty Gives user control over the line type of the significance threshold to be drawn. Default threshold lty is 3.
#' @param thresh.lwd Gives user control over the line width of the significance threshold to be drawn. Default threshold lwd is 2.
#' @param thresh.col Gives user control over the line colour of the significance threshold to be drawn. Default threshold colour is dark red. If plotting multiple analyses with \code{rescale = FALSE},
#' it can be useful to provide the same colours to this argument as to \code{colour}, so that LOD profiles can be linked to their respective LOD thresholds.
#' @param return_plotData Logical, by default \code{FALSE}. If \code{TRUE}, then the x and y coordinates of the plot data are returned when \code{layout = "l"},
#' which can be useful for subsequent plot manipulations and overlays. For other layouts, no plot data is returned.
#' @param show_thresh_CI Logical, by default \code{FALSE}. Should confidence interval bounds around LOD threshold be shown if available? If \code{LOD_data} is a list from multiple analyses,
#' this option is ignored to prevent plot becoming too cluttered.
#' @param use_LG_names Logical, by default \code{TRUE}. Should original character LG names (the names of list \code{LOD_data}) be used as axis labels? If \code{FALSE}, numbering is used instead.
#' @param axis_label.cex Argument to adjust the size of the axis labels. Can be useful if there are many linkage groups to plot
#' @param custom_LG_names Option to specify a vector that contains custom linkage group names. By default \code{NULL}. See previous argument \code{use_LG_names}, which is the usual manner
#' for controlling x-axis labels.
#' @param LGdiv.col Colour of dividing lines between linkage groups when \code{layout = "l"}, by default grey.
#' @param ylab.at Distance from the y-axis to place label (by default at 2.5 points)
#' @param highlight_positions Option to include a (list of) positions to highlight (e.g. peak QTL positions). Each list element should be a 2-column data.frame with columns giving
#' the linkage group numbers (numeric) and the corresponding cM positions (numeric) to highlight. If \code{LOD_data} is the result of a single genome-wide scan, it is also possible to just directly provide the
#' 2-column data.frame (again, with column 1 containing linkage group numbers and column 2 containing corresponding cM positions).
#' If \code{LOD_data} has been provided as a list of multiple analyses, you may wish to highlight different positions from each analysis. Then \code{highlight_positions} should also be a list 
#' of the same length and in the same order as \code{LOD_data}. Each data.frame of positions will be coloured in the same colour as the LOD output.
#' If no position is to be highlighted for some analyses, add the corresponding list element as \code{NULL}. For example, if you wish to highlight positions for analyses 1 and 3 in a 3-analysis overlay, 
#' then use something like \code{highlight_positions = list(data.frame(lg = 1, cM = 50),NULL,data.frame(lg=c(2,3),cM=c(11,99)))}.
#' The default setting is \code{NULL}, meaning no positions are highlighted.
#' @param mainTitle Option to supply vector of plot titles if \code{layout = "s"} or \code{layout = "g"}. Argument ignored if using the default layout. Single character vector also allowed and will be recycled.
#' For no plot titles, leave as default, i.e. \code{FALSE}
#' @param rescale If results from multiple analyses are to be overlaid and different significance thresholds are detected, then by default plots will be rescaled so that threshold lines overlap.
#' This behaviour can be disabled by setting \code{rescale = FALSE}. 
#' @param \dots Arguments passed to \code{\link{plot}}, and \code{\link{lines}} or \code{\link{points}} as appropriate (see argument \code{plot_type}).
#' @aliases plotLinearQTL plotLinearQTL_list
#' @return The plot data, if return_plotData = TRUE. Otherwise \code{NULL}. Output is returned invisibly
#' @examples
#' \dontrun{
#' data("qtl_LODs.4x")
#' plotQTL(LOD_data = qtl_LODs.4x,layout = "l")
#' }
#' @export
plotQTL <- function(LOD_data,
                    layout = "l",
                    inter_chm_gap = 5, 
                    ylimits = NULL,
                    sig.unit = "LOD",
                    plot_type = "lines",
                    colour = c("black","red","dodgerblue","sienna4"),
                    add_xaxis = TRUE,
                    add_rug = TRUE,
                    add_thresh = TRUE,
                    override_thresh = NULL,
                    thresh.lty = 3,
                    thresh.lwd = 2,
                    thresh.col = "darkred",
                    return_plotData = FALSE,
                    show_thresh_CI = FALSE,
                    use_LG_names = TRUE,
                    axis_label.cex = 1,
                    custom_LG_names = NULL,
                    LGdiv.col = "gray42",
                    ylab.at = 2.5,
                    highlight_positions = NULL,
                    mainTitle = FALSE,
                    rescale = TRUE,
                    ...){
  layout <- match.arg(layout, choices = c("s","l","g"))
  
  oldpar <- par(no.readonly = TRUE)    
  on.exit(par(oldpar))            # Thanks to Gregor, CRAN
  
  overlay <- FALSE #whether user wants to add extra plots, by default FALSE
  multiplot <- NULL #just in case
  linear_list <- TRUE #for plotting multiple analyses with layout "l" and sig thresh available for all
  
  if(list.depth(LOD_data) == 3) {
    nran <- length(LOD_data) #number of analyses
    overlay <- TRUE
    
    ## Check whether input data all has the same number of linkage groups (and extract intersect if not)
    nLG <- sapply(LOD_data, function(x) length(unique(x$QTL.res$chromosome)))

    if(length(unique(nLG)) != 1) {
      
      warning("Variable numbers of LGs detected in LOD_data! Attempting to resolve...")
      LG_ls <- lapply(LOD_data, function(x) unique(x$QTL.res$chromosome))
      core_LG <- Reduce(intersect, LG_ls) #gets the intersection of all list elements
      
      for(i in 1:length(LOD_data)){
        LOD_data[[i]]$QTL.res <- LOD_data[[i]]$QTL.res[LOD_data[[i]]$QTL.res$chromosome %in% core_LG,]
      }
      
      #Re-calculate:
      nLG <- sapply(LOD_data, function(x) length(unique(x$QTL.res$chromosome)))
      }
    nLG <- unique(nLG)
    
    ## Disable the linear list approach which rescales:
    if(!rescale) {
      linear_list <- FALSE
      
      #Vectorise thresh.col
      if(length(thresh.col) == 1){
        message("Multiple analyses being overlaid with only a single thresh.col. Suggest to re-plot with different threshold colours!")
        thresh.col <- rep(thresh.col,nran)
      } else if(length(thresh.col) != length(LOD_data)) {
        stop(paste("If plotting multiple threshold colours, please supply as a vector of length",length(LOD_data)))
      }
      
    }
    ## If there are any thresholds missing it is not possible to rescale, so also disable linear list approach:
    if(any(sapply(sapply(LOD_data, function(x) x$Perm.res),is.null)) & is.null(override_thresh)) linear_list <- FALSE
    
    ## Vectorise some arguments as they may be indexed by [i] later on
    if(length(plot_type) == 1){
      plot_type <- sapply(1:nran,function(n) match.arg(plot_type, choices = c("lines","points")))
    } else if(length(plot_type) != length(LOD_data)) {
      stop(paste("If multiple plot_type required, please supply as a vector of length",length(LOD_data)))
    } else{
      plot_type <- sapply(plot_type,match.arg,choices = c("lines","points"))
    }
    
    if(length(colour) == 1){
      colour <- rep(colour,nran)
    } else if(length(colour) < length(LOD_data)) {
      stop(paste("By default only 4 plotting colours are provided. In your case, please supply 'colour' as a vector of length",length(LOD_data)))
    } else{
      colour <- colour[1:length(LOD_data)]
    }
    
  } else{
    ## Results of a single QTLscan to be plotted
    nLG <- length(unique(LOD_data$QTL.res$chromosome)) #there may be chr0 data
    # nLG <- length(LOD_data$LG_names)
  }
  
  if(layout == "g"){
    
    ## Have to define the grid layout
    base <- floor(sqrt(nLG)) #this is the ncols to remain a roughly square grid
    if(nLG - base^2 > 0){
      nr <- base + 1
      nc <- base
      # excess <- (base)*(base + 1) - nLG
    } else{
      nr <- nc <- base
      # excess <- 0
    }
    multiplot <- c(nr,nc)
  }
  
  ## Pass to the 3 possible plotting functions:
  if(layout == "l" & !overlay){
    plotQTL.genome_linear(LOD_data = LOD_data,
                          inter_chm_gap = inter_chm_gap,
                          overlay = FALSE,
                          ylimits = ylimits,
                          sig.unit = sig.unit,
                          plot_type = plot_type,
                          colour = colour,
                          add_xaxis = add_xaxis,
                          add_rug = add_rug,
                          add_thresh = add_thresh,
                          override_thresh = override_thresh,
                          thresh.lty = thresh.lty,
                          thresh.lwd = thresh.lwd,
                          thresh.col = thresh.col,
                          return_plotData = return_plotData,
                          show_thresh_CI = show_thresh_CI,
                          use_LG_names = use_LG_names,
                          axis_label.cex = axis_label.cex,
                          custom_LG_names = custom_LG_names,
                          LGdiv.col = "gray42",
                          ylab.at = ylab.at,
                          highlight_positions = highlight_positions,
                          ...)
    
  } else if(layout == "l" & overlay){
    
    if(linear_list){
      plotQTL.genome_linear_list(LOD_data.ls = LOD_data,
                                 inter_chm_gap = inter_chm_gap,
                                 ylimits = ylimits,
                                 sig.unit = sig.unit,
                                 plot_type = plot_type,
                                 add_xaxis = add_xaxis,
                                 add_rug = add_rug,
                                 colours = colour, #the plotLinearQTL_list function used colours, not colour for this argument
                                 ylab.at = ylab.at, 
                                 thresh.lty = thresh.lty,
                                 thresh.lwd = thresh.lwd,
                                 thresh.col = thresh.col,
                                 return_plotData = return_plotData,
                                 highlight_positions = highlight_positions,
                                 LGdiv.col = LGdiv.col,
                                 use_LG_names = use_LG_names,
                                 axis_label.cex = axis_label.cex,
                                 custom_LG_names = custom_LG_names,
                                 override_thresh = override_thresh,
                                 ...)
    } else{
      ## the user would like to overlay plots but does not want to rescale, or thresholds were absent
      # Determine the plotting region
      if(is.null(ylimits)){
        maxy <- max(unlist(sapply(lapply(LOD_data, "[[","QTL.res"),"[[","LOD")))
        ylimits = c(0,maxy)
      }
      
      plotQTL.genome_linear(LOD_data = LOD_data[[1]],
                            inter_chm_gap = inter_chm_gap,
                            overlay = FALSE,
                            ylimits = ylimits,
                            sig.unit = sig.unit,
                            plot_type = plot_type[1],
                            colour = colour[1],
                            add_xaxis = add_xaxis,
                            add_rug = add_rug,
                            add_thresh = add_thresh,
                            override_thresh = override_thresh,
                            thresh.lty = thresh.lty,
                            thresh.lwd = thresh.lwd,
                            thresh.col = thresh.col[1],
                            return_plotData = return_plotData,
                            show_thresh_CI = show_thresh_CI,
                            use_LG_names = use_LG_names,
                            axis_label.cex = axis_label.cex,
                            custom_LG_names = custom_LG_names,
                            LGdiv.col = "gray42",
                            ylab.at = ylab.at,
                            highlight_positions = highlight_positions,
                            ...)
      
      ## Overlay the rest of the plots
      for(i in 2:length(LOD_data)){
        plotQTL.genome_linear(LOD_data = LOD_data[[i]],
                              inter_chm_gap = inter_chm_gap,
                              overlay = TRUE,
                              ylimits = ylimits,
                              sig.unit = sig.unit,
                              plot_type = plot_type[i],
                              colour = colour[i],
                              add_xaxis = add_xaxis,
                              add_rug = add_rug,
                              add_thresh = add_thresh,
                              override_thresh = override_thresh,
                              thresh.lty = thresh.lty,
                              thresh.lwd = thresh.lwd,
                              thresh.col = thresh.col[i],
                              return_plotData = return_plotData,
                              show_thresh_CI = show_thresh_CI,
                              use_LG_names = use_LG_names,
                              axis_label.cex = axis_label.cex,
                              custom_LG_names = custom_LG_names,
                              LGdiv.col = "gray42",
                              ylab.at = ylab.at,
                              highlight_positions = highlight_positions,
                              ...)
      }
    }
    

  } else if(layout == "s" & !overlay){
    
    plotQTL.single_LG(LOD_data = LOD_data,
                      ylimits = ylimits,
                      multiplot = NULL,
                      plot_type = plot_type,
                      overlay = FALSE,
                      add_xaxis = add_xaxis,
                      add_rug = add_rug,
                      mainTitle = mainTitle,
                      override_thresh = override_thresh,
                      thresh.lty = thresh.lty,
                      thresh.lwd = thresh.lwd,
                      thresh.col = thresh.col,
                      colour = colour[1],
                      axis_label.cex = axis_label.cex,
                      ...)
    
  } else if(layout == "s" & overlay){
    ## Draw first plot:
    if(length(LOD_data[[1]]$LG_names) > 1) stop("It is not possible to overlay multiple analyses and multiple chromosomes with layout = 's'. Use layout 'l' instead.")
    plotQTL.single_LG(LOD_data = LOD_data[[1]],
                      ylimits = ylimits,
                      multiplot = NULL,
                      plot_type = plot_type[1],
                      overlay = FALSE,
                      add_xaxis = add_xaxis,
                      add_rug = add_rug,
                      mainTitle = mainTitle,
                      override_thresh = override_thresh,
                      thresh.lty = thresh.lty,
                      thresh.lwd = thresh.lwd,
                      thresh.col = thresh.col,
                      colour = colour[1],
                      axis_label.cex = axis_label.cex,
                      ...)
    ## Overlay subsequent plots:
    for(i in 2:length(LOD_data)){
      tempdata <- LOD_data[[i]]
      if(length(tempdata$LG_names) > 1) stop("It is not possible to overlay multiple analyses and multiple chromosomes with layout = 's'. Use layout 'l' instead.")
      
      plotQTL.single_LG(LOD_data = tempdata,
                        ylimits = ylimits,
                        multiplot = NULL,
                        plot_type = plot_type[i],
                        overlay = TRUE,
                        add_xaxis = add_xaxis,
                        add_rug = add_rug,
                        mainTitle = mainTitle,
                        override_thresh = override_thresh,
                        thresh.lty = thresh.lty,
                        thresh.lwd = thresh.lwd,
                        thresh.col = thresh.col,
                        colour = colour[i],
                        axis_label.cex = axis_label.cex,
                        ...)
    }
  } else if(layout == "g" & !overlay){
    
    plotQTL.single_LG(LOD_data = LOD_data,
                      ylimits = ylimits,
                      multiplot = multiplot,
                      plot_type = plot_type,
                      overlay = FALSE,
                      add_xaxis = add_xaxis,
                      add_rug = add_rug,
                      mainTitle = mainTitle,
                      override_thresh = override_thresh,
                      thresh.lty = thresh.lty,
                      thresh.lwd = thresh.lwd,
                      thresh.col = thresh.col,
                      colour = colour[1],
                      axis_label.cex = axis_label.cex,
                      ...)
    
  } else if(layout == "g" & overlay){
    stop("It is not possible to overlay multiple analyses and multiple chromosomes with layout = 'g'. Use layout 'l' instead.")
  }
  
} #plotQTL


#' Deprecated function name for plotLinearQTL, still maintained in package version 0.1.0
#' @rdname plotQTL
#' @export
plotLinearQTL <- plotQTL

#' Deprecated function name for plotLinearQTL_list, still maintained in package version 0.1.0
#' @rdname plotQTL
#' @export
plotLinearQTL_list <- plotQTL



#' Plot the recombination landscape across the genome
#' @description Function which visualises the recombination landscape in two ways: per linkage group, and per individual. 
#' For the first analysis, a rudimentary spline is also fitted to estimate the recombination rate along a grid of positions defined by \code{gap},
#' which is also returned by the function.
#' @param recombination_data Data on predicted recombination events, as returned by the function \code{\link{count_recombinations}}
#' @param gap The size (in cM) of the gap used to define the grid of positions to define the window in which to estimate 
#' recombination rate. By default 1 cM. Interpolated positions are taken to be the centre of an interval, so a 1 cM gap would
#' result in predictions for positions 0.5 cM, 1.5 cM etc.
#' @param plot_per_LG Logical argument, plot recombination events per linkage group? By default \code{TRUE}.
#' @param plot_per_ind Logical argument, plot recombination events per individual? By default \code{TRUE}.
#' @param \dots Option to pass extra arguments to the \code{\link{plot}} function for the per_LG plots. This may lead
#' to conflicts with arguments already declared internally (such as \code{main} for example).
#' @return A list with two elements, \code{per_LG} and \code{per_individual}. The first of these is itself a list with the same length as \code{recombination_data}, giving the estimated recombination rates along the linkage group.
#' This rate is simply estimated as the (weighted) count of recombination breakpoints divided by the population size.
#' @examples 
#' data("Rec_Data_4x")
#' plotRecLS(Rec_Data_4x)
#' @export
plotRecLS <- function(recombination_data,
                      plot_per_LG = TRUE,
                      plot_per_ind = TRUE,
                      gap = 1,
                      ...){
  
  ## First, per LG:
  perLG.ls <- setNames(lapply(1:length(recombination_data), function(ch=1){
    
    recdat <- recombination_data[[ch]]$recombinations
    N <- length(recdat) #number of offspring
    
    ## Extract all recombination breakpoints and their probabilities into a 2 column matrix:
    rec.mat <- do.call(rbind,lapply(1:N, function(n){
      # print(n)
      f1 <- recdat[[n]]
      
      if(!is.na(f1[1])){
        val.dat <- lapply(f1, function(val = f1[[1]]){
          temp <- unlist(val$recombination_breakpoints)
          temp <- temp[!is.na(temp)]
          if(length(temp) > 0)  return(as.data.frame(cbind(cM = temp,
                                                           # Pr = val$Valent_probability,
                                                           Ind = names(recdat[n]))))
          
        })
        
        if(!is.null(val.dat)){
          
          if(length(val.dat) > 1){
            # in this case there are multiple possible valency solutions - combine them:
            val.combined <- do.call(rbind,val.dat)
            
            if(is.null(val.combined)){
              val.out <- NULL
            } else if(ncol(val.combined) == 2){
              cM.unique <- sort(unique(as.numeric(val.combined[,"cM"])))
              # val.out <- as.data.frame(matrix(0,nrow = length(cM.unique),ncol = 3,
              # dimnames = list(NULL,c("cM","Pr","Ind"))))
              
              val.out <- as.data.frame(matrix(0,nrow = length(cM.unique),ncol = 2,
                                              dimnames = list(NULL,c("cM","Ind"))))
              
              val.out[,"cM"] <- cM.unique
              # val.out[,"Pr"] <- sapply(cM.unique, function(cm) sum(as.numeric(val.combined[val.combined[,"cM"] == cm,"Pr"])))
              val.out[,"Ind"] <- names(recdat[n])
            } else{
              val.out <- NULL
            }
            
          } else{
            val.out <- val.dat[[1]]
          }
          
          return(val.out)
        }
      }
    })
    )
    
    # rec.mat[,"Pr"] <- as.numeric(rec.mat[,"Pr"])
    rec.mat[,"cM"] <- as.numeric(rec.mat[,"cM"])
    
    # Generate recombination rates per bin of size gap:
    cM_bins <- seq(0, ceiling(max(recombination_data[[ch]]$map$position)), gap)
    rec.mat <- cbind(rec.mat, bin = findInterval(rec.mat[,"cM"],cM_bins))
    
    # plotdata <- tapply(rec.mat[,"Pr"],rec.mat[,"bin"],sum)
    plotdata <- table(rec.mat[,"bin"])
    plotdata <- cbind("cM" = as.numeric(names(plotdata)) - gap/2,
                      "Count" = as.numeric(plotdata))
    
    splined.counts <- predict(smooth.spline(plotdata[,"cM"],plotdata[,"Count"]),cM_bins)$y
    
    if(plot_per_LG){
      plot(plotdata[,"cM"],plotdata[,"Count"],
           xlim = c(0,max(recombination_data[[ch]]$map$position)),
           ylim = c(0,max(plotdata[,"Count"])),
           xlab = "cM", ylab = "Count",
           main = names(recombination_data)[ch],
           ...)
      lines(cM_bins,splined.counts,...)
    }
    
    return(list(rec.mat = cbind(rec.mat,"LG" = ch),
                rec.rates = cbind(cM = cM_bins, rate = round(splined.counts/N,3))))
  }),
  names(recombination_data))
  
  ## Next, summarise the per-individual counts (i.e combined counts across linkage groups)
  N <- length(recombination_data[[1]]$recombinations)
  
  combined_rec.mat <- do.call(rbind,lapply(perLG.ls,"[[","rec.mat"))
  
  ## Generate counts of recombination events per individual:
  offspring_counts <- sort(table(combined_rec.mat$Ind))
  
  if(plot_per_ind){
    plot(NULL,xlim = c(0, length(offspring_counts)),
         ylim = range(offspring_counts),
         xlab = "Individual",ylab = "Count (all LGs)",axes=F,
         main = "Genome-wide count of recombinations per individual")
    
    axis(1,las=2,at = 1:length(offspring_counts),
         labels = names(offspring_counts),
         cex.axis = 0.6)
    axis(2,las = 1);box()
    points(1:length(offspring_counts),offspring_counts,
           ...)
  }
  
  ## Return output of both dimensions:
  invisible(list(per_LG = perLG.ls,
                 per_individual = offspring_counts))
  
} #plotRecLS


#' Function to determine the percentage variance explained (PVE) of a (maximal) QTL model, and explore sub-models.
#' @description This function builds a (maximal) QTL model from previously detected QTL peaks and outputs the percentage variance explained (PVE)
#' of the full QTL model and all sub-models. It uses a similar approach to the fitting of genetic co-factors in the function \code{\link{QTLscan}}.
#' The PVE is very similar to but not exactly equal to the adjusted R2 returned in \code{QTLscan} at each position (and note: in the former case, these
#' R2 values are per-locus, while this function can estimate the PVE combined over multiple loci). The discrepancy has to do with how PVE is calculated 
#' using the formula 100(1 - RSS0/RSS1), where RSS0 and RSS1 are the residual sums of squares of the NULL and QTL models, respectively.
#' @param IBD_list List of IBD probabilities
#' @param Phenotype.df A data.frame containing phenotypic values
#' @param genotype.ID The colname of \code{Phenotype.df} that contains the offspring identifiers (F1 names)
#' @param trait.ID The colname of \code{Phenotype.df} that contains the response variable to use in the model
#' @param block The blocking factor to be used, if any (must be colname of \code{Phenotype.df}). By default \code{NULL}, in 
#' which case no blocking structure (for unreplicated experiments)
#' @param QTL_df A 2-column data frame of previously-detected QTL; column 1 gives linkage group identifiers, 
#' column 2 specifies the cM position of the QTL. If not specified, an error results. It can be convenient to generate a compatible
#' data.frame by first running the function \code{\link{check_cofactors}} to build a multi-QTL model.
#' @param prop_Pheno_rep The minimum proportion of phenotypes represented across blocks. If less than this, the individual is 
#' removed from the analysis. If there is incomplete data, the missing phenotypes are imputed using the mean values across the recorded observations.
#' @param log Character string specifying the log filename to which standard output should be written. If \code{NULL} log is send to stdout.
#' @param verbose Should messages be written to standard output?
#' @return
#' A list with percentage variance explained of maximal QTL model and all sub-models
#' @examples
#' data("IBD_4x","Phenotypes_4x")
#' PVE(IBD_list = IBD_4x,
#'     Phenotype.df = Phenotypes_4x,
#'     genotype.ID = "geno",trait.ID = "pheno",
#'     block = "year",
#'     QTL_df = data.frame(LG=1,cM=12.3))
#' @export
PVE <- function(IBD_list,
                Phenotype.df,
                genotype.ID,
                trait.ID,
                block = NULL,
                QTL_df = NULL,
                prop_Pheno_rep = 0.5,
                log = NULL,
                verbose = FALSE){
  
  IBD_list <- test_IBD_list(IBD_list)
  
  if(is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  ## Run some initial checks:
  if(!is.null(block)) if(!block %in% colnames(Phenotype.df)) stop("Blocking factor incorrectly specified.")
  if(!trait.ID %in% colnames(Phenotype.df)) stop("Trait colname incorrectly specified.")
  
  blockTF <- FALSE
  if(!is.null(block)) blockTF <- TRUE
  if(is.null(QTL_df)) stop("Maximal QTL model to be tested not supplied!")
  if(nrow(QTL_df) == 1) warning("Maximal QTL model has a single QTL!\nA single QTL scan using function QTLscan provides comparable results!")
  
  ploidy <- IBD_list[[1]]$ploidy
  ploidy2 <- IBD_list[[1]]$ploidy2
  
  ## Check the list depth:
  # listdepth <- list.depth(IBD_list) #polyqtlR:::list.depth(IBD_list)
  # if(listdepth != 3) stop("Unexpected input - IBD_list is expected to be a nested list representing 1 or more chromosomes! Please check input.")
  
  IBDarray <- IBD_list[[1]]$IBDarray
  IBDtype <- IBD_list[[1]]$IBDtype
  bivalent_decoding <- IBD_list[[1]]$biv_dec
  gap <- IBD_list[[1]]$gap
  
  if(IBDtype == "genotypeIBD"){
    # if(IBD_list[[1]]$method == "hmm_TO"){
    #   ## We are using the output of TetraOrigin
    #   GenCodes <- t(sapply(IBD_list[[1]]$genocodes,function(x) as.numeric(strsplit(as.character(x),"")[[1]])))
    #   Nstates <- nrow(GenCodes)
    #   indicatorMatrix <- matrix(0,nrow=Nstates, ncol=ploidy + ploidy2)
    #   for(r in 1:nrow(indicatorMatrix)){
    #     for(h in 1:ncol(GenCodes)){
    #       indicatorMatrix[r,GenCodes[r,h]] <- indicatorMatrix[r,GenCodes[r,h]] + 1
    #     }
    #   }
    #   
    # } else{
    ## We are using the output of estimate_IBD
    GenCodes <- t(sapply(sapply(IBD_list[[1]]$genocodes,function(x) strsplit(x,"")), function(y) sapply(y,function(z) which(letters == z))))
    Nstates <- nrow(GenCodes)
    indicatorMatrix <- matrix(0,nrow=Nstates, ncol=ploidy + ploidy2)
    for(r in 1:nrow(indicatorMatrix)){
      for(h in 1:ncol(GenCodes)){
        indicatorMatrix[r,GenCodes[r,h]] <- indicatorMatrix[r,GenCodes[r,h]] + 1
      }
    }
    
    # }
    
  } else if(IBDtype == "haplotypeIBD"){
    Nstates <- ploidy + ploidy2 #this is not actually needed..
    indicatorMatrix <- NULL
  } else{
    stop("Unable to determine IBD type, please check IBD_list elements have a valid $IBDtype tag.")
  }
  
  ## Could just check for equality here also, modulo check was from older version...
  if(dim(IBDarray)[2] %% Nstates != 0) stop("Incompatible input detected..")
  
  ## First get the starting set - the ones that were phenotyped and genotyped:
  if(any(is.na(Phenotype.df[,genotype.ID]))){
    warning("Missing genotype.ID values detected. Removing and proceeding without...")
    Phenotype.df <- Phenotype.df[!is.na(Phenotype.df[,genotype.ID]),]
  }
  
  phenoGeno0 <- intersect(
    dimnames(IBDarray)[[3]],
    unique(Phenotype.df[,genotype.ID]))
  
  if(verbose) write(paste("There are",length(phenoGeno0),"individuals with matching identifiers between phenotypic and genotypic datasets.\n"),log.conn)
  
  nr.block <- 1 #there is always 1 block...
  
  if(blockTF){
    write(paste("Summary of trait:",trait.ID,"over blocks\n"),log.conn)
    summTab <- tapply(Phenotype.df[,trait.ID],Phenotype.df[,block],summary)
    
    ## If there are no Missing values for some blocks, add NA = 0:
    for(i in seq(length(summTab))){
      if(!"NA's" %in% names(summTab[[i]])) {
        summTab[[i]] <- setNames(c(summTab[[i]],0),c(names(summTab[[i]]),"NA's"))
      }
    }
    summTab.df <- do.call("rbind", lapply(seq(length(summTab)), function(l) summTab[[l]]))
    rownames(summTab.df) <- names(summTab)
    write(knitr::kable(summTab.df), file=log.conn)
    write("\n",log.conn)
    df.columns <- c(genotype.ID,trait.ID,block)
    
    ## Make sure there is only a single observation per genotype per block:
    Phenotype.df <- do.call("rbind", lapply(levels(as.factor(Phenotype.df[,block])), function(bl){
      do.call("rbind",
              sapply(unique(Phenotype.df[Phenotype.df[,block]==bl,genotype.ID]), function(gen){
                temp.mean <- mean(Phenotype.df[Phenotype.df[,block]==bl & Phenotype.df[,genotype.ID] == gen,trait.ID], na.rm = TRUE)
                temp.out <- Phenotype.df[Phenotype.df[,block]==bl & Phenotype.df[,genotype.ID] == gen,][1,,drop=FALSE]
                temp.out[,trait.ID] <- round(temp.mean,1)
                return(temp.out)}, simplify = FALSE))
    }))
    
    Phenotype.df[,block] <- as.factor(Phenotype.df[,block])
    
    ## Get the counts of observations of the offspring over blocks:
    count.fun <- function(PGeno){
      colSums(sapply(PGeno, function(gen) sapply(levels(as.factor(Phenotype.df[,block])), function(bl) {
        out <- FALSE
        if(gen %in% Phenotype.df[Phenotype.df[,block]==bl,genotype.ID]){
          if(!is.na(Phenotype.df[Phenotype.df[,block]==bl & Phenotype.df[,genotype.ID] == gen,trait.ID])) out <- TRUE
        }
        return(out)
      }
      )))
    } #count.fun()
    
    counts <- count.fun(phenoGeno0)
    block.lev <- levels(as.factor(Phenotype.df[,block]))
    nr.block <- length(block.lev)
    
    rem <- which(counts < nr.block*prop_Pheno_rep)
    
    if(length(rem)!=0){
      
      write(paste("\nThe following individuals had fewer than",nr.block*prop_Pheno_rep,"observations over",nr.block,"blocks and were removed from the analysis:\n"),log.conn)
      write(knitr::kable(as.data.frame(matrix(c(names(rem),rep(" ",4-length(rem)%%4)),ncol=4))),file = log.conn)
      phenoGeno1 <- setdiff(phenoGeno0,names(rem))
      
    } else{
      phenoGeno1 <- phenoGeno0
    }
    
    ## Next, impute for the counts that are left:
    counts <- count.fun(phenoGeno1)
    imputers <- which(counts != nr.block)
    
    if(length(imputers) != 0) {
      
      #First take out the complete data - will add at the end:
      pheno0 <- Phenotype.df[Phenotype.df[,genotype.ID] %in% setdiff(phenoGeno1,names(imputers)),df.columns]
      
      warning(paste0("Estimating block effects with ",
                     length(setdiff(phenoGeno1,names(imputers))),
                     " of the total population size of ",
                     length(phenoGeno1),"!"))
      
      #Estimate block effects:
      lmDat <- pheno0[,2:3]
      colnames(lmDat) <- c("Pheno","Block")
      lm.block <- lm(Pheno ~ Block, data = lmDat)
      block.effects <- c(0,lm.block$coefficients[2:length(lm.block$coefficients)]) #Add a 0 temporarily, will remove later if not needed..
      
      pheno_impute <- NULL
      
      for(gen in names(imputers)){
        temp <- Phenotype.df[Phenotype.df[,genotype.ID]==gen,df.columns]
        if(nrow(temp)!=nr.block){
          #A bit of overkill, but general for any number of missing rows (previously, only 1 missing row was allowed.)
          miss.bl <- setdiff(block.lev,temp[,3])
          temp <- rbind(temp,setNames(data.frame(
            matrix(c(rep(gen,length(miss.bl)),
                     rep(NA,length(miss.bl)),
                     miss.bl),ncol = 3)),names(temp)))
        }
        temp <- temp[order(temp[,3]),]
        levelNA <- which(is.na(temp[,2]))
        levelOK <- setdiff(seq(nr.block),levelNA)
        
        temp[levelNA,2] <- mean(as.numeric(temp[levelOK,2]) - block.effects[levelOK]) + block.effects[levelNA]
        rownames(temp) <- paste(temp[1,1],temp[,3], sep="_")
        pheno_impute <- rbind(pheno_impute,temp)
      }
      
      pheno <- rbind(pheno0,pheno_impute)
      pheno[,2] <- as.numeric(pheno[,2])
      
    } else{
      pheno <- Phenotype.df[Phenotype.df[,genotype.ID] %in% phenoGeno1,df.columns]
      pheno[,2] <- as.numeric(pheno[,2])
    }
    
    pheno <- pheno[order(pheno[,3],pheno[,1]),]
    colnames(pheno)[3] <- "Block"
  } else{ #No blocks included.
    pheno <- Phenotype.df[Phenotype.df[,genotype.ID] %in% phenoGeno0,c(genotype.ID,trait.ID)]
    
    ## Make sure there is only a single observation per genotype:
    if(length(which(duplicated(pheno[,genotype.ID]))) > 0){
      warning(paste("Duplicate phenotypes identified for",length(which(duplicated(pheno[,genotype.ID]))), "samples. Their mean values were used."))
      pheno <-
        do.call("rbind",
                sapply(unique(pheno[,genotype.ID]), function(gen){
                  temp.mean <- mean(pheno[pheno[,genotype.ID] == gen,trait.ID], na.rm = TRUE)
                  temp.out <- pheno[pheno[,genotype.ID] == gen,][1,,drop=FALSE]
                  temp.out[,trait.ID] <- round(temp.mean,1)
                  return(temp.out)}, simplify = FALSE))
    }
    
    if(length(which(is.na(pheno[,trait.ID]))) > 0){
      write(paste("A further",length(which(is.na(pheno[,trait.ID]))),
                  "matched individuals had missing phenotypic values and were removed.\n"),
            log.conn)
      pheno <- pheno[!is.na(pheno[,trait.ID]),]
    }
    
    pheno <- pheno[order(pheno[,1]),]
    phenoGeno0 <- pheno[,1]
  }
  
  colnames(pheno)[2] <- "Pheno"
  phenoGeno <- as.character(unique(pheno[,1]))
  popSize <- length(phenoGeno)
  
  if(popSize < ploidy + ploidy2) stop("Insufficient population size to fit QTL model with IBD probabilities. Suggest to use single marker approach")
  if(popSize <= 2*(ploidy + ploidy2)) warning("Population size may be too small for accurate model fitting. Results should be interpreted with caution.")
  
  
  lm0form <- "Pheno ~ 1"
  lm0Dat <- as.data.frame(pheno)
  
  if(blockTF) lm0form <- "Pheno ~ Block"
  
  ## Record the RSS from the null model, without QTL co-factors (but with blocks, if they are present)
  lmNULL <- lm(lm0form, data = lm0Dat)
  RSSNULL <- anova(lmNULL)[nrow(anova(lmNULL)),2]
  
  QTLmodel.list <- lapply(1:nrow(QTL_df), function(n) combn(1:nrow(QTL_df),n))
  
  PVE.list <- lapply(QTLmodel.list, function(q_mat = QTLmodel.list[[1]]){
    Qmodels <- setNames(apply(q_mat,2,function(sub_q){
      cof.mat <- do.call("cbind",lapply(sub_q, function(l) {
        lg <- QTL_df[l,1]
        temp.IBD <- IBD_list[[lg]]$IBDarray
        map <- IBD_list[[lg]]$map
        
        if(is.null(gap)){ #no splines were fitted..
          ## Try to find where the QTL should fit:
          minpos <- which.min(abs(QTL_df[l,2] - map$position))
          tempIBD <- t(temp.IBD[minpos,,phenoGeno])
          
          if(min(abs(QTL_df[l,2] - map$position)) > 0)
            
            warning(paste(QTL_df[l,2],"does not specify the co-factor position precisely! Proceeding with nearest alternative position of", map$position[minpos]))
          
        } else{
          tempIBD <- t(temp.IBD[paste0("cM",as.character(QTL_df[l,2])),,phenoGeno])
        }
        
        if(IBDtype == "haplotypeIBD"){
          out <- tempIBD
        } else{
          out <- tempIBD %*% indicatorMatrix
        }
        ## Exert boundary condition:
        out <- out[,setdiff(1:(ploidy + ploidy2),c(1,ploidy + 1))]
        
        ## Need to concatenate if there are blocks:
        out.bl <- do.call(rbind, lapply(1:nr.block, function(x) out))
        
        return(out.bl)
      }))
      
      colnames(cof.mat) <- paste0("C",1:ncol(cof.mat)) #Give them colnames so we can use them in model.
      if(any(is.na(cof.mat))) cof.mat[is.na(cof.mat)] <- 0
      
      
      lm0form <- as.formula(paste(c(lm0form, paste(colnames(cof.mat), collapse = " + ")), collapse = " + "))
      lm0Dat <- as.data.frame(cbind(pheno,cof.mat))
      
      lm0Res <- lm(lm0form, data = lm0Dat)
      RSS0 <- anova(lm0Res)[nrow(anova(lm0Res)),2]
      
      # Return proportion of explained variance of QTL sub-model
      return(100*round(1 - RSS0/RSSNULL,4))
    }),
    apply(q_mat,2,function(x) paste0(paste0("Q",x), collapse = "_")))
  })
  
  names(PVE.list) <- paste(seq(length(PVE.list)),"QTL model")
  
  return(PVE.list)
  
} #PVE


#' General QTL function that allows for co-factors, completely randomised block designs and the possibility to derive LOD thresholds using a permutation test
#' @description Function to run QTL analysis using IBD probabilties given (possibly replicated) phenotypes, assuming randomised experimental design
#' @param IBD_list List of IBD probabilities
#' @param Phenotype.df A data.frame containing phenotypic values
#' @param genotype.ID The colname of \code{Phenotype.df} that contains the offspring identifiers (F1 names)
#' @param trait.ID The colname of \code{Phenotype.df} that contains the response variable to use in the model
#' @param block The blocking factor to be used, if any (must be colname of \code{Phenotype.df}). By default \code{NULL}, in which case no blocking structure (for unreplicated experiments)
#' @param cofactor_df A 3-column data frame of co-factor(s); column 1 gives the numeric linkage group identifier(s), 
#' column 2 specifies the cM position of the co-factor(s), column 3 specifies whether the QTL was fitted using "a" = additive effects or 
#' "f" = full allelic interactions (note that any other symbol for the full model will also be accepted, as long as it is not "a"). 
#' For backward compatibility with package versions <= 0.0.9, it is possible to just supply the first two columns, 
#' in which case an additive-effects model is assumed for each cofactor (so, a third column will be automatically filled with "a"). 
#' By default \code{cofactor_df = NULL}, in which case no co-factors are included in the analysis.
#' @param allelic_interaction The QTL detection model can be for additive main effects only (by default \code{allelic_interaction = FALSE}). If \code{TRUE}, then the full model is used
#' (i.e. all possible genotype combinations are included as predictors in the model). This runs the risk of overfitting, especially if double reduction was also allowed. 
#' Both types of analyses can ideally be performed and compared. Note that if IBD probabilities were estimated using the "heuristic" method rather than the HMM method 
#' (see \code{\link{estimate_IBD}}), then IBDs are actually haplotype probabilities rather than genotype probabilities, meaning that allelic interaction effects cannot be included in the model.
#' @param folder If markers are to be used as co-factors, the path to the folder in which the imported IBD probabilities is contained can be provided here. 
#' By default this is \code{NULL}, if files are in working directory.
#' @param filename.short If TetraOrigin was used and co-factors are being included, the shortened stem of the filename of the \code{.csv} files containing the output of TetraOrigin, 
#' i.e. without the tail "_LinkageGroupX_Summary.csv" which is added by default to all output of TetraOrigin.
#' @param prop_Pheno_rep The minimum proportion of phenotypes represented across blocks. If less than this, the individual is removed from the analysis. If there is incomplete
#' data, the missing phenotypes are imputed using the mean values across the recorded observations.
#' @param perm_test Logical, by default \code{FALSE}. If \code{TRUE}, a permutation test will be performed to determine a
#' genome-wide significance threshold.
#' @param N_perm.max The maximum number of permutations to run if \code{perm_test} is \code{TRUE}; by default this is 1000.
#' @param alpha The P-value to be used in the selection of a threshold if \code{perm_test} is \code{TRUE}, by default 0.05 (i.e. the 0.95 quantile).
#' @param gamma The width of the confidence intervals used around the permutation test threshold using the approach of Nettleton & Doerge (2000), by default 0.05.
#' @param ncores Number of cores to use if parallel computing is required. Works both for Windows and UNIX (using \code{doParallel}).
#' Use \code{parallel::detectCores()} to find out how many cores you have available.
#' @param verbose Logical, by default \code{TRUE}. Should messages be printed during running?
#' @param log Character string specifying the log filename to which standard output should be written. If \code{NULL} log is send to stdout.
#' @param \dots Arguments passed to \code{\link{plot}}
#' @return
#' A nested list; each list element (per linkage group) contains the following items:
#' \describe{
#' \item{QTL.res}{Single matrix of QTL results with columns chromosome, position, LOD, adj.r.squared and PVE (percentage variance explained).}
#' \item{Perm.res}{If \code{perm_test} = \code{FALSE}, this will be \code{NULL}.
#' Otherwise, Perm.res contains a list of the results of the permutation test, with list items
#' "quantile","threshold" and "scores". Quantile refers to which quantile of scores was used to determine the threshold.
#' Note that scores are each of the maximal LOD scores across the entire genome scan per permutation, thus returning a
#' genome-wide threshold rather than a chromosome-specific threshold. If the latter is preferred, restricting the
#' \code{IBD_list} to a single chromosome and re-running the permutation test will provide the desired threshold.}
#' \item{Residuals}{If a blocking factor or co-factors are used, this is the (named) vector of residuals used as input for the
#' QTL scan. Otherwise, this is the set of (raw) phenotypes used in the QTL scan.}
#' \item{Map}{Original map of genetic marker positions upon which the IBDs were based, most often used
#' for adding rug of marker positions to QTL plots.}
#' \item{LG_names}{Names of the linkage groups}
#' \item{allelic_interaction}{Whether argument \code{allelic_interaction} was \code{TRUE} or \code{FALSE} in the QTL scan}
#' }
#' @examples
#' data("IBD_4x","Phenotypes_4x")
#' qtl_LODs.4x <- QTLscan(IBD_list = IBD_4x,
#'                        Phenotype.df = Phenotypes_4x,
#'                        genotype.ID = "geno",
#'                        trait.ID = "pheno",
#'                        block = "year")
#' @export
QTLscan <- function(IBD_list,
                    Phenotype.df,
                    genotype.ID,
                    trait.ID,
                    block = NULL,
                    cofactor_df = NULL,
                    allelic_interaction = FALSE,
                    folder = NULL,
                    filename.short,
                    prop_Pheno_rep = 0.5,
                    perm_test = FALSE,
                    N_perm.max = 1000,
                    alpha = 0.05,
                    gamma = 0.05,
                    ncores = 1,
                    log = NULL,
                    verbose = TRUE,
                    ...){
  
  ## Run some initial checks:
  if(!trait.ID %in% colnames(Phenotype.df)) stop("Trait colname incorrectly specified.")
  if(!genotype.ID %in% colnames(Phenotype.df)) stop("Genotype colname incorrectly specified.")
  IBD_list <- test_IBD_list(IBD_list)
  
  ploidy <- IBD_list[[1]]$ploidy
  ploidy2 <- IBD_list[[1]]$ploidy2
  
  ## In the case of no blocking or no co-factor and permutation test required, pass to fast_permute function:
  if(perm_test && is.null(block) && is.null(cofactor_df)){
    
    return(fast_permute(IBD_list = IBD_list,
                        Phenotype.df = Phenotype.df,
                        genotype.ID = genotype.ID,
                        trait.ID = trait.ID,
                        ploidy = ploidy,
                        ploidy2 = ploidy2,
                        N_perm = N_perm.max,
                        alpha = alpha,
                        allelic_interaction = allelic_interaction,
                        ncores = ncores,
                        verbose = verbose,
                        log = log))
    
  } else {
    
    if(is.null(log)) {
      log.conn <- stdout()
    } else {
      matc <- match.call()
      write.logheader(matc, log)
      log.conn <- file(log, "a")
    }
    
    if(!is.null(block)) if(!block %in% colnames(Phenotype.df)) stop("Blocking factor incorrectly specified.")
    
    blockTF <- cofactTF <- FALSE
    if(!is.null(block)) blockTF <- TRUE
    if(!is.null(cofactor_df)) cofactTF <- TRUE
    
    IBDarray <- IBD_list[[1]]$IBDarray
    IBDtype <- IBD_list[[1]]$IBDtype
    bivalent_decoding <- IBD_list[[1]]$biv_dec
    gap <- IBD_list[[1]]$gap
    
    if(IBDtype == "haplotypeIBD" & allelic_interaction) stop("QTL model with allele interactions is not possible with IBD probabilities generated using the heuristic method. Please re-run estimate_IBD using method = 'HMM'")
    
    if(IBDtype == "genotypeIBD"){
      if(is.numeric(IBD_list[[1]]$genocodes)){
        ## We are using the output of TetraOrigin
        GenCodes <- t(sapply(IBD_list[[1]]$genocodes,function(x) as.numeric(strsplit(as.character(x),"")[[1]])))
        Nstates <- nrow(GenCodes)
        indicatorMatrix <- matrix(0,nrow=Nstates, ncol=ploidy + ploidy2)
        for(r in 1:nrow(indicatorMatrix)){
          for(h in 1:ncol(GenCodes)){
            indicatorMatrix[r,GenCodes[r,h]] <- indicatorMatrix[r,GenCodes[r,h]] + 1
          }
        }
        
      } else{
        ## We are using the output of estimate_IBD
        GenCodes <- t(sapply(sapply(IBD_list[[1]]$genocodes,function(x) strsplit(x,"")), function(y) sapply(y,function(z) which(letters == z))))
        Nstates <- nrow(GenCodes)
        indicatorMatrix <- matrix(0,nrow=Nstates, ncol=ploidy + ploidy2)
        for(r in 1:nrow(indicatorMatrix)){
          for(h in 1:ncol(GenCodes)){
            indicatorMatrix[r,GenCodes[r,h]] <- indicatorMatrix[r,GenCodes[r,h]] + 1
          }
        }
        
      }
      
    } else if(IBDtype == "haplotypeIBD"){
      Nstates <- ploidy + ploidy2 #this is not actually needed..
      indicatorMatrix <- NULL
    } else{
      stop("Unable to determine IBD type, please check IBD_list elements have a valid $IBDtype tag.")
    }
    
    ## Could just check for equality here also, modulo check was from older version...
    if(dim(IBDarray)[2] %% Nstates != 0) stop("Incompatible input detected..")
    
    ## First get the starting set - the ones that were phenotyped and genotyped:
    if(any(is.na(Phenotype.df[,genotype.ID]))){
      warning("Missing genotype.ID values detected. Removing and proceeding without...")
      Phenotype.df <- Phenotype.df[!is.na(Phenotype.df[,genotype.ID]),]
    }
    
    phenoGeno0 <- intersect(
      dimnames(IBDarray)[[3]],
      unique(Phenotype.df[,genotype.ID]))
    
    if(verbose) write(paste("There are",length(phenoGeno0),"individuals with matching identifiers between phenotypic and genotypic datasets.\n"),log.conn)
    
    nr.block <- 1 #there is always 1 block...
    
    if(blockTF){
      write(paste("Summary of trait:",trait.ID,"over blocks\n"),log.conn)
      summTab <- tapply(Phenotype.df[,trait.ID],Phenotype.df[,block],summary)
      
      ## If there are no Missing values for some blocks, add NA = 0:
      for(i in seq(length(summTab))){
        if(!"NA's" %in% names(summTab[[i]])) {
          summTab[[i]] <- setNames(c(summTab[[i]],0),c(names(summTab[[i]]),"NA's"))
        }
      }
      summTab.df <- do.call("rbind", lapply(seq(length(summTab)), function(l) summTab[[l]]))
      rownames(summTab.df) <- names(summTab)
      write(knitr::kable(summTab.df), file=log.conn)
      
      df.columns <- c(genotype.ID,trait.ID,block)
      
      ## Make sure there is only a single observation per genotype per block:
      Phenotype.df <- do.call("rbind", lapply(levels(as.factor(Phenotype.df[,block])), function(bl){
        do.call("rbind",
                sapply(unique(Phenotype.df[Phenotype.df[,block]==bl,genotype.ID]), function(gen){
                  temp.mean <- mean(Phenotype.df[Phenotype.df[,block]==bl & Phenotype.df[,genotype.ID] == gen,trait.ID], na.rm = TRUE)
                  temp.out <- Phenotype.df[Phenotype.df[,block]==bl & Phenotype.df[,genotype.ID] == gen,][1,,drop=FALSE]
                  temp.out[,trait.ID] <- round(temp.mean,1)
                  return(temp.out)}, simplify = FALSE))
      }))
      
      Phenotype.df[,block] <- as.factor(Phenotype.df[,block])
      
      ## Get the counts of observations of the offspring over blocks:
      count.fun <- function(PGeno){
        colSums(sapply(PGeno, function(gen) sapply(levels(as.factor(Phenotype.df[,block])), function(bl) {
          out <- FALSE
          if(gen %in% Phenotype.df[Phenotype.df[,block]==bl,genotype.ID]){
            if(!is.na(Phenotype.df[Phenotype.df[,block]==bl & Phenotype.df[,genotype.ID] == gen,trait.ID])) out <- TRUE
          }
          return(out)
        }
        )))
      } #count.fun()
      
      counts <- count.fun(phenoGeno0)
      block.lev <- levels(as.factor(Phenotype.df[,block]))
      nr.block <- length(block.lev)
      
      rem <- which(counts < nr.block*prop_Pheno_rep)
      
      if(length(rem)!=0){
        
        write(paste("\nThe following individuals had fewer than",nr.block*prop_Pheno_rep,"observations over",nr.block,"blocks and were removed from the analysis:\n"),log.conn)
        write(knitr::kable(as.data.frame(matrix(c(names(rem),rep(" ",4-length(rem)%%4)),ncol=4))),file = log.conn)
        phenoGeno1 <- setdiff(phenoGeno0,names(rem))
        
      } else{
        phenoGeno1 <- phenoGeno0
      }
      
      ## Next, impute for the counts that are left:
      counts <- count.fun(phenoGeno1)
      imputers <- which(counts != nr.block)
      
      if(length(imputers) != 0) {
        
        #First take out the complete data - will add at the end:
        pheno0 <- Phenotype.df[Phenotype.df[,genotype.ID] %in% setdiff(phenoGeno1,names(imputers)),df.columns]
        
        warning(paste0("Estimating block effects with ",
                       length(setdiff(phenoGeno1,names(imputers))),
                       " of the total population size of ",
                       length(phenoGeno1),"!"))
        
        #Estimate block effects:
        lmDat <- pheno0[,2:3]
        colnames(lmDat) <- c("Pheno","Block")
        lm.block <- lm(Pheno ~ Block, data = lmDat)
        block.effects <- c(0,lm.block$coefficients[2:length(lm.block$coefficients)]) #Add a 0 temporarily, will remove later if not needed..
        
        pheno_impute <- NULL
        
        for(gen in names(imputers)){
          temp <- Phenotype.df[Phenotype.df[,genotype.ID]==gen,df.columns]
          if(nrow(temp)!=nr.block){
            #A bit of overkill, but general for any number of missing rows (previously, only 1 missing row was allowed.)
            miss.bl <- setdiff(block.lev,temp[,3])
            temp <- rbind(temp,setNames(data.frame(
              matrix(c(rep(gen,length(miss.bl)),
                       rep(NA,length(miss.bl)),
                       miss.bl),ncol = 3)),names(temp)))
          }
          temp <- temp[order(temp[,3]),]
          levelNA <- which(is.na(temp[,2]))
          levelOK <- setdiff(seq(nr.block),levelNA)
          
          temp[levelNA,2] <- mean(as.numeric(temp[levelOK,2]) - block.effects[levelOK]) + block.effects[levelNA]
          rownames(temp) <- paste(temp[1,1],temp[,3], sep="_")
          pheno_impute <- rbind(pheno_impute,temp)
        }
        
        pheno <- rbind(pheno0,pheno_impute)
        pheno[,2] <- as.numeric(pheno[,2])
        
      } else{
        pheno <- Phenotype.df[Phenotype.df[,genotype.ID] %in% phenoGeno1,df.columns]
        pheno[,2] <- as.numeric(pheno[,2])
      }
      
      pheno <- pheno[order(pheno[,3],pheno[,1]),]
      colnames(pheno)[3] <- "Block"
    } else{ #No blocks included.
      pheno <- Phenotype.df[Phenotype.df[,genotype.ID] %in% phenoGeno0,c(genotype.ID,trait.ID)]
      
      ## Make sure there is only a single observation per genotype:
      if(length(which(duplicated(pheno[,genotype.ID]))) > 0){
        warning(paste("Duplicate phenotypes identified for",length(which(duplicated(pheno[,genotype.ID]))), "samples. Their mean values were used."))
        pheno <-
          do.call("rbind",
                  sapply(unique(pheno[,genotype.ID]), function(gen){
                    temp.mean <- mean(pheno[pheno[,genotype.ID] == gen,trait.ID], na.rm = TRUE)
                    temp.out <- pheno[pheno[,genotype.ID] == gen,][1,,drop=FALSE]
                    temp.out[,trait.ID] <- round(temp.mean,1)
                    return(temp.out)}, simplify = FALSE))
      }
      
      if(length(which(is.na(pheno[,trait.ID]))) > 0){
        write(paste("A further",length(which(is.na(pheno[,trait.ID]))),
                    "matched individuals had missing phenotypic values and were removed.\n"),
              log.conn)
        pheno <- pheno[!is.na(pheno[,trait.ID]),]
      }
      
      pheno <- pheno[order(pheno[,1]),]
      phenoGeno0 <- pheno[,1]
    }
    
    colnames(pheno)[2] <- "Pheno"
    phenoGeno <- as.character(unique(pheno[,1]))
    popSize <- length(phenoGeno)
    
    if(popSize < ploidy + ploidy2) stop("Insufficient population size to fit QTL model with IBD probabilities. Suggest to use single marker approach")
    if(popSize <= 2*(ploidy + ploidy2)) warning("Population size may be too small for accurate model fitting. Results should be interpreted with caution.")
    
    ## Now deal with whether there is a (genetic) co-factor being included
    if(cofactTF){
      if(ncol(cofactor_df) == 2){
        cofactor_df[,3] <- "a" #for backward compatibility with v. <= 0.0.9
      }
      
      cof.mat <- do.call("cbind",lapply(1:nrow(cofactor_df), function(l) {
        
        lg <- cofactor_df[l,1]
        
        temp.IBD <- IBD_list[[lg]]$IBDarray
        map <- IBD_list[[lg]]$map
        
        if(is.null(gap)){
          if(!as.character(cofactor_df[l,2]) %in% dimnames(temp.IBD)[[1]]){
            ## Try to find where the co-factor should fit:
            minpos <- which.min(abs(cofactor_df[l,2] - map$position))
            tempIBD <- t(temp.IBD[minpos,,phenoGeno])
            
            if(min(abs(cofactor_df[l,2] - map$position)) > 0)
              warning(paste(as.character(cofactor_df[l,2]),
                            "does not specify the co-factor position precisely! Proceeding with nearest alternative position of", map$position[minpos]))
            
          } else{
            tempIBD <- t(temp.IBD[as.character(cofactor_df[l,2]),,phenoGeno])
          }
        } else{
          tempIBD <- t(temp.IBD[paste0("cM",as.character(cofactor_df[l,2])),,phenoGeno])
        }
        
        
        if(IBDtype == "haplotypeIBD"){
          out <- tempIBD
        } else{
          if(cofactor_df[l,3] == "a"){
            out <- tempIBD %*% indicatorMatrix
          } else{
            out <- tempIBD #this co-factor was detected with the full QTL model
          }
        }
        
        ## Need to concatenate if there are blocks:
        out.bl <- do.call(rbind, lapply(1:nr.block, function(x) out))
        
        return(out.bl)
      })
      )
      colnames(cof.mat) <- paste0("C",1:ncol(cof.mat)) #Give them colnames so we can use them in model.
      
      if(any(is.na(cof.mat))) cof.mat[is.na(cof.mat)] <- 0
    }
    
    ## The estimation of RSS0 depends on the model being used...
    lm0form <- "Pheno ~ 1"
    lm0Dat <- as.data.frame(pheno)
    
    if(blockTF) lm0form <- "Pheno ~ Block"
    
    if(cofactTF) {
      lm0form <- paste(c(lm0form, paste(colnames(cof.mat), collapse = " + ")), collapse = " + ")
      lm0Dat <- as.data.frame(cbind(pheno,cof.mat))
    }
    
    lm0form <- as.formula(lm0form)
    lm0Res <- lm(lm0form, data = lm0Dat)
    RSS0 <- anova(lm0Res)[nrow(anova(lm0Res)),2]
    
    ## Save the residuals from this model for further analysis:
    resid.pheno <- lm0Res$residuals
    
    if(verbose) write(paste("\nRunning QTL analysis with",popSize,"individuals...\n"),log.conn)
    
    if(length(IBD_list) > 1){
      IBDprobs <- abind::abind(lapply(seq(length(IBD_list)), function(i) IBD_list[[i]]$IBDarray[,,phenoGeno]),
                               along = 1)
    } else{
      IBDprobs <- IBD_list[[1]]$IBDarray[,,phenoGeno]
    }
    
    ## Compile the original map positions
    map <- do.call("rbind",
                   lapply(seq(length(IBD_list)), function(c)
                     IBD_list[[c]]$map[,c("chromosome","position")]
                   ))
    ## Compile the map positions of the IBDd positions also:
    if(!is.null(gap)) {map.1 <- do.call(rbind,
                                        lapply(seq(length(IBD_list)), function(c)
                                          data.frame("chromosome" = c,
                                                     "position" = as.numeric(setdiff(unlist(strsplit(dimnames(IBD_list[[c]]$IBDarray)[[1]],"M")),"c"))
                                          )#Just record chm and position.
                                        ))
    } else{
      map.1 <- map
    }
    
    run_lm <- function(IBDmat = IBDprobs[1,,],
                       residual_vect = resid.pheno,
                       IBDtyp = IBDtype){ #wghts is a column vector of IBDprobs
      
      if(IBDtyp == "genotypeIBD"){
        if(any(is.na(IBDmat))){
          warning(paste(round(100*length(which(is.na(IBDmat)))/length(IBDmat)),
                        "% of IBDs at supplied peak were NA. Proceeding without..."))
          IBDmat[is.na(IBDmat)] <- 0
        }
        if(!allelic_interaction){
          temp.haplo <- t(IBDmat)%*%indicatorMatrix
        } else{
          temp.haplo <- t(IBDmat)
        }
        
      } else{
        temp.haplo <- t(IBDmat)
      }
      
      # colnames(temp.haplo) <- paste0("X", 1:(ploidy + ploidy2))
      colnames(temp.haplo) <- paste0("X", 1:ncol(temp.haplo)) #this is more general
      
      if(blockTF) temp.haplo <- do.call("rbind", lapply(1:nr.block, function(n) temp.haplo))
      
      if(!allelic_interaction){
        modelterms <- paste0("X", setdiff(1:(ploidy + ploidy2),c(1,ploidy + 1))) #boundary condition applied
      } else{
        modelterms <- colnames(temp.haplo)
      }
      
      # modelterms <- paste0("X", setdiff(1:(ploidy + ploidy2),c(1,ploidy + 1))) #boundary condition applied
      lmDat <- as.data.frame(cbind(matrix(residual_vect,ncol=1),temp.haplo))
      colnames(lmDat)[1] <- "Resid.pheno"
      lmform <- as.formula(paste("Resid.pheno ~ ", paste(modelterms, collapse = " + "), collapse = ""))
      
      lmRes <- lm(formula = lmform, data = lmDat)
      
      anv <- anova(lmRes)
      RSS1 <- anv[nrow(anv),2]
      
      if(RSS1 == 0) warning("Model over-fitting encountered. No residual variation left on which to judge model fit! Perhaps increase population size?")
      LOD <- (popSize/2)*log10(RSS0/RSS1)
      adj.r.squared <- max(summary(lmRes)$adj.r.squared,0) #This can be negative..don't want negative R2adj, so set to 0
      PVE <- 100*(1 - 10^(-2*LOD/popSize))
      return(c(LOD,adj.r.squared,PVE))
    } #run_lm
    
    # Apply the regression to all the positions fitted using IBDs:
    LOD_var <- t(apply(IBDprobs,1, run_lm, resid.pheno, IBDtype))
    suppressWarnings(output <- cbind(map.1,LOD_var))
    colnames(output)[3:5] <- c("LOD","adj.r.squared","PVE")
    rownames(output) <- NULL
    
    if(perm_test){
      if(verbose) message("Running permutation test....")
      time_start <- Sys.time()
      
      Nran <- Nmin <- ceiling(5/alpha) #The starting number of permutations. Nran is placeholder for number of perm ran
      
      runPermutations <- function(N, IBDs, blockTF,
                                  ...){ #dots are passed to run_lm.perm
        
        run_lm.perm <- function(IBDmat,
                                residual_vect,
                                IBDtyp,
                                XMatrix,
                                pl,
                                pl2,
                                popN,
                                RSSzero,
                                ai,
                                numblock = 1){
          
          if(IBDtyp == "genotypeIBD"){
            if(!ai){
              temp.haplo <- t(IBDmat)%*%XMatrix
            } else{
              temp.haplo <- t(IBDmat)
            }
            
          } else{
            temp.haplo <- t(IBDmat)
          }
          # colnames(temp.haplo) <- paste0("X", 1:(pl + pl2))
          colnames(temp.haplo) <- paste0("X", 1:ncol(temp.haplo))
          
          if(blockTF) temp.haplo <- do.call("rbind", lapply(1:numblock, function(n) temp.haplo))
          
          if(!ai){
            modelterms <- paste0("X", setdiff(1:(pl + pl2),c(1,pl + 1)))
          } else{
            modelterms <- colnames(temp.haplo)
          }
          lmDat <- as.data.frame(cbind(matrix(residual_vect,ncol=1),temp.haplo))
          colnames(lmDat)[1] <- "Resid.pheno"
          lmform <- as.formula(paste("Resid.pheno ~ ", paste(modelterms, collapse = " + "), collapse = ""))
          
          lmRes <- lm(formula = lmform, data = lmDat)
          
          anv <- anova(lmRes)
          RSS1 <- anv[nrow(anv),2]
          LOD <- (popN/2)*log10(RSSzero/RSS1)
          adj.r.squared <- summary(lmRes)$adj.r.squared
          PVE <- 100*(1 - 10^(-2*LOD/popN))
          return(c(LOD,adj.r.squared,PVE))
        } #run_lm.perm
        
        if(blockTF){
          ## Need to permute within blocks....
          seqs <- sapply(1:nr.block, function(n) ((n-1)*popSize + 1):(n*popSize))
          Pheno.array <- lapply(1:N, function(n) {
            shufTem <- sample(seqs[,1])
            resid.pheno[c(shufTem,apply(seqs[,2:nr.block,drop=FALSE],2,function(x) x[shufTem]))]
          })
          
        } else{
          Pheno.array <- lapply(1:N, function(n) {
            resid.pheno[sample(1:length(resid.pheno))]
          })
        }
        
        if(ncores > 1){
          win <- Sys.info()["sysname"] == "Windows"
          if (win) {
            cl <- parallel::makeCluster(ncores)
            doParallel::registerDoParallel(cl)
          } else {
            doParallel::registerDoParallel(cores = ncores)
          }
          
          perm_list <-  split(seq(N), ceiling(round(seq(N)/(N/ncores),4)))
          
          LOD_vector <-
            foreach::foreach(
              i = 1:length(perm_list),
              .combine = c,
              .inorder = FALSE
            ) %dopar% {
              do.call("c",
                      lapply(perm_list[[i]], function(n)
                        max(apply(IBDs,1,run_lm.perm,Pheno.array[[n]],...)[1,]) #replaced IBDtype with ...
                      )
              )
            }
          
          if(win) parallel::stopCluster(cl)
          
        } else{ #If there are problems with parallel computing, offer an alternative.
          LOD_vector <- sapply(1:N, function(n) max(apply(IBDs,1,run_lm.perm,Pheno.array[[n]],...)[1,])) #replaced IBDtype with ...
        }
        
        return(LOD_vector)
      } #runPermutations
      
      perm.LODs <- runPermutations(N = Nmin,IBDs = IBDprobs,
                                   blockTF = blockTF,
                                   IBDtyp = IBDtype,XMatrix = indicatorMatrix,
                                   pl = ploidy, pl2 = ploidy2,
                                   popN = popSize, RSSzero = RSS0,
                                   numblock = nr.block,
                                   ai = allelic_interaction)
      
      if(any(is.infinite(perm.LODs))){
        warning("Unable to run permutation tests!")
        thresholds <- NA
      } else{
        thresholds <- sort(perm.LODs)[NettletonDoerge(N=Nmin, alpha = alpha, gamma = gamma)]
        
        ## Check whether the QTL position has been resolved, if not run more tests:
        while(any(LOD_var[,1] > thresholds[1] & LOD_var[,1] < thresholds[2]) & Nran < N_perm.max){
          ## Here we are unresolved => repeat the procedure, adding N_perm.min permutations each time
          message(paste0("There are still ",
                         length(which(LOD_var[,1] > thresholds[1] & LOD_var[,1] < thresholds[2])),
                         " unresolved LOD scores (from ",nrow(LOD_var),
                         "). Running a further ",Nmin," permutations to attempt to resolve these..."))
          perm.LODs <- c(perm.LODs, runPermutations(N = Nmin,IBDs = IBDprobs,
                                                    blockTF = blockTF, 
                                                    IBDtyp = IBDtype,XMatrix = indicatorMatrix,
                                                    pl = ploidy, pl2 = ploidy2,
                                                    popN = popSize, RSSzero = RSS0,
                                                    numblock = nr.block,
                                                    ai = allelic_interaction))
          thresholds <- sort(perm.LODs)[NettletonDoerge(N=length(perm.LODs),alpha = alpha, gamma = gamma)]
          Nran <- Nran + Nmin
        }
        
        if(all(LOD_var[,1] > thresholds[1] & LOD_var[,1] < thresholds[2])){
          message("Finally, all loci have a resolved LOD score!")
        } else{
          message(paste("Warning, there were",length(which(LOD_var[,1] > thresholds[1] & LOD_var[,1] < thresholds[2])),
                        "unresolved LOD scores after",Nran,"runs. It may be an idea to increase N_perm.max..."))
        }
        
      }
      
      time_end <- Sys.time()
      timediff <- as.numeric(time_end - time_start, units = "mins")
      
      write(
        paste(Nran,"permutations on a",Sys.info()["sysname"],"machine using",ncores,
              ifelse(ncores>1,"cores","core"),"took",round(timediff, 2),"minutes."),
        log.conn
      )
    }
    
    if(!is.null(log)) close(log.conn)
    
    out.resid <- pheno
    out.resid$Pheno <- resid.pheno
    
    if(!perm_test){
      return(list(QTL.res = output,
                  Perm.res = NULL,
                  Residuals = out.resid,
                  Map = map,
                  LG_names = names(IBD_list),
                  allelic_interaction = allelic_interaction))
    } else{
      return(list(QTL.res = output,
                  Perm.res = list(quantile = 1 - alpha,
                                  threshold = mean(thresholds), #take the mean as the threshold
                                  threshold.bounds = thresholds,
                                  scores = perm.LODs,
                                  Nperm = Nran),
                  Residuals = out.resid,
                  Map = map,
                  LG_names = names(IBD_list),
                  allelic_interaction = allelic_interaction))
    }
    
  }
} #QTLscan

#' Create a list of possible QTL segregation types
#' @description  Function to generate list of segregation types for the \code{\link{exploreQTL}} function
#' @param ploidy The ploidy of the population. Currently assumed to be an even number for this function.
#' @param segtypes List of QTL segregation types to consider,
#' so e.g. c(1,0) would mean all possible simplex x nulliplex QTL (ie. 4 QTL, on each of homologues 1 - 4 of parent 1).
#' Note that symmetrical QTL types that cannot be distinguished are not automatically removed and need to be manually identified.
#' If this is an issue, use the inbuilt list for tetraploids provided with the package to search the full model space.
#' Such an inbuilt list is currently only available for tetraploids, and is available from the \code{\link{exploreQTL}} function.
#' @param modes Character vector of modes of QTL action to consider, with options "a" for "additive" and "d" for dominant QTL action.
#' @export
segMaker <- function(ploidy,
                     segtypes,
                     modes=c("a","d")){
  
  modes <- match.arg(modes,several.ok = TRUE)
  len <- length(modes)*sum(sapply(seq(length(segtypes)),
                                  function(l)
                                    choose(ploidy,segtypes[[l]][1])*choose(ploidy,segtypes[[l]][2])))
  
  expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
  
  outlist <- lapply(1:len, function(x) NULL)
  counter <- 1
  
  for(l in seq(length(segtypes))){
    
    if(segtypes[[l]][1] %in% c(0,ploidy)){
      outcombs <- t(combn((ploidy + 1):(2*ploidy),segtypes[[l]][2]))
    } else if (segtypes[[l]][2] %in% c(0,ploidy)){
      outcombs <- t(combn(1:ploidy,segtypes[[l]][1]))
    } else{
      outcombs <- expand.grid.df(t(combn(1:ploidy,segtypes[[l]][1])),
                                 t(combn((ploidy + 1):(2*ploidy),segtypes[[l]][2])))
    }
    
    for(md in modes){
      for(r in 1:nrow(outcombs)){
        outlist[[counter]] <- list("homs"=outcombs[r,],"mode"=md)
        counter <- counter + 1
      }
    }
  }
  return(outlist)
} #segMaker

#' Run a single marker regression using marker dosages
#' @description Function to run a single marker regression using marker dosages
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns.
#' All markers in this matrix will be tested for association with the trait.
#' @param Phenotype.df A data.frame containing phenotypic values
#' @param genotype.ID The colname of \code{Phenotype.df} that contains the population identifiers (F1 names) (must be a colname of \code{Phenotype.df})
#' @param trait.ID The colname of \code{Phenotype.df} that contains the response variable to use in the model (must be a colname of \code{Phenotype.df})
#' @param maplist Option to include linkage map in the format returned by \code{MDSMap_from_list} from \code{polymapR}. If maplist is
#' not specified (by default \code{NULL}) then no ordering of markers from dosage-matrix is performed. Note that all markers in dosage_matrix are tested;
#' markers with dosages that were not on the maplist will be assigned unordered to linkage group 0 with dummy cM positions 1,2,3 etc.
#' @param perm_test Logical, by default \code{FALSE}. If \code{TRUE}, a permutation test will be performed to determine a
#' genome-wide significance threshold.
#' @param N_perm Integer. The number of permutations to run if \code{perm_test} is \code{TRUE}; by default this is 1000.
#' @param alpha Numeric. The P-value to be used in the selection of a threshold if \code{perm_test} is \code{TRUE};
#'  by default 0.05 (i.e. the 0.95 quantile).
#' @param ncores Number of cores to use if parallel processing required. Works both for Windows and UNIX (using \code{doParallel}).
#' Use \code{parallel::detectCores()} to find out how many cores you have available.
#' @param return_R2 Should the (adjusted) R2 of the model fit also be determined? 
#' @param log Character string specifying the log filename to which standard output should be written. If \code{NULL} log is send to stdout.
#' @return A list containing the following components:
#' \describe{
#' \item{QTL.res}{ The -log(p) of the model fit per marker are returned as "LOD" scores, although "LOP" would have been a better description.
#' If requested, R2 values are also returned in column "R2adj"}
#' \item{Perm.res}{ The results of the permutation test if performed, otherwise \code{NULL}}
#' \item{Map}{ The linkage map if provided, otherwise \code{NULL}}
#' \item{LG_names}{Names of the linkage groups, if a map was provided, otherwise \code{NULL}}
#' }
#' @examples
#' data("SNP_dosages.4x","BLUEs.pheno")
#' Trait_1.smr <- singleMarkerRegression(dosage_matrix = SNP_dosages.4x,
#' Phenotype.df = BLUEs.pheno,genotype.ID = "Geno",trait.ID = "BLUE")
#' @export
singleMarkerRegression <- function(dosage_matrix,
                                   Phenotype.df,
                                   genotype.ID,
                                   trait.ID,
                                   maplist = NULL,
                                   perm_test = FALSE,
                                   N_perm = 1000,
                                   alpha = 0.05,
                                   ncores = 1,
                                   return_R2 = FALSE,
                                   log = NULL) {
  
  dosage_matrix <- test_dosage_matrix(dosage_matrix)
  
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  if(any(is.na(Phenotype.df[,genotype.ID]))){
    warning("Missing genotype.ID values detected. Removing and proceeding without...")
    Phenotype.df <- Phenotype.df[!is.na(Phenotype.df[,genotype.ID]),]
  }
  
  phenoGeno <- intersect(colnames(dosage_matrix),unique(Phenotype.df[,genotype.ID]))
  
  if(length(phenoGeno) != nrow(Phenotype.df))
    warning(paste("Of the",nrow(Phenotype.df),"phenotypes, proceeding with core set of",
                  length(phenoGeno),"values. Screening for NA values..."))
  
  Phenotype.df <- Phenotype.df[match(phenoGeno,Phenotype.df[,genotype.ID]),]
  
  ## Single marker analysis cannot handle missing phenotypic data, remove these offspring as well:
  if(any(is.na(Phenotype.df[,trait.ID]))) phenoGeno <- phenoGeno[-which(is.na(Phenotype.df[,trait.ID]))]
  
  PhenotypeData <- Phenotype.df[match(phenoGeno,Phenotype.df[,genotype.ID]),trait.ID]
  GenotypeData <- dosage_matrix[,match(phenoGeno,colnames(dosage_matrix))]
  
  # popNum <- length(phenoGeno)
  # RSS0 <- sum((PhenotypeData - mean(PhenotypeData))^2)
  
  # run_lm.vect <- function(X,phenos){
  #   # X is a vector of genotype (dosage) data
  #   # phenos is a matrix of phenotypes, with at least 1 column
  #   lmRes <- lm(phenos ~ X,na.action = na.exclude)
  #   temppopN <- length(which(!is.na(X)))
  #   RSS1 <- apply(as.matrix(lmRes$residuals),2,function(y) sum(y^2))
  #   
  #   LOD <- ifelse(RSS1 == 0,0,(temppopN/2)*log10(RSS0/RSS1))
  #   return(LOD)
  # }
  
  run_lm.vect <- function(X,phenos){
    # X is a vector of genotype (dosage) data
    # phenos is a matrix of phenotypes, with at least 1 column
    aovRes <- aov(phenos ~ X,na.action = na.exclude)
    
    SS_res <- apply(as.matrix(aovRes$residuals),2,function(y) sum(y^2))
    msRes <- SS_res/aovRes$df.residual
    # Only 1 df in model, so SS_regression = MS_regression:
    msReg <- apply(as.matrix(aovRes$model$phenos),2,var)*(nrow(aovRes$model) - 1) - SS_res
    TT <- msReg/msRes
    
    p <- pf(TT, 1, aovRes$df.residual, lower.tail = FALSE)
    
    return(p)
  }
  
  # run_lm.vect <- function(X, phenos) {
  #   lmRes <- lm(phenos ~ X, na.action = na.exclude, singular.ok=T)
  #   p <- tryCatch(anova(lmRes)$'Pr(>F)'[1], error=function(e) NA)
  #   return(p)
  # }
  
  determine_R2 <- function(X,phenos){
    lmRes <- lm(phenos ~ X,na.action = na.exclude)
    R2 <- max(summary(lmRes)$adj.r.squared,0)
    return(R2)
  }
  
  ## Check if parallel processing required:
  if(ncores > 1){
    
    win <- Sys.info()["sysname"] == "Windows"
    if (win) {
      cl <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
    } else {
      doParallel::registerDoParallel(cores = ncores)
    }
    
    # perm.ls <- lapply(1:ncores, function(x) matrix(PhenotypeData,ncol=1))
    pheno.mat <- matrix(PhenotypeData,ncol=1)
    
    row.split <- split(1:nrow(GenotypeData), ceiling(round((1:nrow(GenotypeData))/(nrow(GenotypeData)/ncores),4)))
    
    if(perm_test){
      perms <- sapply(1:N_perm, function(n) sample(PhenotypeData))
      pheno.mat <- cbind(pheno.mat,perms)
    }
    
    ## R CMD CHECK work-around:
    i <- NULL
    
    LOP.res <-
      foreach::foreach(
        i = 1:length(row.split),
        # .combine = c,
        .inorder = FALSE
      ) %dopar% {
        list(set = row.split[[i]],
             res = apply(GenotypeData[row.split[[i]],,drop=FALSE],1,
                         function(x) run_lm.vect(X=x, phenos = pheno.mat))
        )
      }
    
    if(win) parallel::stopCluster(cl)
    
    extracted.order <- do.call(c,lapply(LOP.res,function(x) x$set))
    
    if(perm_test){
      LOPResults <- -log10(do.call(c,lapply(LOP.res,function(x) x$res[1,])))
    } else{
      LOPResults <- -log10(do.call(c,lapply(LOP.res,function(x) x$res)))
    }
    
    #Re-order (just in case order lost in parallel call)
    LOPResults <- LOPResults[order(extracted.order)]
    
    if(perm_test) {
      permLOPs <- -log10(do.call(cbind, lapply(1:ncores, function(l) LOP.res[[l]]$res[-1,])))
      
      threshold <- quantile(apply(permLOPs,1,max), probs = 1 - alpha)
      
      out.list <- list("QTL.res" = LOPResults,
                       "Perm.res" = list(
                         quantile = 1 - alpha,
                         threshold = threshold,
                         threshold.bounds = c(threshold,threshold),#for now, no Nettleton Doerge approximation
                         scores = apply(permLOPs,1,max))
      )
      
    } else{
      out.list <- list("QTL.res" = LOPResults,
                       "Perm.res" = NULL)
    }
    
  } else{ #No parallel
    if(perm_test) stop("Permutation testing without ncores > 1 is not implemented. Please try again.")
    
    out.list <- list("QTL.res" = -log10(apply(GenotypeData,1,run_lm.vect,matrix(PhenotypeData,ncol=1))),
                     "Perm.res" = NULL)
  }
  
  ## Run R2 analysis separately:
  if(return_R2) R2res <- apply(GenotypeData,1,determine_R2,PhenotypeData)
  
  ## Order the markers if a maplist is provided
  if(!is.null(maplist)){
    
    mapmat <- matrix(c(rep(seq(length(maplist)),sapply(maplist,nrow)),
                       do.call(rbind,maplist)[,"position"]),
                     ncol=2,dimnames = list(do.call(rbind,maplist)[,"marker"],
                                            c("chromosome","position")))
    
    mapped_checked <- intersect(rownames(GenotypeData),rownames(mapmat))
    
    if(length(mapped_checked) < length(out.list$QTL.res)){
      warning("Unmapped markers detected in dosage_matrix were also analysed. Assigning to chromosome 0 with dummy positions...")
      mapmat <- rbind(matrix(c(rep(0,length(setdiff(names(out.list$QTL.res),mapped_checked))),
                               1:length(setdiff(names(out.list$QTL.res),mapped_checked))),
                             ncol=2,
                             dimnames=list(setdiff(names(out.list$QTL.res),mapped_checked),
                                           c("chromosome","position"))),
                      mapmat[mapped_checked,])
      
      mapped_checked <- intersect(rownames(GenotypeData),rownames(mapmat))
    }
    
    if(!return_R2){
      out.list$QTL.res <- cbind(mapmat[mapped_checked,],
                                out.list$QTL.res[match(mapped_checked,names(out.list$QTL.res))])
    } else{
      out.list$QTL.res <- cbind(mapmat[mapped_checked,],
                                out.list$QTL.res[match(mapped_checked,names(out.list$QTL.res))],
                                R2res[match(mapped_checked,names(R2res))])
    }
    
    out.list$QTL.res <- out.list$QTL.res[order(out.list$QTL.res[,"chromosome"],out.list$QTL.res[,"position"]),]
    
    if(return_R2){
      colnames(out.list$QTL.res)[3:4] <- c("LOD","R2adj")
    } else{
      colnames(out.list$QTL.res)[3] <- "LOD" 
    }
    
    out.list <- c(out.list,list("Map" = mapmat, "LG_names" = names(maplist)))
    
  } else{ #No maplist
    
    if(return_R2){
      out.list$QTL.res <- matrix(c(out.list$QTL.res,R2res[match(names(out.list$QTL.res),names(R2res))]),ncol=2,
                                 dimnames = list(names(out.list$QTL.res),c("LOD","R2adj")))
    } else{
      out.list$QTL.res <- matrix(out.list$QTL.res,ncol=1,
                                 dimnames = list(names(out.list$QTL.res),"LOD"))
    }
    
    out.list <- c(out.list,list("Map" = NULL, "LG_names" = NULL))
  }
  
  ## For compatibility with IBD output
  out.list$QTL.res <- as.data.frame(out.list$QTL.res)
  
  return(out.list)
} #singleMarkerRegression

#' Fit splines to IBD probabilities
#' @description Fits splines to IBD probabilities at a grid of positions at user-defined spacing.
#' @param IBD_list List of IBD probabilities
#' @param gap The size (in centiMorgans) of the gap between splined positions
#' @param method One of two options, either "linear" or "cubic". The default method (cubic) fits cubic splines, and although more accurate,
#' becomes computationally expensive in higher-density data-sets, where the linear option may be preferable.
#' @param ncores Number of cores to use, by default 1 only. Works both for Windows and UNIX (using \code{doParallel}).
#' Use \code{parallel::detectCores()} to find out how many cores you have available. Note that with large datasets, using multiple cores
#' will use large amounts of memory (RAM). Single-core or e.g. 2-core evaluations, although slower, is less memory-intensive.
#' @param log Character string specifying the log filename to which standard output should be written. If \code{NULL} log is send to stdout.
#' @return Returns a list of similar format as IBD_list, with a splined \code{IBD_array} in place of the original \code{IBD_array}
#' @examples
#' data("IBD_4x")
#' IBD_4x.spl <- spline_IBD(IBD_list = IBD_4x, gap = 1)
#' @export
spline_IBD <- function(IBD_list,
                       gap,
                       method = "cubic",
                       ncores = 1,
                       log = NULL){
  
  IBD_list <- test_IBD_list(IBD_list)
  
  if(!is.numeric(gap)) stop("gap must be numeric (e.g. gap = 1)")
  
  method <- match.arg(method,choices = c("linear","cubic"))
  
  if(is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  spline_LG <- function(l){
    map.pos <- IBD_list[[l]]$map$position
    gridPos <- mapseq(map.pos, gap = gap)
    
    ## Generate an empty array for the splined IBDs:
    IBDarray.spl <- array(data = 0,dim = c(length(gridPos),dim(IBD_list[[l]]$IBDarray)[2:3]),
                          dimnames = c(list(paste0("cM",gridPos)),dimnames(IBD_list[[l]]$IBDarray)[2:3]))
    
    predict.linear <- function(ibd, mp, grd){
      grd_snip <- grd[2:(length(grd)-1)] #this is because the end-points need not be interpolated
      hits <- as.factor(findInterval(grd_snip,mp))
      
      pred <- function(x0,x1,y0,y1,x_targ) {
        return(y0 + (y1-y0)*(x_targ - x0)/(x1-x0))
      }
      
      out <- do.call(c,lapply(as.numeric(levels(hits)),function(n){
        pred(x0 = mp[n], x1 = mp[n+1], 
             y0 = ibd[n], y1 = ibd[n+1],
             x_targ = grd_snip[hits==n])
      }))
      
      ## return first and last positions (that were snipped out)
      return(c(ibd[1],out,ibd[length(ibd)]))
    }
    
    ## Splining is per individual. Fill in a loop
    for(f1 in 1:dim(IBDarray.spl)[3]){
      IBDarray.spl[,,f1] <- apply(IBD_list[[l]]$IBDarray[,,f1],2,function(xm){
        if(method == "linear"){
          return(predict.linear(ibd = xm,mp = map.pos,grd = gridPos))
        } else{
          return(predict(smooth.spline(map.pos,xm),gridPos)$y)
        }
      })
      IBDarray.spl[,,f1][IBDarray.spl[,,f1] < 0] <- 0
      # normalise:
      IBDarray.spl[,,f1] <- IBDarray.spl[,,f1]/rowSums(IBDarray.spl[,,f1])
    }
    
    outlist <- IBD_list[[l]]
    outlist$IBDarray <- IBDarray.spl
    outlist$gap <- gap
    
    return(outlist)
  } #spline_LG
  
  if(ncores > 1){
    win <- Sys.info()["sysname"] == "Windows"
    if (win) {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
    } else {
      doParallel::registerDoParallel(cores = ncores)
    }
    
    ## R CMD CHECK work-around:
    lgp <- NULL
    
    IBD_out <- 
      foreach::foreach(
        lgp = seq(length(IBD_list)),
        .inorder = TRUE
      ) %dopar% {
        spline_LG(lgp)
      }
    
    if(win) parallel::stopCluster(cl)
    
  } else{ #no parallel processing
    IBD_out <- lapply(seq(length(IBD_list)), spline_LG)
  }
  
  if(!is.null(names(IBD_list))){
    names(IBD_out) <- names(IBD_list)
  } else{
    warning(paste("IBD_list had no LG names assigned. Assigning names",
                  paste0("'LG",seq(length(IBD_out)),"'", collapse = ", "),"to output list."))
    names(IBD_out) <- paste0("LG",seq(length(IBD_out)))
  }
  
  
  return(IBD_out)
}


#' Thin out map data 
#' @description \code{thinmap} is a function for thinning out an integrated map, in order that IBD estimation runs more quickly. Especially
#' useful for maps with very high marker densities for which the \code{\link{estimate_IBD}} function is to be used.
#' @param maplist A list of maps. In the first column marker names and in the second their position.
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns.
#' @param bin_size Numeric. Size (in cM) of the bins to include. By default, a bin size of 1 cM is used. Larger \code{bin_size} results in 
#' fewer markers being left on the resulting map. 
#' @param bounds Numeric vector. If \code{NULL} (by default) then all positions are included, however if specified then output
#' is limited to a specific region, which may be useful if fine-mapping a region of interest.
#' @param remove_markers Optional vector of marker names to remove from the maps. Default is \code{NULL}.
#' @param plot_maps Logical. Plot the marker positions of the selected markers using \code{polymapR::plot_map}.
#' @param use_SN_phase Logical, by default \code{FALSE}. If \code{TRUE}, then 1x0 and 0x1 are binned per phase, to increase coverage of these
#' marker types across parental homologues. If not, at most one of each are retained per bin.
#' @param parent1 Identifier of parent 1, by default assumed to be \code{"P1"}
#' @param parent2 Identifier of parent 2, by default assumed to be \code{"P2"}
#' @param log Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.
#' @return A maplist of the same structure as the input maplist, but with fewer markers based on the bin_size.
#' @examples
#' data("phased_maplist.4x","SNP_dosages.4x")
#' maplist_thin<-thinmap(maplist=phased_maplist.4x,dosage_matrix=SNP_dosages.4x)
#' @export
thinmap <- function(maplist,
                    dosage_matrix,
                    bin_size = 1,
                    bounds = NULL,
                    remove_markers = NULL,
                    plot_maps = TRUE,
                    use_SN_phase = FALSE,
                    parent1 = "P1",
                    parent2 = "P2",
                    log = NULL) {
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  dosage_matrix <- test_dosage_matrix(dosage_matrix)
  
  if(is.null(names(maplist))){
    if(length(maplist) == 1) stop("maplist provided was an unnamed list of length 1. Chromosome number not found.\n Assign names to maplist using names() function first!")
    warning("maplist provided was an unnamed list. Proceeding with automatic names")
    names(maplist) <- paste0("LG",seq(length(maplist)))
  }
  
  outlist <- lapply(seq(length(maplist)),function(mapn){
    
    map <- maplist[[mapn]]
    
    chrn <- mapn
    if(length(maplist) == 1) chrn <- gsub("(.*)([0-9]+)(.*)","\\2", names(maplist)[mapn])
    
    if(!is.null(remove_markers)) map <- map[!map$marker %in% remove_markers,]
    
    if(!is.null(bounds)){
      if(length(bounds) != 2 | bounds[2]<=bounds[1]) stop("Incorrect bounds specified. This should be a vector of two distinct cM positions.")
      map <- map[map[,"position"] >= bounds[1] & map[,"position"] <= bounds[2],]
    }
    
    bins <- seq(0,ceiling(max(map["position"])) + bin_size,by = bin_size)
    
    binned.markers <- lapply(1:(length(bins)-1), function(n) as.character(map[,"marker"])[map[,"position"] >= bins[n] & map[,"position"] < bins[n+1]])
    
    
    remove.marks <- function(merkers){
      numNA <- sapply(merkers, function(m) length(which(is.na(dosage_matrix[match(m,rownames(dosage_matrix)),3:ncol(dosage_matrix)]))))
      
      ParSeg <- as.factor(sapply(merkers, function(m) paste0(dosage_matrix[match(m,rownames(dosage_matrix)),1],
                                                             dosage_matrix[match(m,rownames(dosage_matrix)),2]))
      )
      
      whittled<-do.call(c,lapply(levels(ParSeg), function(ps) setdiff(names(ParSeg[ParSeg==ps]),names(which.min(numNA[ParSeg==ps])))
      ))
      
      return(whittled)
    }
    
    # if(!use_SN_phase){
    removed.markers <- do.call(c,lapply(seq(length(binned.markers)),function(l)
      if(length(binned.markers[[l]]) > 0) remove.marks(binned.markers[[l]])))
    
    # } else{
    
    # }
    
    
    map <- map[!map$marker%in%removed.markers,]
    
    if(length(which(map[,"marker"] %in% rownames(dosage_matrix))) != nrow(map)) {
      warning("Not all mapped markers have corresponding dosages. Attempting to continue with those that do...")
      map <- map[map[,"marker"] %in% rownames(dosage_matrix),]
    }
    
    if(nrow(map) < 2) stop("Insufficient map information to proceed. Check maplist and/or dosage_matrix.")
    
    write(paste0("\n",nrow(map)," markers from a possible ",
                 nrow(maplist[[mapn]])," on LG",mapn," were included."),
          log.conn)
    
    return(map)
  })
  
  if(!is.null(log)) close(log.conn)
  
  if(plot_maps){
    polymapR::plot_map(maplist=outlist,
                       color_by_type = TRUE,
                       dosage_matrix = dosage_matrix,
                       bg_col = "white",
                       parent1 = parent1,
                       parent2 = parent2)}
  
  if(is.null(names(outlist))){
    names(outlist) <- names(maplist)
  }
  
  return(outlist)
}

#' Visualise Genotypic Information Coefficient
#' @description Function to visualise the GIC of a certain region
#' @param GIC_list List of GIC data, the output of \code{\link{estimate_GIC}}
#' @param add_rug Should original marker positions be added to the plot?
#' @param add_leg Should a legend be added to the plot?
#' @param ylimits Optional argument to control the plotting area, by default \code{NULL}
#' @param gic.cex Option to increase the size of the GIC
#' @param show_markers Should markers be shown?
#' @param add.mainTitle Should a main title be added to the plot?
#' @param plot.cols Optional argument to specify plot colours, otherwise suitable contrasting colours are chosen
#' @return The phased map data for the specified region, recoded into 1's and 0's.
#' @examples
#' data("GIC_4x")
#' visualiseGIC(GIC_list = GIC_4x)
#' @export
visualiseGIC <- function(GIC_list,
                         add_rug = TRUE,
                         add_leg = FALSE,
                         ylimits = NULL,
                         gic.cex = 1,
                         show_markers = TRUE,
                         add.mainTitle=TRUE,
                         plot.cols = NULL){
  
  ploidy <- GIC_list[[1]]$ploidy
  ploidy2 <- GIC_list[[1]]$ploidy2
  
  if(ploidy == 2 & ploidy2 == 2) warning("GIC per homologue is identical within each parent! GIC lines overlap and only a single line is visible")
  
  #Currently can handle up to 12 colours.. If different ploidy levels, plotting colours will be different
  if(is.null(plot.cols)) {
    plot.cols <- RColorBrewer::brewer.pal(max(2*ploidy, 2*ploidy2),"Paired")[c(seq(2, 2*ploidy,2),seq(2,2*ploidy2,2))]
  } else{
    ## Make sure plot.cols is the expected length:
    plot.cols <- cbind(plot.cols, 1:(ploidy+ploidy2))[,1] #use R's vector recycling
  }
  cols.rgb <- t(col2rgb(plot.cols))/255
  
  for(j in seq(length(GIC_list))){
    
    GICdata <- GIC_list[[j]]
    
    if(!show_markers){
      
      if(is.null(ylimits)) ylimits <- c(0,1)
      
      plot(NULL,xlim=range(as.numeric(rownames(GICdata$GIC))),
           ylim=ylimits,xlab="cM",ylab="GIC",
           axes=FALSE,cex.lab=1.1,
           main = ifelse(add.mainTitle,paste("LG", GICdata$map$chromosome[1]),""))
      axis(1)
      axis(2,las=1)
      if(add_rug) rug(GICdata$map[,"position"])
      
      for(h in 1:(ploidy + ploidy2)) lines(rownames(GICdata$GIC),GICdata$GIC[,h], col = plot.cols[h], lwd=2)
      
      if(add_leg) legend(x = max(GICdata$map[,"position"])/4, y = 0.5,
                         lwd=2,col=plot.cols,
                         legend = paste0("h",1:(ploidy + ploidy2)),
                         title="Parental homologues",
                         bty="n",cex=0.9, ncol=2)
      
    } else{ #Show markers as well
      if(!all(rownames(GICdata$parental_phase) == GICdata$map$marker)) stop("Differences in integrated map and parental map detected.")
      
      round.down <- function(x, interval = 0.2) x - (x %% interval)
      
      ## I will split the parents for clarity:
      p1_gic1 <- seq(round.down(min(GICdata$GIC[,1:ploidy])),1,0.2)
      p2_gic1 <- seq(round.down(min(GICdata$GIC[,(ploidy + 1):(ploidy + ploidy2)])),1,0.2)
      
      yrange <- c(1, 1 + ploidy + ploidy2 + length(c(p1_gic1,p2_gic1)) + 0.25)
      ticposns <- c(1:ploidy,(ploidy+1):(ploidy+1+length(p1_gic1) - 1),
                    (ploidy+1+length(p1_gic1)+1):(ploidy+1+length(p1_gic1)+ploidy2),
                    (yrange[2] - 0.25 - length(p2_gic1) + 1):(yrange[2] - 0.25))
      
      p1transform <- function(x) {
        slope <- (length(p1_gic1)-1)/(p1_gic1[length(p1_gic1)] - p1_gic1[1])
        return(slope*(x - p1_gic1[1]) + ploidy + 1)
      }
      
      p2transform <- function(x) {
        slope <- (length(p2_gic1)-1)/(p2_gic1[length(p2_gic1)] - p2_gic1[1])
        return(slope*(x - 1) + ticposns[length(ticposns)])
      }
      
      plot(NULL, xlim=range(GICdata$map$position),ylim=yrange,
           xlab="cM", ylab = "", axes=FALSE,
           main = ifelse(add.mainTitle,paste("LG", GICdata$map$chromosome[1]),""))
      axis(1); Hmisc::minor.tick(nx=4,ny=1)
      axis(2, las=1, at = ticposns,
           labels=c(paste0("h",1:ploidy),p1_gic1,
                    paste0("h",(ploidy + 1):(ploidy + ploidy2)),p2_gic1))
      abline(h = ploidy + 1, col = "gray65")
      abline(h = setdiff(ticposns[1]:ticposns[length(ticposns)],ticposns))
      abline(h = (yrange[2] - 0.25 - length(p2_gic1) + 1), col = "gray65")
      box()
      
      ## Render P1 data first:
      for(i in 1:ploidy){
        hits <- which(GICdata$parental_phase[,i]==1)
        points(GICdata$map$position[hits],rep(i, length(hits)), pch = 19,
               col=rgb(red = cols.rgb[i,1],
                       green = cols.rgb[i,2],
                       blue = cols.rgb[i,3],
                       alpha = 0.5))
        lines(rownames(GICdata$GIC),p1transform(GICdata$GIC[,i]), col = plot.cols[i], lwd = 2)
      }
      ## Now render P2 data:
      for(i in (ploidy + 1):(ploidy + ploidy2)) {
        hits <- which(GICdata$parental_phase[,i]==1)
        points(GICdata$map$position[hits],rep(i+length(p1_gic1)+1, length(hits)), pch = 19,
               col=rgb(red = cols.rgb[i,1],
                       green = cols.rgb[i,2],
                       blue = cols.rgb[i,3],
                       alpha = 0.5))
        lines(rownames(GICdata$GIC),p2transform(GICdata$GIC[,i]), col = plot.cols[i], lwd = 2)
      }
      ## Add some text:
      text(x = par()$usr[1] + 0.15*(par()$usr[2] - par()$usr[1]), y = ploidy + 1.35,pos = 1,labels = "Parent 1", font = 3,cex=0.95)
      text(x = par()$usr[1] + 0.15*(par()$usr[2] - par()$usr[1]), y = (yrange[2] - length(p2_gic1) + 1.05),pos = 1,labels = "Parent 2", font = 3,cex=0.95)
      
      mtext(text = "GIC",side = 2,line = 3,
            at = c(mean((ploidy+1):(ploidy+length(p1_gic1))),
                   mean((yrange[2] - length(p2_gic1) + 0.75):(yrange[2] - 0.25))), cex=gic.cex)
      
    }
  }
} #visualiseGIC

#' Visualise haplotypes in certain individuals in a certain region
#' @description Function to visualise the haplotypes of a certain region in certain individuals
#' @param IBD_list List of IBD probabilities
#' @param display_by Option to display a subset of the population's haplotypes either by \code{"phenotype"} or \code{"name"}.
#' If \code{"phenotype"} is supplied, then \code{Phenotype.df},\code{genotype.ID},\code{trait.ID} and \code{pheno_range} must also be specified.
#' if \code{"name"} is supplied, then \code{select_offspring} must be specified.
#' @param linkage_group Numeric identifier of the linkage group being examined, based on the order of \code{IBD_list}.
#' Only a single linkage group is allowed. If \code{IBD_list} corresponds to a single linkage group, default value
#' of \code{NULL} will suffice
#' @param Phenotype.df A data.frame containing phenotypic values, which can be used to select a subset of the population
#' to visualise (with extreme phenotypes for example). By default \code{NULL}, in which case a subset of the
#' population may be selected using the \code{select_offspring} argument.
#' @param genotype.ID The colname of \code{Phenotype.df} that contains the population identifiers (F1 names) (must be a colname of \code{Phenotype.df})
#' @param trait.ID The colname of \code{Phenotype.df} that contains the response variable to use in the model (must be a colname of \code{Phenotype.df})
#' @param pheno_range Vector of numeric bounds of the phenotypic scores to include (offspring selection).
#' @param cM_range Vector of numeric bounds of the genetic region to be explored. If none are specified, the default of \code{"all"} means all cM positions will be included.
#' @param highlight_region Option to hightlight a particular genetic region on the plot; can be a single position or a vector of 2 positions. By default \code{NULL}.
#' @param select_offspring Vector of offspring identifiers to visualise, must be supplied if \code{display_by} = \code{"name"}. Specifying \code{"all"} will result in
#' all offspring haplotypes being visualised.
#' @param recombinant_scan Vector of homologue numbers between which to search for recombinant offspring in the visualised region and selected individuals.
#' By default \code{NULL}, in which case no search is preformed.
#' @param allele_fish Vector of homologue numbers of interest, for which to search for offspring that carry these homologues (in the visualised
#' region). By default \code{NULL}, in which case no search ("fishing") is performed.
#' @param presence_threshold Numeric. The minimum probability used to declare presence of a homologue in an individual. This is only needed if a \code{recombinant_scan} is
#' performed. By default a value of 0.95 is used. When searching for recombinants, this value is also used to denote the proportion of loci carrying the required number
#' of homologues (i.e. by default 95 per cent of loci should have between 0.95 and 1.1 copies of the specified recombinant homologues).
#' @param xlabl Logical, by default \code{TRUE}. Should an x-axis label be used?
#' @param ylabl Logical, by default \code{TRUE}. Should a y-axis label be used?
#' @param mainTitle Option to override default plot titles with a (vector of) captions. By default \code{NULL}.
#' @param multiplot Vector of integers. By default \code{NULL} so haplotypes are plotted singly;
#' otherwise a vector specifying the number of rows and columns in the plot layout.
#' @param append Option to allow user to append new plots to spaces generated by \code{multiplot}, otherwise these are
#' filled with blank plots. By default \code{FALSE}. If \code{TRUE}, then a large enough \code{multiplot} grid should be
#' generated to make this option meaningful.
#' @param colPal Colour palette to use in the visualisation (best to provide 3 colours).
#' @param hap.wd The width of the haplotype tracks to be plotted, generally recommended to be about 0.4 (default value)
#' @param recombination_data List object as returned by the function \code{count_recombinations}. By default \code{NULL}, 
#' in which case no overlay of predicted recombination events is performed. However, it can be useful
#' to visualise predicted recombination events, particularly as this might help inform the choice of argument \code{plausible_pairing_prob}
#' of that function. See \code{\link{count_recombinations}} for more details.
#' @param reset_par By default \code{TRUE}, reset par on exit.
#' @param log Character string specifying the log filename to which standard output should be written. If \code{NULL} log is send to stdout.
#' @return If \code{recombinant_scan} vector is supplied, a vector of recombinant offspring ID in the region of interest (otherwise \code{NULL}).
#' @examples
#' data("IBD_4x")
#' visualiseHaplo(IBD_list = IBD_4x,
#'                display_by = "name",
#'                linkage_group = 1,
#'                select_offspring = "all",
#'                multiplot = c(3,3))
#' @export
visualiseHaplo <- function(IBD_list,
                           display_by = c("phenotype","name"),
                           linkage_group = NULL,
                           Phenotype.df = NULL,
                           genotype.ID = NULL,
                           trait.ID = NULL,
                           pheno_range = NULL,
                           cM_range = "all",
                           highlight_region = NULL,
                           select_offspring = NULL,
                           recombinant_scan = NULL,
                           allele_fish = NULL,
                           presence_threshold = 0.95,
                           xlabl=TRUE,
                           ylabl=TRUE,
                           mainTitle = NULL,
                           multiplot = NULL,
                           append = FALSE,
                           colPal = c("white","navyblue","darkred"),
                           hap.wd = 0.4,
                           recombination_data = NULL,
                           reset_par = TRUE,
                           log = NULL) {
  
  oldpar <- par(no.readonly = TRUE)
  if(reset_par) on.exit(par(oldpar)) 
  
  # Check input
  IBD_list <- test_IBD_list(IBD_list)
  if(hap.wd >= 1) warning("Selecting hap.wd >= 1 will result in overlapping tracks! Suggest to reduce value to 0.4.")
  
  display_by = match.arg(display_by)
  
  ploidy <- IBD_list[[1]]$ploidy
  ploidy2 <- IBD_list[[1]]$ploidy2
  
  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }
  
  if(length(linkage_group)!=1 | !is.numeric(linkage_group)) {
    if(length(IBD_list) != 1) stop("linkage_group should be a single numeric identifier.")
    linkage_group <- 1
  }
  
  IBDarray <- IBD_list[[linkage_group]]$IBDarray
  IBDtype <- IBD_list[[linkage_group]]$IBDtype
  bivalent_decoding <- IBD_list[[linkage_group]]$biv_dec
  gap <- IBD_list[[linkage_group]]$gap
  
  if(IBDtype == "genotypeIBD"){
    
    # if(IBD_list[[1]]$method == "hmm_TO"){
    #   ## We are using the output of TetraOrigin
    #   # Nstates <- Nstates.fun(biv_dec = bivalent_decoding, pl = ploidy, pl2 = ploidy2)
    #   # mname <- paste0("GenotypeMat",Nstates)
    #   # indicatorMatrix <- as.matrix(get(mname, envir = getNamespace("polyqtlR")))
    #   
    #   GenCodes <- t(sapply(IBD_list[[1]]$genocodes,function(x) as.numeric(strsplit(as.character(x),"")[[1]])))
    #   Nstates <- nrow(GenCodes)
    #   indicatorMatrix <- matrix(0,nrow=Nstates, ncol=ploidy + ploidy2)
    #   for(r in 1:nrow(indicatorMatrix)){
    #     for(h in 1:ncol(GenCodes)){
    #       indicatorMatrix[r,GenCodes[r,h]] <- indicatorMatrix[r,GenCodes[r,h]] + 1
    #     }
    #   }
    #   
    # } else{
    ## We are using the output of estimate_IBD or PolyOrigin
    GenCodes <- t(sapply(sapply(IBD_list[[1]]$genocodes,function(x) strsplit(x,"")), function(y) sapply(y,function(z) which(letters == z))))
    Nstates <- nrow(GenCodes)
    indicatorMatrix <- matrix(0,nrow=Nstates, ncol=ploidy + ploidy2)
    for(r in 1:nrow(indicatorMatrix)){
      for(h in 1:ncol(GenCodes)){
        indicatorMatrix[r,GenCodes[r,h]] <- indicatorMatrix[r,GenCodes[r,h]] + 1
      }
    }
    
    # }
    
  } else if(IBDtype == "haplotypeIBD"){
    Nstates <- ploidy + ploidy2
    indicatorMatrix <- NULL
  } else{
    stop("Unable to determine IBD type, please check IBD_list elements have a valid $IBDtype tag.")
  }
  
  if(dim(IBDarray)[2] %% Nstates != 0) stop("Incompatible input detected..")
  
  if(!is.null(gap)) {
    cMpositions <- as.numeric(unlist(lapply(strsplit(dimnames(IBDarray)[[1]],"cM"),function(x) x[2])))
  } else{
    cMpositions <- IBD_list[[linkage_group]]$map$position
    
    # mapdata <- IBD_list[[linkage_group]]$map
    # if("position" %in% colnames(mapdata)){
    #   cMpositions <- mapdata[,"position"]
    # } else if("cM" %in% colnames(mapdata)){
    #   cMpositions <- mapdata[,"cM"]
    # } else{
    #   stop("Could not extract map positions using header 'position' or 'cM' from $map of IBD_list. Suggest to check and manually re-name!")
    # }
    
  }
  
  subcM <- 1:length(cMpositions) #for subsetting later..
  
  if(cM_range[1] != "all"){ #Used cM_range[1] to avoid warnings
    if(min(cM_range) < min(cMpositions) | max(cM_range) > max(cMpositions)){
      stop(paste0("cM_range [",paste(cM_range,collapse = ","),
                  "] does not fall within bounds of IBD input: [",
                  paste(range(cMpositions),collapse=","),"]"))
      
    } else{
      
      subcM <- which(cMpositions >= min(cM_range) & cMpositions <= max(cM_range))
      cMpositions <- cMpositions[cMpositions >= min(cM_range) & cMpositions <= max(cM_range)]
      
    }
  }
  
  if(display_by == "phenotype"){
    if(is.null(Phenotype.df) | is.null(genotype.ID) | is.null(trait.ID) | is.null(pheno_range))
      stop('If displaying haplotypes by using "phenotype", then phenotypes must be correctly supplied. Consult help using ?visualiseHaplo')
    
    # Match phenos and genotypes..
    if(any(is.na(Phenotype.df[,genotype.ID]))){
      warning("Missing genotype.ID values detected. Removing and proceeding without...")
      Phenotype.df <- Phenotype.df[!is.na(Phenotype.df[,genotype.ID]),]
    }
    
    phenoGeno <- intersect(
      dimnames(IBDarray)[[3]],
      unique(Phenotype.df[,genotype.ID]))
    
    phenoGeno <- Phenotype.df[Phenotype.df[,genotype.ID] %in% phenoGeno &
                                !is.na(Phenotype.df[,trait.ID]) &
                                Phenotype.df[,trait.ID] >= pheno_range[1] &
                                Phenotype.df[,trait.ID] <= pheno_range[2],genotype.ID]
    
    if(length(phenoGeno) != length(unique(phenoGeno))) stop("Check input - unequal dimensions of phenotypes and genotype data.\nSuggest to use polyqtlR:::BLUE() to generate consensus phenotypes first.")
    
    plotorder <- order(Phenotype.df[match(phenoGeno,Phenotype.df[,genotype.ID]),trait.ID])
    
  } else{
    if(is.null(select_offspring)) stop('If displaying haplotypes by using "name", then select_offspring must be correctly supplied.')
    
    suppressWarnings(
      if(select_offspring[1] == "all") {
        phenoGeno <- dimnames(IBDarray)[[3]]
      } else {
        phenoGeno <- intersect(dimnames(IBDarray)[[3]],select_offspring)
        if(length(phenoGeno) != length(select_offspring)) warning("Merely ",paste(length(phenoGeno),"of the",length(select_offspring),"offspring names could be matched!"))
      }
    )
    
    plotorder <- 1:length(phenoGeno)
  }
  
  popSize <- length(phenoGeno)
  
  genoProbs <- IBDarray[dimnames(IBDarray)[[1]][subcM],,phenoGeno,drop=FALSE]
  
  if(IBDtype == "genotypeIBD"){
    haploProbs <- array(NA, dim=c(dim(genoProbs)[1], dim(indicatorMatrix)[2], dim(genoProbs)[3]))
    haploProbs[] <- apply(genoProbs,3,function(x) x %*% indicatorMatrix)
    
    dimnames(haploProbs)[c(1,3)] <- dimnames(genoProbs)[c(1,3)]
    dimnames(haploProbs)[[2]] <- paste0("h", 1:(ploidy + ploidy2))
    
  } else{
    haploProbs <- genoProbs
  }
  
  if(!is.null(multiplot) & length(multiplot)==2) par(mfrow = multiplot)
  
  colfunc <- colorRampPalette(colPal)
  colours <- colfunc(201)
  interv <- seq(0,2,0.01)
  
  if(xlabl & ylabl) par(mar=c(4,4,2,0.5))
  
  if(!is.null(mainTitle)){
    if(length(mainTitle) != length(plotorder)){
      warning("mainTitle not same length as number of plots.")
      mainTitle <- suppressWarnings(cbind(mainTitle,plotorder))[,1] #use R's recycling
    }
  } else{
    mainTitle <- as.character(phenoGeno)
  }
  
  for(i in plotorder){
    
    plot(NULL,xlim=range(cMpositions),
         ylim=c(0.75,ploidy+ploidy2),
         xlab=ifelse(xlabl,"cM",""),ylab=ifelse(ylabl,"homologue",""),
         main=mainTitle[i],
         axes=FALSE)
    
    axis(1)
    axis(2,las=1,at = 1:(ploidy + ploidy2))
    
    abline(h = ploidy + 0.5, lwd=2,col="gray50")
    
    ## This part could be improved in a future version to give a background effect:
    # ## Add a layer of background probabilities that appear more smoothed
    # for(j in 1:dim(haploProbs)[1]) {
    #   col.temp <- colours[findInterval(haploProbs[j,,i],interv)]
    #   col.rgb <- col2rgb(col.temp)
    #   
    #   points(rep(cMpositions[j],ploidy+ploidy2),1:(ploidy+ploidy2),
    #          col=rgb(red = col.rgb["red",]/265,
    #                  green = col.rgb["green",]/265,
    #                  blue = col.rgb["blue",]/265,
    #                  alpha = 0.1),
    #          pch=15,cex = 2)
    # }
    
    ## Add the specific IBD haplotype probabilities at the estimated positions:
    
    for(j in 1:dim(haploProbs)[1]) {
      segments(x0 = rep(cMpositions[j],ploidy+ploidy2),y0 = 1:(ploidy+ploidy2) - hap.wd/2,
               x1 = rep(cMpositions[j],ploidy+ploidy2),y1 = 1:(ploidy+ploidy2) + hap.wd/2,
               col=colours[findInterval(haploProbs[j,,i],interv)],lwd=2)
    }
    
    
    
    if(!is.null(highlight_region)){
      if(length(highlight_region) > 2) stop("At most 2 positions can be highlighted on the plot")
      if(!is.numeric(highlight_region)) stop("highlight_region should contain numeric values")
      if(!all(findInterval(x = highlight_region,vec = range(cMpositions)) == 1)) 
        warning("Not all supplied values of highlight_region are contained in the plot range")
      if(length(highlight_region) == 1) 
        abline(v = highlight_region, col = "red")
      else
        rect(highlight_region[1],0.8,highlight_region[2],ploidy+ploidy2+0.2, density = 0,border = "red")
    }
    
    if(!is.null(recombination_data)){
      rec_dat <- recombination_data[[linkage_group]]$recombinations[[phenoGeno[i]]]
      
      if(list.depth(rec_dat) > 1){ #this might need tweaking
        for(v in 1:length(rec_dat)){ #handles multiple (probabilistic) valencies
          bp <- rec_dat[[v]]$recombination_breakpoints
          
          sapply(1:length(bp),function(k){
            if(!is.na(bp[k])){
              rh <- match(strsplit(names(bp[k]),"")[[1]],LETTERS)
              for(m in 1:length(bp[[k]])){
                segments(x0 = bp[[k]][m],y0 = rh[1],x1 = bp[[k]][m],y1 = rh[2],
                         col = rgb(1,0,0,rec_dat[[v]]$Valent_probability),
                         lwd = 2)
                points(rep(bp[[k]][m],length(rh)),rh,pch=19,
                       col = rgb(1,0,0,rec_dat[[v]]$Valent_probability))
              }
              
            }
          })
        }
      }
    }
    
    
  }
  
  # Fill in any remaining empty plot spaces:
  if(!is.null(multiplot) & !append){
    over <- 0
    
    if(prod(multiplot) > popSize) over <- prod(multiplot) - popSize
    if(prod(multiplot) < popSize) {
      if(popSize %% prod(multiplot) != 0)
        over <- prod(multiplot) - (popSize%%prod(multiplot))
    }
    
    if(over > 0) for(i in 1:over) plot(NULL,xlim=c(0,1),ylim=c(0,1),
                                       xlab=NA,ylab=NA,axes=F)
    par(mfrow=c(1,1))
  }
  
  if (!is.null(log))
    close(log.conn)
  
  recombinants <- NULL #Assign NULL to initialise
  
  if(!is.null(recombinant_scan)){
    if(!is.numeric(recombinant_scan) | !(all(recombinant_scan %in% 1:ploidy) | all(recombinant_scan %in% (ploidy+1):(ploidy+ploidy2))))
      stop("Inappropriate input detected for recombinant_scan. Vector must be numeric, and recombinations can only occur between homologues within 1 parent.")
    
    if(length(recombinant_scan) == ploidy & all(recombinant_scan %in% 1:ploidy)){
      #want more than ploidy/2 haplotypes to have a probability > presence_threshold
      recombinants <- names(apply(haploProbs,3,function(x = haploProbs[,,1]){
        length(which(apply(x[,paste0("h",recombinant_scan)] > presence_threshold,2,any))) > ploidy/2
      }))
      
    } else if(length(recombinant_scan) == ploidy2 & all(recombinant_scan %in% (ploidy+1):(ploidy+ploidy2))){
      
      recombinants <- names(apply(haploProbs,3,function(x = haploProbs[,,1]){
        length(which(apply(x[,paste0("h",recombinant_scan)] > presence_threshold,2,any))) > ploidy2/2
      }))
      
    } else{
      #need to identify recombinants between a subset of parental homologues, using an approximate approach
      recombinants <- names(which(apply(haploProbs,3,function(x){
        all(length(which(apply(x[,paste0("h",recombinant_scan)] > presence_threshold,2,any))) == length(recombinant_scan) &
              length(which(findInterval(rowSums(x[,paste0("h",recombinant_scan)]), c(presence_threshold, 3 - 2*presence_threshold)) == 1)) >=
              presence_threshold*nrow(x)) #recycle presence_threshold to also mean the propn of loci carrying approx. 1 copy of these 2 homologues.
      })))
      
    }
  }
  
  alleles_fished <- NULL
  
  if(!is.null(allele_fish)){
    
    if(!all(allele_fish %in% 1:(ploidy+ploidy2))) stop(paste("Only alleles", paste(1:(ploidy+ploidy2),collapse = ","), "are allowed in allele_fish"))
    if(length(allele_fish) > (ploidy+ploidy2)/2) stop(paste("Only", (ploidy+ploidy2)/2, "alleles can possibly be inherited in one individual!"))
    
    carrying_alleles.comb <- dimnames(haploProbs)[[3]][which(apply(haploProbs,3,function(x) all(colMeans(x[,paste0("h",allele_fish),drop = FALSE]) >= presence_threshold)))]
    
    if(length(allele_fish) > 1){
      
      carrying_alleles.sing <- setNames(lapply(allele_fish, function(a) {
        dimnames(haploProbs)[[3]][which(apply(haploProbs,3,function(x) all(colMeans(x[,paste0("h",a),drop=FALSE]) >= presence_threshold)))]
      }),paste0("h",allele_fish))
      
      alleles_fished <- c(setNames(list(carrying_alleles.comb),paste0("h",paste0(allele_fish,collapse = "_"))),
                          carrying_alleles.sing)
      
    } else{
      alleles_fished <- setNames(list(carrying_alleles.comb),paste0("h",paste0(allele_fish,collapse = "_")))
    }
  }
  
  invisible(list(recombinants = recombinants,
                 allele_fish = alleles_fished))
  
} #visualiseHaplo

#' Visualise pairing of parental homologues
#' @description Function to visualise the pairing of parental homologues across the population using graph, with nodes 
#' to denote parental homologues and edges to denote deviations from expected proportions under a polysomic model of inheritance
#' @param meiosis_report.ls List output of function \code{\link{meiosis_report}}
#' @param pos.col Colour corresponding to excess of pairing associations predicted (positive deviations), by default red 
#' @param neg.col Colour corresponding to lack of pairing associations predicted (negative deviations), by default blue
#' @param parent The parent, either "P1" (mother) or "P2 (father)
#' @param max.lwd Maximum line width, by default 20
#' @param datawidemax This argument is currently a work-around to allow multiple plots to have the same scale (line thicknesses consistent).
#' No default is provided. To estimate this value, simply set argument \code{return.data = TRUE}, and record the
#' maximum absolute value over columns 'count', which are the deviations from random expectations. This should be done
#' over multiple function calls if e.g. comparing both P1 and P2 values. When a global maximum (absolute) deviation is known,
#' re-run the function with this value for \code{datawidemax}. The line width specified by \code{max.lwd} will then be
#'  used for this, and all other line widths re-scaled accordingly.
#' @param add.label Should a label be applied, giving the maximum deviation in the plot? By default \code{TRUE}
#' @param return.data Should plot data be returned? By default \code{FALSE}
#' @param \dots Optional arguments passed to \code{\link[igraph]{plot.igraph}}
#' @return If \code{return.data = TRUE}, the values for pairwise deviations from the expected numbers are 
#' returned, useful for determining the value \code{datawidemax} to provide consistent scaling across multiple plots
#' @examples
#' data("mr.ls")
#' visualisePairing(meiosis_report.ls = mr.ls,
#'                  parent = "P1",
#'                  datawidemax = 3)
#' @export
visualisePairing <- function(meiosis_report.ls, 
                             pos.col = "red", 
                             neg.col = "blue",
                             parent,
                             max.lwd = 20,
                             datawidemax,
                             add.label = TRUE,
                             return.data = FALSE,
                             ...){
  
  parent <- match.arg(parent, choices = c("P1","P2"))
  
  out.ls <- setNames(lapply(seq(length(meiosis_report.ls)), function(l){
    ploidy <- ifelse(parent == "P1",
                     meiosis_report.ls[[l]]$ploidy,
                     meiosis_report.ls[[l]]$ploidy2)
    
    pairmat <- meiosis_report.ls[[l]][[paste0(parent,"_pairing")]]
    
    countmat <- pairmat[,1:ploidy]
    
    if(ncol(countmat) != nrow(countmat)) stop("pairmat should be checked.")
    cs <- t(combn(rownames(countmat),2))
    edges <- as.data.frame(cbind(cs,countmat[lower.tri(countmat)]))
    colnames(edges) <- c("h1","h2","count")
    
    ## Update 09/23: Add the chm name to the homologue name for vertex labels
    edges$h1 <- paste0(as.character(l),edges$h1)
    edges$h2 <- paste0(as.character(l),edges$h2)
    
    edges$count <- as.numeric(edges$count)
    ## express the counts as deviations from the mean:
    m.cnt <- mean(countmat[1,2:ploidy])
    edges$count <- edges$count - m.cnt
    edges$sign <- edges$count > 0
    edges$count <- abs(edges$count)
    edges$col <- rep(neg.col,nrow(edges))
    edges$col[edges$sign] <- pos.col
    
    ## Need to re-scale so that weights are equal across plots.
    edges$rescaled <- max.lwd*edges$count/datawidemax
    
    g <- igraph::graph_from_data_frame(edges[,c("h1","h2")], directed = FALSE)
    l <- igraph::layout_in_circle(g)
    
    xtrm <- edges[which.max(edges$count),]
    x <- unlist(round(xtrm[3]))
    if(!xtrm$sign) x <- -x
    
    homs <- paste0("(",paste(xtrm[1:2], collapse = " & "),")")
    
    igraph::plot.igraph(g,layout = l, 
                        vertex.label.color = "black",
                        vertex.label.family = "sans",
                        edge.width = edges$rescaled,
                        edge.color = edges$col,
                        vertex.color = "white",
                        ...)
    
    if(add.label) text(x = 0, y = 1.2*par("usr")[4], 
                       labels = bquote(delta[max] ~ .(homs) == .(x)), cex = 1.2)
    
    if(return.data) {
      return(edges)
    } else{
      return(NULL)
    }
    
  }),
  names(meiosis_report.ls))
  
  return(invisible(out.ls))
  
} #visualisePairing

#' Visualise QTL homologue effects around a QTL position
#' @description Function to visualise the effect of parental homologues around a QTL peak across the population.
#' @param IBD_list List of IBD probabilities
#' @param Phenotype.df A data.frame containing phenotypic values
#' @param genotype.ID The colname of \code{Phenotype.df} that contains the population identifiers (F1 names) (must be a colname of \code{Phenotype.df})
#' @param trait.ID The colname of \code{Phenotype.df} that contains the response variable to use in the model (must be a colname of \code{Phenotype.df})
#' @param linkage_group Numeric identifier of the linkage group being tested, based on the order of \code{IBD_list}.
#' Only a single linkage group is allowed.
#' @param LOD_data  Output of \code{\link{QTLscan}} function
#' @param cM_range If required, the plotting region can be restricted to a specified range of centiMorgan positions
#' (provided as a vector of start and end positions).
#' @param col.pal Vector of colours to use in the visualisations (it is best to provide two or three colours for simplicity). By default,
#' effects will be coloured from purple to green through white.
#' @param point.density Parameter to increase the smoothing of homologue effect tracks
#' @param zero.sum How allele substitution effect should be defined. If \code{FALSE} (by default), the effect of each homologue
#' is computed relative to the overall phenotypic mean, otherwise contrasts (against offspring without the inherited homologue) are used.
#' @param allelic_interaction By default \code{FALSE}, in which case the additive effects of parental alleles are visualised. If \code{TRUE}, a plot of
#' the mean effect of combinations of parental alleles is visualised instead. \code{exploreQTL_output} is required in this case.
#' @param exploreQTL_output If \code{allelic_interaction = TRUE}, the output of the function \code{\link{exploreQTL}} must be provided.
#' @param return_plotData Logical, by default \code{FALSE}. If \code{TRUE}, plot data is returned, otherwise \code{NULL}.
#' @return The estimated effects of the homologues, used in the visualisation
#' @examples
#' data("IBD_4x","BLUEs.pheno","qtl_LODs.4x")
#' visualiseQTLeffects(IBD_list = IBD_4x,
#'                     Phenotype.df = BLUEs.pheno,
#'                     genotype.ID = "Geno",
#'                     trait.ID = "BLUE",
#'                     linkage_group = 2,
#'                     LOD_data = qtl_LODs.4x)
#' @export
visualiseQTLeffects <- function(IBD_list,
                                Phenotype.df,
                                genotype.ID,
                                trait.ID,
                                linkage_group,
                                LOD_data,
                                cM_range = NULL,
                                col.pal = c("purple4","white","seagreen"),
                                point.density = 50,
                                zero.sum = FALSE,
                                allelic_interaction = FALSE,
                                exploreQTL_output = NULL,
                                return_plotData = FALSE){
  
  oldpar <- par(no.readonly = TRUE)    
  on.exit(par(oldpar))            # Thanks Gregor
  
  ploidy <- IBD_list[[1]]$ploidy
  ploidy2 <- IBD_list[[1]]$ploidy2
  
  ## Test input:
  if(length(linkage_group)!=1 | !is.numeric(linkage_group))
    stop("linkage_group should be a single numeric identifier.")
  
  IBD_list <- test_IBD_list(IBD_list)
  
  IBDarray <- IBD_list[[linkage_group]]$IBDarray
  IBDtype <- IBD_list[[linkage_group]]$IBDtype
  bivalent_decoding <- IBD_list[[linkage_group]]$biv_dec
  gap <- IBD_list[[linkage_group]]$gap
  
  ## Update v.0.1.0, add visualisation of allele interactions at QTL peak (different plot type)
  if(allelic_interaction){
    ## Initial checks
    if(is.null(exploreQTL_output)) stop("Argument exploreQTL_output must also be provided")
    if(IBDtype!= "genotypeIBD") stop("Cannot visualise allelic interactions using haplotype probabilities! Re-run estimate_IBD with HMM method.")
    if(!"genotype.means" %in% names(exploreQTL_output)) stop("Incompatible exploreQTL_output detected. Make sure this is output from the function exploreQTL (>= v.0.1.0)")
    
    means <- exploreQTL_output$genotype.means
    # means <- means[!is.na(means),1] #convert to a named vector
    
    genos <- lapply(rownames(means), function(x) match(strsplit(x,"")[[1]],letters))
    
    ## Split this out into parent1 and parent2 genos:
    pars <- lapply(genos, function(x) list(P1 = paste0(x[x %in% 1:ploidy],collapse = ""),
                                           P2 = paste0(x[x %in% (ploidy+1):(ploidy + ploidy2)],collapse = "")))
    
    P1labels <- sort(unique(sapply(pars,'[[',"P1")))
    P2labels <- sort(unique(sapply(pars,'[[',"P2")))
    
    ## The dimensions are defined by ploidy, ploidy2 and bivdec:
    P1dim <- Nstates.onepar(biv_dec = bivalent_decoding,pl = ploidy)
    P2dim <- Nstates.onepar(biv_dec = bivalent_decoding,pl = ploidy2)
    
    if(nrow(means) != P1dim*P2dim) stop(paste0("Expected nrow of exploreQTL_output$genotype.means was ",P1dim*P2dim,". Detected ",nrow(means)," rows!"))
    
    ## Plot effects rather than means of the phenotypes in each class
    means <- means - mean(means, na.rm = TRUE) #i.e. centred at 0
    
    ## Generate colour matrix for allele effects:
    colours <- colorRampPalette(col.pal)(100)
    
    colmin <- range(means,na.rm = TRUE)[1] #These need to be defined for later use by color.bar
    colmax <- range(means,na.rm = TRUE)[2]
    
    colbreaks <- seq(colmin, colmax, (colmax - colmin)/100)
    
    cols <- apply(means, 1, function(coln) {
      if(is.na(coln)){
        return("white")
      } else{
        return(colours[findInterval(x = coln,vec = colbreaks)])
      }
    })
    
    cols[is.na(cols)] <- colours[100] #the max is undefined.
    
    plot_title <- paste("Genotype class effects at",
                        exploreQTL_output$cM,"cM on LG",exploreQTL_output$linkage_group)
    
    layout(matrix(1:2,ncol=2),
           widths = c(1,0.11))
    plot(NULL, xlim = c(0.75,P2dim+0.25),ylim = c(0.75,P1dim+0.25),
         xlab = "P2 alleles",ylab = "P1 alleles",axes = FALSE,
         main = plot_title)
    axis(1,labels = P2labels,at = 1:P2dim)
    axis(2,labels = P1labels,at = 1:P1dim)
    box()
    
    for(i in 1:length(pars)){
      xpos <- match(pars[[i]]$P2,P2labels)
      ypos <- match(pars[[i]]$P1,P1labels)
      
      rect(xleft = xpos - 0.25,
           ybottom = ypos - 0.25,
           xright = xpos + 0.25,
           ytop = ypos + 0.25,
           density = NA,
           col = cols[i],lwd = 1)
    }
    
    ## Now add the colour bar on the right:
    par(mar=c(0,2,0,0))
    colour.bar(col.data = colorRampPalette(col.pal)(10),
               min=colmin,max=colmax,nticks=8,
               ticks=round(seq(colmin, colmax, len=8),2),
               cex.ticks=0.8)
    
    
  } else{
    
    if(IBDtype == "genotypeIBD"){
      ## We are using the output of estimate_IBD
      GenotypeCodes <- t(sapply(sapply(IBD_list[[1]]$genocodes,function(x) strsplit(x,"")), function(y) sapply(y,function(z) which(letters == z))))
      Nstates <- nrow(GenotypeCodes)
      indicatorMatrix <- matrix(0,nrow=Nstates, ncol=ploidy + ploidy2)
      for(r in 1:nrow(indicatorMatrix)){
        for(h in 1:ncol(GenotypeCodes)){
          indicatorMatrix[r,GenotypeCodes[r,h]] <- indicatorMatrix[r,GenotypeCodes[r,h]] + 1
        }
      }
    } else if(IBDtype == "haplotypeIBD"){
      Nstates <- ploidy + ploidy2 #this is not actually needed..
      indicatorMatrix <- NULL
    } else{
      stop("Unable to determine IBD type, please check IBD_list elements have a valid $IBDtype tag.")
    }
    
    if(dim(IBDarray)[2] %% Nstates != 0) stop("Incompatible input detected..")
    
    if(!is.null(gap)) {
      cMpositions <- as.numeric(unlist(lapply(strsplit(dimnames(IBDarray)[[1]],"cM"),function(x) x[2])))
    } else{
      cMpositions <- IBD_list[[linkage_group]]$map$position
    }
    
    subcM <- 1:length(cMpositions) #for subsetting later..
    
    if(is.null(cM_range)){
      cM_range <- cMpositions
      index_range <- dimnames(IBDarray)[[1]]
    } else if(min(cM_range) < min(cMpositions) | max(cM_range) > max(cMpositions)) {
      message("Incompatible cM_range detected. Proceeding with the following positions:")
      print(cMpositions)
      cM_range <- cMpositions
      index_range <- dimnames(IBDarray)[[1]]
    } else{ #cM_range is properly specified. Fish out identifiers..
      cM_range <- cMpositions[cMpositions >= min(cM_range) & cMpositions <= max(cM_range)]
      index_range <- dimnames(IBDarray)[[1]][cMpositions >= min(cM_range) & cMpositions <= max(cM_range)]
    }
    
    ## First get the starting set - the ones that were phenotyped and genotyped:
    if(any(is.na(Phenotype.df[,genotype.ID]))){
      warning("Missing genotype.ID values detected. Removing and proceeding without...")
      Phenotype.df <- Phenotype.df[!is.na(Phenotype.df[,genotype.ID]),]
    }
    
    phenoGeno0 <- intersect(
      dimnames(IBDarray)[[3]],
      unique(Phenotype.df[,genotype.ID]))
    
    pheno <- Phenotype.df[Phenotype.df[,genotype.ID] %in% phenoGeno0 &
                            !is.na(Phenotype.df[,trait.ID]),c(genotype.ID,trait.ID)]
    pheno <- pheno[order(pheno[,1]),]
    
    phenoGeno <- as.character(unique(pheno[,1]))
    popSize <- length(phenoGeno)
    Phenotypes <- pheno[,2]
    
    if(popSize != nrow(pheno)) stop("Check input - unequal dimensions of phenotypes and genotype data.\nSuggest to use polyqtlR:::BLUE() to generate consensus phenotypes first.")
    
    prob_matched <- match(phenoGeno,dimnames(IBDarray)[[3]])
    
    calc.AlEffects <- function(cM.pos=index_range[1], tag){
      
      if(tag == "genotypeIBD"){
        haploProbs <- t(IBDarray[cM.pos,,prob_matched]) %*% indicatorMatrix
      } else{ #tag == "haplotypeIBD"
        haploProbs <- t(IBDarray[cM.pos,,prob_matched])
      }
      
      mean_sd_weights <- function(wghts,min.wghts, zerosum = zero.sum){
        
        plus.mean <- sum(wghts*Phenotypes)/sum(wghts)
        minus.mean <- sum(min.wghts*Phenotypes)/sum(min.wghts)
        if(zerosum){
          return(plus.mean - minus.mean) #here, only interested in the difference between presence and absense
        } else{
          return(plus.mean - mean(Phenotypes)) #just return the mean of presence
        }
      } #mean_sd_weights
      
      output <- setNames(
        sapply(1:(ploidy + ploidy2), function(sub.config){
          weights <- haploProbs[,sub.config]
          min.weights <- 1-haploProbs[,sub.config]
          mean_sd_weights(weights,min.weights)
        }),
        paste0("H",1:(ploidy + ploidy2))
      )
    } #calc.AlEffects
    
    AlEffects <- do.call(rbind,lapply(index_range,calc.AlEffects, IBDtype))
    rownames(AlEffects) <- cM_range
    ## Generate colour matrix for allele effects:
    colours <- colorRampPalette(col.pal)(100)
    
    colmin <- range(AlEffects)[1] #These need to be defined for later use by color.bar
    colmax <- range(AlEffects)[2]
    
    colbreaks <- seq(colmin, colmax, (colmax - colmin)/100)
    
    ## Now visualise:
    layout(matrix(c(1:(ploidy+ploidy2+1),
                    rep(ploidy+ploidy2+2,ploidy+ploidy2+1)),
                  ncol=2),
           heights = c(1, rep(0.1,ploidy + ploidy2 - 1), 0.3),
           widths = c(1,0.11)) #afterwards add a side panel for the colour bar...
    par(mar=c(1,4,1,0.5)) #adjust this later
    
    plotData <- LOD_data$QTL.res[LOD_data$QTL.res$chromosome == linkage_group &
                                   LOD_data$QTL.res$position %in% cM_range,]
    with(plotData,
         plot(position,LOD,type="l",lwd=2,
              xlab="",ylab="LOD",axes=FALSE,cex.lab=1.2))
    axis(2,las=1); box()
    
    ## Add a threhold if available:
    if(!is.null(LOD_data$Perm.res$threshold))
      abline(h = LOD_data$Perm.res$threshold, lwd = 2, lty = 3, col = "darkred")
    
    ## I want to spline the colours to get a nicer spread on the plot..
    span.param <- ceiling(point.density*(par("usr")[2]-par("usr")[1]))
    spreadX <- seq(cM_range[1],cM_range[length(cM_range)],
                   (cM_range[length(cM_range)] - cM_range[1])/span.param)
    
    PlotAlEffects <- apply(AlEffects, 2, function(coln) predict(smooth.spline(y = coln,x = cM_range),spreadX)$y)
    
    ## Any over-fitting outside the range of AlEffects needs to be suppressed (only interpolation allowed!)
    PlotAlEffects[PlotAlEffects>max(AlEffects)] <- max(AlEffects)
    PlotAlEffects[PlotAlEffects<min(AlEffects)] <- min(AlEffects)
    
    cols <- apply(PlotAlEffects, 2, function(coln) colours[findInterval(x = coln,vec = colbreaks)])
    cols[is.na(cols)] <- colours[100] #the max is undefined.
    
    par(mar=c(0,4,0,0.5))
    for(h in (ploidy + ploidy2):2) {
      plot(spreadX,rep(0.5,nrow(PlotAlEffects)),
           ylim=c(0,1),col=cols[,h],pch=15,xlab="",ylab="", axes=FALSE)
      mtext(side=2,text = paste("H",h),line = 0,las=1)}
    # Add a bottom axis:
    par(mar=c(4,4,0,0.5))
    plot(spreadX,rep(0.5,nrow(PlotAlEffects)),
         ylim=c(0,1),col=cols[,1],pch=15,xlab="cM",ylab="", cex.lab=1.2,
         axes=FALSE)
    mtext(side=2,text = "H 1",line = 0,las=1)
    axis(1)
    ## Now add the colour bar on the right:
    par(mar=c(0,2,0,0))
    colour.bar(col.data = colorRampPalette(col.pal)(10),
               min=colmin,max=colmax,nticks=8,
               ticks=round(seq(colmin, colmax, len=8),2),
               cex.ticks=0.8)
    
    if(return_plotData){
      return(AlEffects)
    } else {
      return(NULL)
    }
  }
} #visualiseQTLeffects

## This needs to be manually added to the NAMESPACE after each call to devtools::document:
# importFrom(Rcpp, evalCpp)
# useDynLib(polyqtlR, .registration = TRUE)
