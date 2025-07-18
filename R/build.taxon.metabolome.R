#' build.taxon.metabolome
#'
#' utilizes downloaded and properly formatted local pubchem data created by 'get.pubchem.ftp' function to filter a dataset created by 'build.pubchem.bio' function
#' @details utilizes downloaded and properly formatted local pubchem data created by 'get.pubchem.ftp' function
#' @param pc.directory directory from which to load pubchem .Rdata files
#' @param taxid integer vector of integer NCBI taxonomy IDs.  i.e.  c(9606, 1425170 ) for Homo sapiens and Homo heidelbergensis.    
#' @param get.properties logical. if TRUE, will return rcdk calculated properties:  XLogP, TPSA, HBondDonorCount and HBondAcceptorCount.
#' @param full.scored logincal.  default = FALSE.  When false, only metabolites which map to the taxid(s) are returned.  When TRUE, all metabolites are returned, with scores assigned based on the distance of non-mapped metabolites to the root node.  i.e. specialized metabolites from distantly related species are going to be scored at or near zero, specialized metabolites of mores similar species higher, and more conserved metabolites will score higher than ore specialized. 
#' @param aggregation.function function. default = max.  can use mean, median, min, etc, or a custom function.  Defines how the aggregate score will be calculated when multiple taxids are used.
#' @param threads integer. how many threads to use when calculating rcdk properties.  parallel processing via DoParallel and foreach packages.  
#' @param db.name character. what do you wish the file name for the saved version of this database to be?  default = 'custom.metabolome', but could be 'taxid.4071' or 'Streptomyces', etc.  Saved as an .Rdata file in the 'pc.directory' location. 
#' @param rcdk.desc vector. character vector of valid rcdk descriptors.  default = rcdk.desc <- c("org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor", "org.openscience.cdk.qsar.descriptors.molecular.AcidicGroupCountDescriptor", "org.openscience.cdk.qsar.descriptors.molecular.BasicGroupCountDescriptor", "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor"). To see descriptor categories: 'dc <- rcdk::get.desc.categories(); dc' .  To see the descriptors within one category: 'dn <- rcdk::get.desc.names(dc\[4\]); dn'. Note that the four default parameters are relatively fast to calculate - some descriptors take a very long time to calculate.  you can calculate as many as you wish, but processing time will increase the more descriptors are added.   
#' @return a data frame containing pubchem CID ('cid'), and lowest common ancestor ('lca') NCBI taxonomy ID integer. will also save to pc.directory as .Rdata file.
#' @author Corey Broeckling
#' 
#' @export
#' 
build.taxon.metabolome <- function(
    pc.directory = NULL,
    taxid = c(),
    get.properties = TRUE,
    full.scored = FALSE,
    aggregation.function = max,
    threads = 8,
    db.name = "custom.metabolome", 
    rcdk.desc = c(
      "org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor",
      "org.openscience.cdk.qsar.descriptors.molecular.AcidicGroupCountDescriptor",
      "org.openscience.cdk.qsar.descriptors.molecular.BasicGroupCountDescriptor",
      "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor"
    )
) {
  
  if(length(taxid) == 0) {
    stop("please list at least one integer taxid, i.e. 'taxid = c(4071, 4081)'", '\n')
  }
  
  if(file.exists(paste0(pc.directory, "/cid.lca.Rdata"))) {
    load(paste0(pc.directory, "/cid.lca.Rdata"))
  } else {
    stop(paste0(pc.directory, "/cid.lca.Rdata"), "does not exist", '\n')
  }
  
  if(file.exists(paste0(pc.directory, "/taxid.hierarchy.Rdata"))) {
    load(paste0(pc.directory, "/taxid.hierarchy.Rdata"))
  } else {
    stop(paste0(pc.directory, "/taxid.hierarchy.Rdata"), "does not exist", '\n')
  }
  
  if(file.exists(paste0(pc.directory, "/pc.bio.Rdata"))) {
    load(paste0(pc.directory, "/pc.bio.Rdata"))
  } else {
    stop(paste0(pc.directory, "/pc.bio.Rdata"), "does not exist", '\n')
  }
  
  taxid.hierarchy <- taxid.hierarchy
  cid.lca <- cid.lca
  pc.bio <- pc.bio
  
  `%dopar%` <- foreach::`%dopar%`
  
  out <- pc.bio
  
  for(i in 1:length(taxid)) {
    tax.match <- which(taxid.hierarchy == taxid[i], arr.ind = TRUE)
    taxid.v <- as.vector(t(data.frame(taxid.hierarchy[tax.match[1,1],])))
    metabolome <- unique(cid.lca$cid[cid.lca$lca %in% taxid.v])
    keep <- which(metabolome %in% pc.bio$cid)
    metabolome <- metabolome[keep]
    ## metabolome <- metabolome[keep]
    ## which(metabolome == 1548943)
    
    if(full.scored) {
      ## get taxid hierarchy. store as vector.  
      ## compare to each out taxid vector by lca
      cat(taxid[i], " ", length(keep), 'mapped metabolites')
      taxid.row <- tax.match[1,1]
      taxid.column <- as.integer(tax.match[1,2])
      taxid.vector <- as.vector(unlist(taxid.hierarchy[taxid.row, 1:ncol(taxid.hierarchy)]))
      
      ## tmp is the index to the correct cid.lca row for each pubchem row
      th.ind <- which(names(cid.lca) == "species"):ncol(cid.lca)
      # th.ind <- th.ind[taxid.column:length(th.ind)]
      tmp <- match(pc.bio$cid, cid.lca$cid)
      cat(" -- calculating similarities", '\n')
      
      do.sim <- which(!is.na(tmp))
      
      error <- NA
      j <- 1
      
      doParallel::registerDoParallel(cl <- parallel::makeCluster(threads))
      results <- foreach::foreach(j = do.sim) %dopar% {
        tryCatch(
          #this is the chunk of code we want to run
          {
            mtch.col <- suppressWarnings((which(taxid.vector == cid.lca[tmp[j], th.ind])[1]) - taxid.column + 1)
            mtch.col
            #when it throws an error, the following block catches the error
          }, error = function(msg){
            stop("error on", j, '\n')
            return(NA)
          }
        )
      }
      
      results <- unlist(results)
      tax.lca.sim <- rep(NA, length(tmp))
      tax.lca.sim[do.sim] <- results
      tax.lca.sim <- round((max(tax.lca.sim, na.rm = TRUE) - tax.lca.sim)/length(th.ind), 4)
      tax.lca.sim[pc.bio$cid %in% metabolome] <- 1
      
      
      
      # tax.lca.sim <- rep(NA, length(tmp))
      # tmp.ind <- 1:length(tmp)
      # tmp.chunks <- split(tmp, ceiling(seq_along(tmp.ind)/100000))
      # out.chunks <- as.list(rep(NA, length(tmp.chunks)))
      # for(y in 1:length(tmp.chunks)) {
      #   for(x in 1:length(tmp.chunks[[y]])) {
      #     out.vec <- rep(NA, length(tmp.chunks))
      #     if(!is.na(tmp.chunks[[y]][x])) {
      #       out.vec[x] <- tryCatch(
      #         #this is the chunk of code we want to run
      #         {
      #           mtch.col <- which(taxid.vector == cid.lca[tmp.chunks[[y]][x], th.ind])[1]
      #           mtch.col
      #           #when it throws an error, the following block catches the error
      #         }, error = function(msg){
      #           stop("error on", i, '\n')
      #           return(NA)
      #         }
      #       )
      #     }
      #   }
      # }
      
      # do.sim <- which(!is.na(tmp))
      # # do.sim <- do.sim[151600:length(do.sim)]
      # do.sim.sim <- sapply((do.sim), FUN = function(x) {
      # # tax.lca.sim <- sapply((800000:length(tmp)), FUN = function(x) {
      #  cat(x, ' ')
      #     tryCatch(
      #       #this is the chunk of code we want to run
      #       {
      #         mtch.col <- which(taxid.vector == cid.lca[tmp[x], th.ind])[1]
      #         mtch.col
      #         #when it throws an error, the following block catches the error
      #       }, error = function(msg){
      #         stop("error on", i, '\n')
      #         return(NA)
      #       }
      #       )
      # }
      # )
      # tax.lca.sim <- rep(NA, length(tmp))
      # tax.lca.sim[do.sim] <- do.sim.sim
      # tax.lca.sim <- (max(tax.lca.sim, na.rm = TRUE) - tax.lca.sim)/max(tax.lca.sim, na.rm = TRUE)
      # tax.lca.sim[keep] <- 1
      
      out[,paste0("taxonomy.lca.similarity.", taxid[i])] <- tax.lca.sim
      
      doParallel::stopImplicitCluster()
    }
  }
  
  if(!full.scored) {
    cat("keeping", length(keep), "metabolites", '\n')
    out <- pc.bio[keep, ]
    if(nrow(out) == 0) {
      error("no metabolites found for taxid(s):", paste0(taxid, collapse = ", "))
    }
  } else {
    
    ## aggregate taxonomy.lca.similarity scores using assigned function
    use.cols <- which(grepl("taxonomy.lca.similarity.", names(out)))
    suppressWarnings(agg.sim <- apply(out[,use.cols, drop = FALSE], 1, aggregation.function, na.rm = TRUE, simplify = TRUE))
    agg.sim[is.infinite(agg.sim)] <- NA
    # agg.sim <- agg.sim[unique(c(which(is.infinite(agg.sim)), which(is.na(agg.sim))))] <- NA
    out[,paste0("taxonomy.lca.similarity.", "aggregate")] <- agg.sim
  }
  
  if(get.properties) {
    cat(" - calclulating rcdk properties",  format(Sys.time()), '\n')
    cid.list <- as.list(out$cid)
    sm.list <- as.list(out$smiles)
    doParallel::registerDoParallel(cl <- parallel::makeCluster(threads))
    results <- foreach::foreach(i = 1:(length(cid.list))) %dopar% {
      desc <- rcdk.desc
      mol <- rcdk::parse.smiles(sm.list[[i]])
      
      names(mol) <- cid.list[[i]]
      if(is.null(mol)) {
        descs <- rep(NA, length(desc))
      } else {
        descs <- rcdk::eval.desc(mol, desc)
      }
      
      descs
    }
    doParallel::stopImplicitCluster()
    
    
    results.df <- do.call("rbind", results)
    out <- out[order(out$cid),]
    results.df <- results.df[order(as.numeric(row.names(results.df))),]
    
    out <- data.frame(
      out,
      results.df
    )
    doParallel::stopImplicitCluster()
  }
  
  return(out)
  doParallel::stopImplicitCluster()
  save(out, file = paste0(pc.directory, "/", db.name, "Rdata"))
}

# pc.bio.sub <- build.taxon.metabolome(taxid = c(4072, 4107, 4047), pc.directory = "C:/Temp/20250703", get.properties = FALSE, full.scored = TRUE)
# pc.bio.sub[pc.bio.sub$cid %in% c(311, 174174, 1548943, 21585658, 139590519),
#            c("cid", "name", "taxonomy.lca.similarity.4072", "taxonomy.lca.similarity.4107", "taxonomy.lca.similarity.aggregate")]
#           ##  citric acid, atropine, capsaicin, daptomycin, saccharomonopyrone C
# load("C:/Temp/20250703/cid.lca.Rdata")
# cid.lca[cid.lca$cid %in% 1548943,]
# load("C:/Temp/20250703/taxid.hierarchy.Rdata")
# sub.taxid.hierarchy <- taxid.heirarchy[taxid.heirarchy$species %in% c(1173, 4072, 4081, 4232, 4932)]
# save(sub.taxid.hierarchy, file = "//csunts.acns.colostate.edu/arc/cbroeckl/Documents/GitHub/pubchem.bio/inst/extdata/sub.taxid.hierarchy.Rda")
# pc.bio.subset <- pc.bio.sub[pc.bio.sub$cid %in% c(311, 174174, 1548943, 21585658, 139590519),]
# save(pc.bio.subset, file = "//csunts.acns.colostate.edu/arc/cbroeckl/Documents/GitHub/pubchem.bio/inst/extdata/pc.bio.tax.scored.subset.Rda")
