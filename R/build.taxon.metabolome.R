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
#' @param pubchem.bio.object R data.table, generally produced by build.pubchem.bio; preferably, define pc.directory
#' @param cid.lca.object R data.table, generally produced by build.cid.lca; preferably, define pc.directory
#' @param taxid.hierarchy.object R data.table, generally produced by get.pubchem.ftp; preferably, define pc.directory
#' @param output.directory directory to which the pubchem.bio database is saved.  If NULL, will try to save in pc.directory (if provided), else not saved. 
#' @param keep.scored.only logical.  If TRUE, biological metabolites with NA for the taxonomy score are removed before returning.  
#' @return a data frame containing pubchem CID ('cid'), and lowest common ancestor ('lca') NCBI taxonomy ID integer. will also save to pc.directory as .Rdata file.
#' @examples
#' data('cid.lca', package = "pubchem.bio")
#' data('pubchem.bio', package = "pubchem.bio")
#' data('taxid.hierarchy', package = "pubchem.bio")
#' my.taxon.db <- build.taxon.metabolome(
#' pubchem.bio.object = pubchem.bio,
#' cid.lca.object = cid.lca, taxid.hierarchy.object = taxid.hierarchy,
#' get.properties = FALSE, threads = 1, taxid = c(1))
#' head(my.taxon.db)
#' @author Corey Broeckling
#' 
#' @export
#' 
#' 
build.taxon.metabolome <- function(
    pc.directory = NULL,
    taxid = c(),
    get.properties = FALSE,
    full.scored = TRUE,
    keep.scored.only = FALSE,
    aggregation.function = max,
    threads = 8,
    db.name = "custom.metabolome", 
    rcdk.desc = c(
      "org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor",
      "org.openscience.cdk.qsar.descriptors.molecular.AcidicGroupCountDescriptor",
      "org.openscience.cdk.qsar.descriptors.molecular.BasicGroupCountDescriptor",
      "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor"
    ),
    pubchem.bio.object = NULL,
    cid.lca.object = NULL, 
    taxid.hierarchy.object = NULL,
    output.directory = NULL
) {
  
  # loadNamespace("data.table")
  
  out.dir <- pc.directory
  if(is.null(out.dir)) out.dir <- output.directory
  
  if(is.null(pc.directory) & is.null(cid.lca.object)) {
    stop("if you opt to note define the pc.directory, you must provide ALL of 'cid.lca.object', 'taxid.hierarchy.object', 'pubchem.bio.object' variables", '\n')
  }
  
  if(length(taxid) == 0) {
    stop("please list at least one integer taxid, i.e. 'taxid = c(4071, 4081)'", '\n')
  }
  
  if(is.null(cid.lca.object)) {
    load(paste0(pc.directory, "/cid.lca.Rdata"))
    cid.lca <- cid.lca
  } else {
    cid.lca <- cid.lca.object
  }
  
  if(is.null(taxid.hierarchy.object)) {
    load(paste0(pc.directory, "/taxid.hierarchy.Rdata"))
    taxid.hierarchy <- taxid.hierarchy
  } else {
    taxid.hierarchy <- taxid.hierarchy.object
  }
  
  if(is.null(pubchem.bio.object)) {
    load(paste0(pc.directory, "/pc.bio.Rdata"))
    pc.bio <- pc.bio
  } else {
    pc.bio <- pubchem.bio.object
  }
  
  `%dopar%` <- foreach::`%dopar%`
  
  out <- pc.bio
  
  cl <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  base::on.exit({
    parallel::stopCluster(cl)
    rm(cl)
  }) 

  ## store maximum taxid value so we can use it later for eliminating NA values for fast matching.
  first.tax.col <- grep("species", names(cid.lca))[1]
  tmp.df <- as.numeric(as.matrix(cid.lca[,1:ncol(cid.lca)]))
  max.taxid <- max(tmp.df, na.rm = TRUE)
  rm(tmp.df); gc()
  
  for(i in 1:length(taxid)) {
    tax.match <- which(taxid.hierarchy == taxid[i], arr.ind = TRUE)
    taxid.v <- as.vector(t(data.frame(taxid.hierarchy[tax.match[1,1],])))
    metabolome <- unique(cid.lca$cid[cid.lca$lca %in% taxid.v])
    keep <- which(metabolome %in% pc.bio$cid)
    metabolome <- metabolome[keep]
    ## metabolome <- metabolome[keep]
    ## which(metabolome == 1548943)
    message(taxid[i], " ", length(keep), ' mapped metabolites')
    
    if(full.scored) {
      ## get taxid hierarchy. store as vector.  
      ## compare to each out taxid vector by lca
      message(" -- calculating similarities", '\n')
      taxid.row <- tax.match[1,1]
      taxid.column <- as.integer(tax.match[1,2])
      taxid.dt <- taxid.hierarchy[taxid.row, 1:ncol(taxid.hierarchy)]
      taxid.vector <- as.numeric(as.vector(unlist(taxid.dt)))
      
      ## replace NA values with large numbers which will not match any taxids later on.
      ## useful since using %in% for matching.  incomparables option only available for match.
      taxid.na <- which(is.na(taxid.vector))
      if(length(taxid.na) > 0) {
        taxid.vector[taxid.na] <- max.taxid + 1:length(taxid.na)
      }

      ## tmp is the index to the correct cid.lca row for each pubchem row
      th.ind <- which(names(cid.lca) == "species"):ncol(cid.lca)
      # th.ind <- th.ind[taxid.column:length(th.ind)]
      
      
      ## create vector of indices on which to perform calculations
      ## all CID which are in pc.bio which are also in cid.lca (taxonomy mapped)
      lca.cid <- cid.lca$cid
      tmp <- list()
      current.iteration <- 1
      tmp[[current.iteration]] <- match(pc.bio$cid, lca.cid)
      lca.cid[tmp[[current.iteration]]] <- NA
      l.is.na <- length(which(is.na(lca.cid)))
      delta.is.na <- 1
      while(delta.is.na > 0) {
        current.iteration <- current.iteration + 1
        tmp[[current.iteration]] <- match(pc.bio$cid, lca.cid)
        lca.cid[tmp[[current.iteration]]] <- NA
        delta.is.na <- length(which(is.na(lca.cid))) - l.is.na
        l.is.na <- length(which(is.na(lca.cid)))
      }
      tmp <- data.frame(tmp)
      names(tmp) <- NULL
      # system.time(tmp <- lapply(pc.bio$cid, FUN = function(x) pc.bio$cid[x] %in% cid.lca$cid))
      # tmp <- match(pc.bio$cid, cid.lca$cid)
      ## only calculate similarities when there is at least one !NA value (would be in first column, if present)
      do.sim <- which(!is.na(tmp[,1]))
      j <- NULL
      cid.lca.df <- data.frame(cid.lca)
      results <- foreach::foreach(j = do.sim) %dopar% {
        tryCatch({
          # requireNamespace('data.table', quietly = TRUE)
          tmp.j <- tmp[j,]
          tmp.j <- tmp.j[!is.na(tmp.j)]
          mtchs <- sapply(1:length(tmp.j), FUN = function(k) {
            tmp.k <- cid.lca.df[tmp.j[k], th.ind]
            min(which(taxid.vector %in% tmp.k))
          })
          mtch.col <- min(mtchs)
          mtch.col
          #when it throws an error, the following block catches the error
        }, error = function(msg){
          stop("error on ", j, '\n')
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
      #  message(x, ' ')
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
      
    }
  }
  
  if(!full.scored) {
    message("keeping ", length(keep), " metabolites", '\n')
    out <- pc.bio[keep, ]
    if(nrow(out) == 0) {
      stop("no metabolites found for taxid(s):", paste0(taxid, collapse = ", "))
    }
  } else {
    
    ## aggregate taxonomy.lca.similarity scores using assigned function
    
    use.cols <- which(grepl("taxonomy.lca.similarity.", names(out)))
    ..use.cols <- use.cols
    suppressWarnings(agg.sim <- apply(out[, ..use.cols, drop = FALSE], 1, aggregation.function, na.rm = TRUE, simplify = TRUE))
    agg.sim[is.infinite(agg.sim)] <- NA
    # agg.sim <- agg.sim[unique(c(which(is.infinite(agg.sim)), which(is.na(agg.sim))))] <- NA
    out[,paste0("taxonomy.lca.similarity.", "aggregate")] <- agg.sim
  }
  
  if(keep.scored.only) {
    if(any(names(out) == "taxonomy.lca.similarity.aggregate")) {
      keep <- !is.na(out$taxonomy.lca.similarity.aggregate)
      out <- out[keep,]
    }
  }
  
  if(get.properties) {
    message(" - calclulating rcdk properties ",  format(Sys.time()), '\n')
    cid.list <- as.list(out$cid)
    sm.list <- as.list(out$smiles)
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
    
    results.df <- do.call("rbind", results)
    out <- out[order(out$cid),]
    results.df <- results.df[order(as.numeric(row.names(results.df))),]
    
    out <- data.frame(
      out,
      results.df
    )
  }
  

  if(!is.null(out.dir)) {
    save(out, file = paste0(out.dir, "/", db.name, ".Rdata"))
  }
  
  return(out)
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
