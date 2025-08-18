#' build.primary.metabolome
#'
#' utilizes downloaded and properly formatted local pubchem data created by 'get.pubchem.ftp' function to filter a dataset created by 'build.pubchem.bio' function
#' @details utilizes downloaded and properly formatted local pubchem data created by 'get.pubchem.ftp' function
#' @param pc.directory directory from which to load pubchem .Rdata files
#' @param get.properties logical. if TRUE, will return rcdk calculated properties:  XLogP, TPSA, HBondDonorCount and HBondAcceptorCount.
#' @param threads integer. how many threads to use when calculating rcdk properties.  parallel processing via DoParallel and foreach packages.  
#' @param db.name character. what do you wish the file name for the saved version of this database to be?  default = 'primary.metabolome.'  Saved as an .Rdata file in the 'pc.directory' location. 
#' @param rcdk.desc vector. character vector of valid rcdk descriptors.  default = rcdk.desc <- c("org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor", "org.openscience.cdk.qsar.descriptors.molecular.AcidicGroupCountDescriptor", "org.openscience.cdk.qsar.descriptors.molecular.BasicGroupCountDescriptor", "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor"). To see descriptor categories: 'dc <- rcdk::get.desc.categories(); dc' .  To see the descriptors within one category: 'dn <- rcdk::get.desc.names(dc\[4\]); dn'. Note that the four default parameters are relatively fast to calculate - some descriptors take a very long time to calculate.  you can calculate as many as you wish, but processing time will increase the more descriptors are added.   
#' @param pubchem.bio.object R data.table, generally produced by build.pubchem.bio; preferably, define pc.directory
#' @param output.directory directory to which the pubchem.bio database is saved.  If NULL, will try to save in pc.directory (if provided), else not saved. 
#' @param keep.primary.only logical.  If TRUE, only biological metabolites scored as 'primary' are returned. If FALSE, full dataset of metabolites is returned, with new logical column, 'primary' 
#' @param min.tax.ct integer.  if assigned an integer value, only those metabolites with at least min.tax.ct unique taxonomy assigments are considered 'primary'.  default = 3. 
#' @return a data frame containing pubchem CID ('cid'), and lowest common ancestor ('lca') NCBI taxonomy ID integer. will also save to pc.directory as .Rdata file.
#' @examples
#' data('cid.lca', package = "pubchem.bio")
#' data('pubchem.bio', package = "pubchem.bio")
#' data('taxid.hierarchy', package = "pubchem.bio")
#' my.taxon.db <- pubchem.bio::build.taxon.metabolome(
#' pubchem.bio.object = pubchem.bio,
#' get.properties = FALSE, threads = 1, taxid = c(1))
#' head(my.taxon.db)
#' @author Corey Broeckling
#' 
#' @export
#' 
#' 
build.primary.metabolome <- function(
    pc.directory = NULL,
    get.properties = FALSE,
    threads = 8,
    db.name = "primary.metabolome", 
    rcdk.desc = c(
      "org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor",
      "org.openscience.cdk.qsar.descriptors.molecular.AcidicGroupCountDescriptor",
      "org.openscience.cdk.qsar.descriptors.molecular.BasicGroupCountDescriptor",
      "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor"
    ),
    pubchem.bio.object = NULL,
    output.directory = NULL,
    keep.primary.only = TRUE,
    min.tax.ct = 3
) {
  
  out.dir <- pc.directory
  if(is.null(out.dir)) out.dir <- output.directory
  
  if(is.null(pc.directory) & is.null(pubchem.bio.object)) {
    stop("if you opt to note define the pc.directory, you must provide ALL of 'pubchem.bio.object', 'taxid.hierarchy.object', 'pubchem.bio.object' variables", '\n')
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
  
  keep <- (out$lca == 1 & out$taxonomy.ct >= min.tax.ct)
  
  if(keep.primary.only) {
    out <- out[keep,]
  } else {
    primary <- rep(FALSE, nrow(out))
    primary[keep] <- TRUE
    out$primary <- primary
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
  
  
  
  return(out)
  
  if(!is.null(out.dir)) {
    save(out, file = paste0(out.dir, "/", db.name, ".Rdata"))
  }
  
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
