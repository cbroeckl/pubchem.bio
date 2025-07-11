#' build.taxonomy.similarity.scores
#'
#' utilizes downloaded and properly formatted local pubchem data created by 'get.pubchem.ftp' and 'build.pubchem.bio' functions to append taxid-lca similiarty scores for each CID from the dataset created by 'build.pubchem.bio' function
#' @details utilizes downloaded and properly formatted local pubchem data created by 'get.pubchem.ftp' function
#' @param pc.directory directory from which to load pubchem .Rdata files
#' @param taxid integer vector of integer NCBI taxonomy IDs.  i.e.  c(9606, 1425170 ) for Homo sapiens and Homo heidelbergensis.    
#' @param get.properties logical. if TRUE, will return rcdk calculated properties:  XLogP, TPSA, HBondDonorCount and HBondAcceptorCount.
#' @param threads integer. how many threads to use when calculating rcdk properties.  parallel processing via DoParallel and foreach packages.  
#' @param aggregation.function character or function.  how should similarity scores be aggregated if multiple taxids are provided?  'mean', 'median', 'max', 'min' are all options, for example.  Default = 'max'. can also use a custom function.
#' @return a data frame containing pubchem CID ('cid') for all pubchem.bio (pc.bio) compounds, with additional data including one 'taxonomy.lca.similarity' column for each taxid provided.
#' @author Corey Broeckling
#' 
#' @export
#' 
build.taxonomy.similarity.scores <- function(
    pc.directory = NULL,
    taxid = c(),
    get.properties = TRUE,
    threads = 8,
    aggregation.function = 'max'
) {
  
  if(length(taxid) == 0) {
    stop("please list at least one integer taxid, i.e. 'taxid = c(4071, 4081)'", '\n')
  }
  
  if(file.exists(paste0(pc.directory, "/cid.lca.Rdata"))) {
    load(paste0(pc.directory, "/cid.lca.Rdata"))
  } else {
    error(paste0(pc.directory, "/cid.lca.Rdata"), "does not exist", '\n')
  }
  
  if(file.exists(paste0(pc.directory, "/taxid.heirarchy.Rdata"))) {
    load(paste0(pc.directory, "/taxid.heirarchy.Rdata"))
  } else {
    error(paste0(pc.directory, "/taxid.heirarchy.Rdata"), "does not exist", '\n')
  }
  
  if(file.exists(paste0(pc.directory, "/pc.bio.Rdata"))) {
    load(paste0(pc.directory, "/pc.bio.Rdata"))
  } else {
    error(paste0(pc.directory, "/pc.bio.Rdata"), "does not exist", '\n')
  }
  
  cid.lca <- cid.lca
  taxid.heirarchy <- taxid.heirarchy
  pc.bio <- pc.bio
  
  # data.table::setkey(d, "cid")
  # tmp <- match(d$cid, cid.preferred$cid)
  # use <- which(!is.na(tmp))
  data.table::setkey(cid.lca, "cid")
  
  for(i in 1:length(taxid)) {
    taxid.match <- which(taxid.heirarchy == taxid[i], arr.ind = TRUE)
    taxid.row <- taxid.match[1,1]
    taxid.column <- taxid.match[1,2]
    taxid.vector <- as.vector(unlist(taxid.heirarchy[taxid.row, taxid.column:ncol(taxid.heirarchy)]))
    tmp <- match(pc.bio$cid, cid.lca$cid)
    lca <- cid.lca$lca[tmp]
    
    rm(tmp); gc()
    # position of lca in target taxid heirarchy
    # probably will not need this - use simiarity = 1 for all metabolites
    # which map to target taxon id, as for 'build.taxon.metabolome'.
    tax.lca.sim.1 <- unlist(sapply((1:length(lca)), FUN = function(x) {
      if(is.na(lca[x])) {
        NA
      } else {
        mtch.col <- as.vector(taxid.vector == which(lca[x])[1])
        if(is.na(mtch.col)) {
          mtch.col <- ncol(taxid.heirarchy)
        }
        mtch.col
      }
      
    }, simplify = TRUE))
    

    
    # position of lca from root
    tax.lca.sim.2 <- match(lca, cid.lca$lca)
    
    tax.lca.sim.2 <- gsub('[[:digit:]]+', '', tax.lca.sim.2)
    tax.lca.sim.2 <- match(tax.lca.sim.2, names(taxid.heirarchy))
    tax.lca.sim.2 <- max(tax.lca.sim.2, na.rm = TRUE) - tax.lca.sim.2
    tax.lca.sim[is.na(tax.lca.sim)] <- ncol(taxid.heirarchy)
    tax.lca.sim <- 1-(tax.lca.sim - min(tax.lca.sim))/(max(tax.lca.sim)-min(tax.lca.sim))

    pc.bio[,paste0("taxonomy.lca.similarity.", taxid[i])] <- tax.lca.sim
    rm(tax.lca.sim)
  }
  
  ## aggregate to a single score, if multiple taxids provides
  if(length(taxid) > 1) {
    sim.cols <- grep("taxonomy.lca.similarity", colnames(pc.bio))
    taxonomy.lca.similarity <- apply(pc.bio[,sim.cols], 1, FUN = aggregation.function)
    pc.bio[,"taxonomy.lca.similarity.aggregate"] <- taxonomy.lca.similarity
  }
  
  out <- pc.bio
  rm(pc.bio)
  gc()
  
  ## create stable dopar function:
  `%dopar%` <- foreach::`%dopar%`
  
  if(get.properties) {
    cat(" - calclulating rcdk properties",  format(Sys.time()), '\n')
    cid.list <- as.list(out$cid)
    sm.list <- as.list(out$smiles)
    doParallel::registerDoParallel(cl <- parallel::makeCluster(threads))
    results <- foreach::foreach(i = 1:(length(cid.list))) %dopar% {
      desc <- c(
        "org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.AcidicGroupCountDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.BasicGroupCountDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor"
      )
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
    parallel::stopCluster(cl)
  }
  
  return(out)
  
}

# pc.bio.tax.scores <- build.taxonomy.similarity.scores(taxid = 4081, pc.directory = "C:/Temp/20250703", get.properties = FALSE)
