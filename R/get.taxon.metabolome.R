#' get.taxon.metabolome
#'
#' utilizes downloaded and properly formatted local pubchem data created by 'get.pubchem.ftp' function to filter a dataset created by 'build.pubchem.bio' function
#' @details utilizes downloaded and properly formatted local pubchem data created by 'get.pubchem.ftp' function
#' @param pc.directory directory from which to load pubchem .Rdata files
#' @param taxid integer vector of integer NCBI taxonomy IDs.  i.e.  c(9606, 1425170 ) for Homo sapiens and Homo heidelbergensis.    
#' @param get.properties logical. if TRUE, will return rcdk calculated properties:  XLogP, TPSA, HBondDonorCount and HBondAcceptorCount.
#' @param threads integer. how many threads to use when calculating rcdk properties.  parallel processing via DoParallel and foreach packages.  
#' @return a data frame containing pubchem CID ('cid'), and lowest common ancestor ('lca') NCBI taxonomy ID integer. will also save to pc.directory as .Rdata file.
#' @author Corey Broeckling
#' 
#' @export
#' 
get.taxon.metabolome <- function(
    pc.directory = NULL,
    taxid = c(),
    get.properties = TRUE,
    threads = 4
) {
  
  if(length(taxid) == 0) {
    stop("please list at least one integer taxid, i.e. 'taxid = c(4071, 4081)'", '\n')
  }
  
  load(paste0(pc.directory, "/cid.lca.Rdata"))
  cid.lca <- cid.lca
  load(paste0(pc.directory, "/taxid.heirarchy.Rdata"))
  taxid.heirarchy <- taxid.heirarchy
  
  if(!all(names(cid.lca) == c("cid", "lca"))) {
    stop("cid.lca object does not look correct, column names should be exactly 'cid' and 'lca'", '\n')
  }
  
  metabolome <- vector(mode = "integer", length = 0)
  for(i in 1:length(taxid)) {
    taxid.i <- taxid[i]
    tax.match <- which(taxid.heirarchy == taxid.i, arr.ind = TRUE)
    taxid.v <- as.vector(t(data.frame(taxid.heirarchy[tax.match[1,1],])))
    metabolome <- unique(c(metabolome, cid.lca$cid[cid.lca$lca %in% taxid.v]))
  }
  
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
  
  return(metabolome)
  
  
  
}
