#' build.cid.lca
#'
#' utilizes downloaded and properly formatted local pubchem data created by 'get.pubchem.ftp' as input to generate a relationship between pubchem CID and the lowest common ancestor NCBI taxid
#' @details utilizes downloaded and properly formatted local pubchem data created by 'get.pubchem.ftp' function
#' @param pc.directory directory from which to load pubchem .Rdata files
#' @param tax.sources character vector.  containing taxonmy source values (from NCBI taxonomy database). i.e. c("LOTUS - the natural products occurrence database", "The Natural Products Atlas").  if NULL (default), use all values.  
#' @param use.pathways logical.  default = TRUE, should pathway data be used in building lowest common ancestor, when taxonomy is associated with a pathway?
#' @param use.conserved.pathways logical. default = FALSE, should 'conserved' pathways be used?  when false, only pathways with an assigned taxonomy are used. 
#' @param threads integer.  number of threads to use when finding lowest common ancestor.  parallel processing via DoParallel and foreach packages.   
#' @return a data frame containing pubchem CID ('cid'), and lowest common ancestor ('lca') NCBI taxonomy ID integer. will also save to pc.directory as .Rdata file.
#' @author Corey Broeckling
#' 
#' @export 
#' @importFrom foreach '%dopar%'
#' 
#' @param tax.sources vector. which taxonomy sources should be used?  defaults to c("LOTUS - the natural products occurrence database", "The Natural Products Atlas", "KNApSAcK Species-Metabolite Database", "Natural Product Activity and Species Source (NPASS)").


build.cid.lca <- function(
    pc.directory = 'R:/RSTOR-PMF/Software/db/met.db/20241216/',
    tax.sources = "LOTUS - the natural products occurrence database",
    use.pathways = TRUE,
    use.conserved.pathways = FALSE,
    threads = 8
) {
  
  ## load necessary files
  load(paste0(pc.directory, "/cid.taxid.Rdata"))
  data.table::setkey(cid.taxid, "cid")
  load(paste0(pc.directory, "/taxid.heirarchy.Rdata"))
  taxid.heirarchy <- taxid.heirarchy
  cat(" -" , nrow(cid.taxid), "taxonomy-cid associations from cid.taxid.Rdata file", '\n')
  
  # load(paste0(pc.directory, "/cid.pwid.Rdata"))
  
  if(!is.null(tax.sources)) {
    if(tax.sources[1] == "interactive") {
      source <- unique(cid.taxid$data.source)
      source.number <- 1:length(source)
      for(i in 1:length(source)) {
        cat(paste(source.number[i], source[i], '\n'))
      }
      use <- readline("enter source.number values for all sources, separated by a space:  ")
      use <- sort(as.numeric(unlist(strsplit(use, " "))))
      tax.sources <- source[use]
    } 
    cid.taxid <- cid.taxid[cid.taxid$data.source %in% tax.sources]
  }
  cat(" -" , nrow(cid.taxid), "taxonomy-cid associations after filtering by source", '\n')
  
  if(use.pathways) {
    load(paste0(pc.directory, "/cid.pwid.Rdata"))
    cid.pwid <- cid.pwid
    data.table::setkey(cid.pwid, "cid")
    sp.spec <- which(!is.na(cid.pwid$taxid))
    
    cid.taxid.2 <- data.frame(
      taxid = cid.pwid$taxid[sp.spec],
      cid = cid.pwid$cid[sp.spec], 
      data.source = cid.pwid$source[sp.spec]
    )
    cid.taxid.2 <- data.table::as.data.table(cid.taxid.2)
    cid.taxid.2 <- cid.taxid.2[!duplicated(cid.taxid.2), ]
    
    ## conserved pathways
    
    if(use.conserved.pathways) {
      con.path <- which(cid.pwid$pwtype == 'conserved')
      cids <- unique(cid.pwid$cid[con.path])
      taxids <- 33090
      cid.taxid.3 <- expand.grid(taxid = taxids, cid = cids, data.source = cid.pwid$source[con.path])
      cid.taxid.2 <- rbind(cid.taxid.2, cid.taxid.3)
      rm(cid.taxid.3)
      rm(taxids)
      rm(cids)
      rm(con.path)
      gc()
    }
    
    cid.taxid <- rbind(
      cid.taxid,
      cid.taxid.2
    )
    
    rm(cid.taxid.2)
    rm(sp.spec)
    rm(cid.pwid)
    gc()
    
    dups <- duplicated(cid.taxid[,1:2])
    cid.taxid <- cid.taxid[!duplicated(cid.taxid), ]
    cat(" -" , nrow(cid.taxid), "taxonomy-cid associations after adding pathway data", '\n')
  }
  
  cid <- table(cid.taxid$cid)
  n.tax <- as.vector(cid)
  cid <- as.numeric(names(cid))
  lca <-  vector(mode = 'integer', length(cid))
  th.mat <- as.matrix(taxid.heirarchy)
  th.vec <- as.vector(th.mat)
  th.vec <- data.table::data.table(
    "taxid" = as.integer(th.vec)
  )
  th.convert <- rep(0, nrow(th.mat))
  for(i in 2:ncol(th.mat)) {
    th.convert <- c(th.convert, rep(((i-1)*nrow(th.mat)), nrow(th.mat)))
  }
  
  cat(" -" , "findling lowest common ancestor for each cid", '\n')
  
  # # for each cid find lca
  # # return cid, lca vector
  # for(i in 1:length(cid)) {
  #   # for(i in 1:100) {
  #   taxids <- unique(cid.taxid$taxid[cid.taxid$cid == cid[i]])
  #   if(length(taxids) == 1) {
  #     lca[i] <- taxids
  #     next
  #   }
  #   
  #   ## convert multicolumn data table to single column data table so single column matching
  #   ## then convert single column rows back to multicolumn row index
  #   mtch <- which(th.vec$taxid %in% taxids)
  #   ## back convert to row numbers of original matrix
  #   tar.rows <- mtch - th.convert[mtch]
  #   
  #   sub.taxid.heirarchy <- taxid.heirarchy[sort(unique(tar.rows)),]
  #   sub.taxid.heirarchy <- data.frame(sub.taxid.heirarchy)
  #   for(j in 1:ncol(sub.taxid.heirarchy)) {
  #     taxids <- (unique(sub.taxid.heirarchy[,j]))
  #     if(length(taxids) == 1) {
  #       if(is.na(taxids[1])) {next}
  #       lca[i] <- sub.taxid.heirarchy[1,j]
  #       break
  #     }
  #   }
  # }
  # 
  # cid.lca <- data.table::data.table(
  #   cid,
  #   lca
  # )
  
  ## create foreach loop intead
  cid.list <- as.list(cid)
  Sys.time()
  # threads = 8
  doParallel::registerDoParallel(cl <- parallel::makeCluster(threads))
  results <- foreach::foreach(i = 1:(length(cid.list))) %dopar% {
    taxids <- unique(cid.taxid$taxid[cid.taxid$cid == cid.list[[i]]])
    out <- c(cid.list[[i]], NA)
    if(length(taxids) == 1) {
      out[2] <- taxids
    } else {
      ## convert multicolumn data table to single column data table so single column matching
      ## then convert single column rows back to multicolumn row index
      mtch <- which(th.vec$taxid %in% taxids)
      ## back convert to row numbers of original matrix
      tar.rows <- mtch - th.convert[mtch]
      
      sub.taxid.heirarchy <- taxid.heirarchy[sort(unique(tar.rows)),]
      sub.taxid.heirarchy <- data.frame(sub.taxid.heirarchy)
      for(j in 1:ncol(sub.taxid.heirarchy)) {
        taxids <- (unique(sub.taxid.heirarchy[,j]))
        if(length(taxids) == 1) {
          if(is.na(taxids[1])) {next}
          out[2] <- sub.taxid.heirarchy[1,j]
          break
        }
      }
    }
    return(out)
  }
  parallel::stopCluster(cl)
  Sys.time()
  cid.lca <- do.call("rbind", results)
  dimnames(cid.lca)[[2]] <- c("cid", "lca")
  cid.lca <- cid.lca[order(cid.lca[,1]),]
  cid.lca <- data.table::data.table(cid.lca)
  data.table::setkey(cid.lca, "cid")
  save(cid.lca, file = paste0(pc.directory, "/cid.lca.Rdata"))
  
  return(cid.lca)
}




