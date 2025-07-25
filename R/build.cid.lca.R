#' build.cid.lca
#'
#' utilizes downloaded and properly formatted local pubchem data created by 'get.pubchem.ftp' as input to generate a relationship between pubchem CID and the lowest common ancestor NCBI taxid
#' @details utilizes downloaded and properly formatted local pubchem data created by 'get.pubchem.ftp' function
#' @param pc.directory directory from which to load pubchem .Rdata files. alternatively provide cid.taxid.object, taxid.hierarchy.object, and cid.pwid.object as data.table R objects. 
#' @param tax.sources vector. which taxonomy sources should be used?  defaults to c("LOTUS - the natural products occurrence database", "The Natural Products Atlas", "KNApSAcK Species-Metabolite Database", "Natural Product Activity and Species Source (NPASS)").
#' @param use.pathways logical.  default = TRUE, should pathway data be used in building lowest common ancestor, when taxonomy is associated with a pathway?
#' @param use.conserved.pathways logical. default = FALSE, should 'conserved' pathways be used?  when false, only pathways with an assigned taxonomy are used. 
#' @param threads integer.  number of threads to use when finding lowest common ancestor.  parallel processing via DoParallel and foreach packages.   
#' @param cid.taxid.object R data.table, generally produced by get.pubchem.ftp; alternatively, define pc.directory
#' @param taxid.hierarchy.object R data.table, generally produced by get.pubchem.ftp; alternatively, define pc.directory
#' @param cid.pwid.object R data.table, generally produced by get.pubchem.ftp; alternatively, define pc.directory
#' @param output.directory directory to which the pubchem.bio database is saved.  If NULL, will try to save in pc.directory (if provided). If both directories are NULL, not saved, only returned as in memory 
#' @return a data frame containing pubchem CID ('cid'), and lowest common ancestor ('lca') NCBI taxonomy ID integer. will also save to pc.directory as .Rdata file.
#' @author Corey Broeckling
#' 
#' @examples
#' data('cid.taxid', package = "pubchem.bio")
#' data('taxid.hierarchy', package = "pubchem.bio")
#' cid.lca.out <- build.cid.lca(tax.sources =  "LOTUS - the natural products occurrence database", 
#' use.pathways = FALSE, 
#' threads = 1, cid.taxid.object = cid.taxid, 
#' taxid.hierarchy.object = taxid.hierarchy)
#' head(cid.lca.out)
#' @export 
#' @importFrom foreach '%dopar%'
#' 


build.cid.lca <- function(
    pc.directory = NULL,
    tax.sources = "LOTUS - the natural products occurrence database",
    use.pathways = TRUE,
    use.conserved.pathways = FALSE,
    threads = 8,
    cid.taxid.object = NULL,
    taxid.hierarchy.object = NULL,
    cid.pwid.object = NULL,
    output.directory = NULL
) {
  
  out.dir <- pc.directory
  if(is.null(out.dir)) out.dir <- output.directory
  
  if(is.null(pc.directory) & is.null(cid.taxid.object)) {
    stop("if you opt to note define the pc.directory, you must provide ALL of 'cid.taxid.object', 'taxid.hierarchy.object', 'cid.pwid.object' variables", '\n')
  }
  
  ## load necessary files
  
  if(!is.null(pc.directory)) {
    
    load(paste0(pc.directory, "/cid.taxid.Rdata"))
    data.table::setkey(cid.taxid, "cid")
    cid.taxid <- cid.taxid
    
    load(paste0(pc.directory, "/taxid.hierarchy.Rdata"))
    taxid.hierarchy <- taxid.hierarchy
    
    if(use.pathways) {
      load(paste0(pc.directory, "/cid.pwid.Rdata"))
      cid.pwid <- cid.pwid
      data.table::setkey(cid.pwid, "cid")
    }
    
  } else {
    cid.taxid = cid.taxid.object
    taxid.hierarchy = taxid.hierarchy.object
    if(use.pathways) {
      cid.pwid = cid.pwid.object
    }
    
  }
  
  message(" - " , nrow(cid.taxid), " taxonomy-cid associations from cid.taxid.Rdata file", '\n')
  
  if(!base::is.null(tax.sources)) {
    
    if(tax.sources[1] == "interactive") {
      source <- unique(cid.taxid$data.source)
      source.number <- 1:length(source)
      for(i in 1:length(source)) {
         message(paste(source.number[i], source[i], '\n'))
      }
      use <- readline("enter source.number values for all sources, separated by a space:  ")
      use <- sort(as.numeric(unlist(strsplit(use, " "))))
      tax.sources.internal <- source[use]
    }
    tax.sources.internal <- tax.sources
    keep <- cid.taxid$data.source %in% tax.sources.internal
    cid.taxid <- cid.taxid[keep,]
  }

  
  message(" - " , nrow(cid.taxid), " taxonomy-cid associations after filtering by source", '\n')

  
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
    message(" - " , nrow(cid.taxid), " taxonomy-cid associations after adding pathway data", '\n')
  }
  
  cid <- table(cid.taxid$cid)
  n.tax <- as.vector(cid)
  cid <- as.numeric(names(cid))
  lca <-  vector(mode = 'integer', length(cid))
  th.mat <- as.matrix(taxid.hierarchy)
  th.vec <- as.vector(th.mat)
  th.vec <- data.table::data.table(
    "taxid" = as.integer(th.vec)
  )
  th.convert <- rep(0, nrow(th.mat))
  for(i in 2:ncol(th.mat)) {
    th.convert <- c(th.convert, rep(((i-1)*nrow(th.mat)), nrow(th.mat)))
  }
  
  message(" - " , "finding lowest common ancestor for each cid", '\n')
  
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
  #   sub.taxid.hierarchy <- taxid.hierarchy[sort(unique(tar.rows)),]
  #   sub.taxid.hierarchy <- data.frame(sub.taxid.hierarchy)
  #   for(j in 1:ncol(sub.taxid.hierarchy)) {
  #     taxids <- (unique(sub.taxid.hierarchy[,j]))
  #     if(length(taxids) == 1) {
  #       if(is.na(taxids[1])) {next}
  #       lca[i] <- sub.taxid.hierarchy[1,j]
  #       break
  #     }
  #   }
  # }
  # 
  # cid.lca <- data.table::data.table(
  #   cid,
  #   lca
  # )
  
  ## create stable dopar function:
  `%dopar%` <- foreach::`%dopar%`
  
  ## create foreach loop instead
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
      
      sub.taxid.hierarchy <- taxid.hierarchy[sort(unique(tar.rows)),]
      sub.taxid.hierarchy <- data.frame(sub.taxid.hierarchy)
      for(j in 1:ncol(sub.taxid.hierarchy)) {
        taxids <- (unique(sub.taxid.hierarchy[,j]))
        if(length(taxids) == 1) {
          if(is.na(taxids[1])) {next}
          out[2] <- sub.taxid.hierarchy[1,j]
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
  
  ## remove any rows with lca = NA.  i think these primarily derive from taxa-lca relationships in which
  ## the rank is not in our taxid.hierarchy levels (i.e. subspecies)
  if(any(is.na(cid.lca$lca))) {
    cid.lca <- cid.lca[!is.na(cid.lca$lca),]
  }
  
  
  ## and now record the column location in the taxid.hierarchy for each lca
  ## turn taxid.hierarchy into a vector for faster matching
  th.vec <- unlist(taxid.hierarchy)
  
  # position of lca in column position
  lca <- cid.lca$lca
  # lca[is.na(lca)] <- max(th.vec, na.rm = TRUE)+max(lca, na.rm = TRUE) + 1
  taxid.dt <- data.table::data.table('taxid' = th.vec, 'level' = names(th.vec))
  # data.table::setkey(taxid.dt, "taxid")
  tmp <- match(lca, taxid.dt$taxid)
  tmp <- names(th.vec)[tmp]
  tmp <- gsub('[[:digit:]]+', '', tmp)
  
  ## at the column number of each lca.  will be used for taxon metabolomics
  hierarchy.column <- match(tmp, names(taxid.hierarchy))
  rm(tmp); gc()
  cid.lca$hierarchy.column <- hierarchy.column
  
  ## add taxid.hierarchy to cid.lca, with values below the hierarchy column set to NA. will be used for taxon metabolome
  tmp <- match(lca, taxid.dt$taxid)
  tmp <- names(th.vec)[tmp]
  tmp <- gsub('[[:alpha:]]+', '', tmp)
  tmp <- as.numeric(tmp)
  if(any(is.na(tmp))) {
    tmp.rm <- which(is.na(tmp))
    cid.lca <- cid.lca[-tmp.rm,]
    tmp <- tmp[-tmp.rm]
  }
  cid.lca.h <- taxid.hierarchy[tmp,]
  for(i in (1:nrow(cid.lca))) {
    cid.lca.h[i, 1:max((cid.lca$hierarchy.column[i]-1),1)] <- NA
  }
  
  cid.lca <- data.frame(
    cid.lca,
    cid.lca.h
  )

  
  if(!is.null(out.dir)) {
    save(cid.lca, file = paste0(out.dir, "/cid.lca.Rdata"))
  }
  
  
  return(cid.lca)
}





