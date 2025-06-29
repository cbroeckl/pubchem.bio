#' tax.filter.pubchem.bio
#'
#' utilizes downloaded and properly formatted local pubchem data created by 'get.pubchem.ftp' function to filter a dataset created by 'build.pubchem.bio' function
#' @details utilizes downloaded and properly formatted local pubchem data created by 'get.pubchem.ftp' function
#' @param pc.directory directory from which to load pubchem .Rdata files
#' @param tax.sources character vector.  containing taxonmy source values (from NCBI taxonomy database). i.e. c("LOTUS - the natural products occurrence database", "The Natural Products Atlas").  if NULL (default), use all values.  
#' @param use.pathways character vector.  default = NULL, all sources used.  can set to 'interactive' to enable you to choose the sources at the console.  or set to source name vector, i.e. c("LOTUS - the natural products occurrence database", "The Natural Products Atlas").
#' @param use.parent.cid logical. should CIDs be replaced with parent CIDs? 
#' @param threads integer. how many threads to use when calculating lowest common ancestor.  parallel processing via DoParallel and foreach packages.  
#' @return a data frame containing pubchem CID ('cid'), and lowest common ancestor ('lca') NCBI taxonomy ID integer. will also save to pc.directory as .Rdata file.
#' @author Corey Broeckling
#' 
#' @export 

build.taxid.lca <- function(
    pc.directory = NULL,
    tax.sources = NULL,
    use.pathways = TRUE,
    use.parent.cid = TRUE,
    threads = 8
) {
  
  ## load necessary files
  load(paste0(pc.directory, "/cid.taxid.Rdata"))
  cid.taxid <- cid.taxid
  load(paste0(pc.directory, "/taxid.heirarchy.Rdata"))
  taxid.heirarchy <- taxid.heirarchy
  
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
  }
  
  if(use.pathways) {
    load(paste0(pc.directory, "/cid.pwid.Rdata"))
    cid.pwid <- cid.pwid
    sp.spec <- which(!is.na(cid.pwid$taxid))
    cid.taxid.2 <- data.frame(
      taxid = cid.pwid$taxid[sp.spec],
      cid = cid.pwid$cid[sp.spec],
      data.source = cid.pwid$source[sp.spec]
    )
    cid.taxid.2 <- data.table::as.data.table(cid.taxid.2)
    cid.taxid.2 <- cid.taxid.2[!duplicated(cid.taxid.2), ]
    
    ## conserved plantcyc pathways
    con.path <- which(cid.pwid$pwtype == 'conserved' & cid.pwid$source == "PlantCyc")
    cids <- unique(cid.pwid$cid[con.path])
    taxids <- 33090
    cid.taxid.3 <- expand.grid(taxid = taxids, cid = cids)
    cid.taxid.3 <- data.frame(cid.taxid.3, data.source = "PlantCyc")
    
    cid.taxid <- rbind(
      cid.taxid,
      cid.taxid.2,
      cid.taxid.3
    )
    
    rm(cid.taxid.2)
    rm(cid.taxid.3)
    rm(cids)
    rm(taxids)
    rm(con.path)
    rm(sp.spec)
    rm(cid.pwid)
    gc()
    
    dups <- duplicated(cid.taxid)
    cid.taxid <- cid.taxid[!duplicated(cid.taxid), ]
  }
  
  # <- data.frame(
  #   clade = c("root", "metazoa", "bacteria", "viridiplantae", "fungi", "archaea", "eukaryota", "mammalia"),
  #   level = c("root", "kingdom", "superkingdom", "kingdom", "kingdom", "superkingdom", "superkingdom", "class"),
  #   taxid = c(1, 33208, 2, 33090, 4751, 2157, 2759, 40674)
  # )
  
  
  cid <- table(cid.taxid$cid)
  n.tax <- as.vector(cid)
  cid <- as.numeric(names(cid))
  # lca <-  vector(mode = 'integer', length(cid))
  ## converted multicolumn data table to single column data table so single column matching
  th.mat <- as.matrix(taxid.heirarchy)
  th.vec <- as.vector(th.mat)
  # th.vec.ch <- as.character(th.vec)
  th.dt <- data.table::data.table(
    "taxid" = as.integer(th.vec)
  )
  th.convert <- rep(0, nrow(th.mat))
  for(i in 2:ncol(th.mat)) {
    th.convert <- c(th.convert, rep(((i-1)*nrow(th.mat)), nrow(th.mat)))
  }
  
  ## replace cid with parent
  if(use.parent.cid) {
    load(paste0(pc.directory, "/cid.parent.Rdata"))
    cid.parent <- cid.parent
    data.table::setkey(cid.parent, "cid")
    m <- match(cid, cid.parent$cid)
    parent.cid <- cid.parent$parent.cid[m]
    parent.cid[which(is.na(parent.cid))] <- cid[which(is.na(parent.cid))]
    cid <- parent.cid
    # cid <- sort(unique(cid))
    # cat(" - after replacing CID with parent CID, current unique cid count:" , length(cid), '\n')
    rm(parent.cid); rm(m); gc()
  }
  
  # for each cid find lca
  # return cid, lca vector
  # cid.l <- split(cid, ceiling(seq_along(cid)/20))
  # for(i in 1:length(cid)) {
  doParallel::registerDoParallel(cl <- parallel::makeCluster(threads))
  set.seed(123)
  # inds <- c(531, sample(1:length(cid), 5000))
  lca <- foreach::foreach(i = 1:length(cid)) %dopar% {
    i <- i
    # set.seed(123); system.time(for(i in sample(1:length(cid), 100)) {
    
    taxids <- cid.taxid$cid == cid[i]
    taxids <- cid.taxid$taxid[taxids]
    taxids <- unique(taxids)
    taxids <- taxids[!is.na(taxids)]
    if(length(taxids) > 1) {
      # taxids.ch <- as.character(taxids)
      # tmp <- sapply(1:length(taxids), FUN = function(x) which(th.mat == taxids[x], arr.ind = TRUE))
      ## converted multicolumn data table to single column data table so single column matching
      ## then convert single column rows back to multicolumn row index
      mtch <- th.dt$taxid %in% taxids # about 14 sec with %fin% for 100 
      # system.time(mtch <- fmatch(taxids, th.dt$taxid)) # 
      mtch <- which(mtch)
      ## back convert to row numbers of original matrix
      tar.rows <- mtch - th.convert[mtch]
      tar.rows <- collapse::funique(tar.rows, method = 3) ## 0.3 second for one item
      sub.taxid.heirarchy <- taxid.heirarchy[tar.rows,]
      sub.taxid.heirarchy <- data.frame(sub.taxid.heirarchy)
      all.equal <- sapply(1:ncol(sub.taxid.heirarchy), FUN = function(x) stats::sd(sub.taxid.heirarchy[,x]))  ## 0.01 second for one item
      if(all(is.na(all.equal))) {
        ret <- taxids[which(taxids %in% sub.taxid.heirarchy[1,])]
      } else {
        ret <- sub.taxid.heirarchy[1,min(which(all.equal == 0))]
      }
    }
    if(length(taxids) == 1) {ret <- taxids}

    return(ret)
  }
  parallel::stopCluster(cl)
  
  lca <- unlist(lca)
  
  cid.lca <- data.table::data.table(
    cid,
    lca
  )
  
  save(cid.lca, file = paste0(pc.directory, "/cid.lca.Rdata"))
  
  return(cid.lca)
}
