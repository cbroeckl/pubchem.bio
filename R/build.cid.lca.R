#' build.cid.lca
#'
#' utilizes downloaded and properly formatted local pubchem data created by 'get.pubchem.ftp' as input to generate a relationship between pubchem CID and the lowest common ancestor NCBI taxid
#' @details utilizes downloaded and properly formatted local pubchem data created by 'get.pubchem.ftp' function
#' @param pc.directory directory from which to load pubchem .Rdata files. alternatively provide cid.taxid.object, taxid.hierarchy.object, and cid.pwid.object as data.table R objects. 
#' @param tax.sources vector. which taxonomy sources should be used?  defaults to c("LOTUS - the natural products occurrence database", "The Natural Products Atlas", "KNApSAcK Species-Metabolite Database", "Natural Product Activity and Species Source (NPASS)").
#' @param use.pathways logical.  default = TRUE, should pathway data be used in building lowest common ancestor, when taxonomy is associated with a pathway?
#' @param use.conserved.pathways logical. default = FALSE, should 'conserved' pathways be used?  when false, only pathways with an assigned taxonomy are used. 
#' @param threads integer.  number of threads to use when finding lowest common ancestor.  parallel processing via DoParallel and foreach packages.   
#' @param cid.taxid.object R data.table, generally produced by get.pubchem.ftp; preferably, define pc.directory
#' @param taxid.hierarchy.object R data.table, generally produced by get.pubchem.ftp; preferably, define pc.directory
#' @param cid.pwid.object R data.table, generally produced by get.pubchem.ftp; preferably, define pc.directory
#' @param min.taxid.table.length integer.  when there are few taxa reported to synthesize a particular compound, and those few taxa are spread widely across biology, the LCA concept breaks down.  This value controls the decision as to whether to determine LCA within taxonomic ranks, rather within the full taxonomy hierarchy.  see details. 
#' @param output.directory directory to which the pubchem.bio database is saved.  If NULL, will try to save in pc.directory (if provided). If both directories are NULL, not saved, only returned as in memory 
#' @return a data frame containing pubchem CID ('cid'), and lowest common ancestor ('lca') NCBI taxonomy ID integer. will also save to pc.directory as .Rdata file.
#' @author Corey Broeckling
#' 
#' @examples
#' data('cid.taxid', package = "pubchem.bio")
#' data('taxid.hierarchy', package = "pubchem.bio")
#' data('cid.pwid', package = "pubchem.bio")
#' cid.lca.out <- build.cid.lca(
#' tax.sources =  "LOTUS - the natural products occurrence database",
#' use.pathways = FALSE, 
#' threads = 1, cid.taxid.object = cid.taxid,
#' taxid.hierarchy.object = taxid.hierarchy,
#' cid.pwid.object = cid.pwid)
#' head(cid.lca.out)
#' @export 
#' @importFrom foreach '%dopar%'
#' 
#' @details
#' Some metabolism is highly conserved - all species perform those reactions.  
#' Other metabolism is highly specific - there is one know species to produce 
#' that metabolite. Sometimes, it is in between.  The lowest common ancestor 
#' approach allows us to analyze these patterns and put them to use to 
#' generalize metabolites for metabolomics across species. 
#' 
#' Biology is more complex than that though.  Natural products are often 
#' reported as being synthesized by an organism which is in symbiosis with 
#' a second organism.  The taxonomic assignment is sometimes both organisms, 
#' even if neither would create that product in isolation, or if only one is 
#' actually capable of producing that metabolite.  In these situations, the 
#' LCA approach can break down.  For example, if a bacteria is in symbiosis 
#' with an algae, and each is listed as producing the metabolite, the LCA will 
#' be assigned as '1' - the root of all biology, since we have to go back to 
#' the base of the taxonomic tree to find the common taxonomic ancestor of 
#' prokaryotes and eukaryotes.  In this example, there are two unique species, 
#' genera, families, orders, etc listed in the full taxonomic
#' hierarchy for this metabolite.  
#' 
#' The 'min.unique.taxid.ct' variable controls 
#' sensitivity to this phenomenon in assigning LCA.  The number of unique taxa 
#' which are mapped to each metabolite varies by taxonomic level.  it may map 
#' to two species, but only one genus.  in that case, the genus is assigned as 
#' the LCA.  However, if the metabolite maps to two unique species, 
#' two unique genera, two unique families, two unique kingdoms, and one unique 
#' domain, we should ask ourselves whether this sparse patterns supports that 
#' this metabolite should be marked as conserved' or 'primary.'  What makes 
#' more intuitive sense is to conclude that there are may be extenuating 
#' circumstances which have resulted from unique biology.  For example, 
#' Ceratodictyol B is reported from _Haliclona cymaeformis_ and _Ceratodictyon_ 
#' _spongiosum_, one of which is a red algal symbiont of the other. At each 
#' taxonomic level, there are either 0, 1, or 2 unique taxonomy IDs. 0 unique 
#' levels is uninteresting - that just reflects that there is no taxonomy 
#' assigned for those lineages at that level.  
#' 
#' What is more interesting is the number of unique levels of the number of 
#' unique taxonomy ids.  in the case of Ceratodictyol B, the only other value is 
#' '2'.  There are 2 unique taxonomy IDs at each level species, genus, order, 
#' class, and phylum.  So there are five taxonomic levels that have exactly 2 
#' unique taxonomy IDs, and there are no taxonomic levels which have more than 2 
#' unique taxids.  We will call this the taxid.ct.table length, where the 
#' taxid.ct.table is the table of frequencies of the number of unique taxids at 
#' each taxonomic level.  the length is the number of unique values when 
#' IGNORING '0' or '1'.  When the taxid.ct.table length is less than or equal
#' to min.taxid.table.length, the lca is calcluated within the lowest taxonomic 
#' level that has the most frequent unique taxonomy ID count.  
#' 
#' For the Ceratodictyol B example, this would mean that we would find that '2' 
#' was the most common number of unique taxids reported, so we find that the 
#' lowest taxonomic level which reports two unique taxids is 'species'.  LCA is 
#' for assigned to those two species.  If however, there were two _Ceratodicyon_ 
#' spp reported, then the species level would have 3 unique taxids, and there 
#' would be 4 levels (rather than five) which have 2unique taxids.  the lowest 
#' taxonomic level with 2 unique taxids, the most frequent count observed, 
#' would now be 'genus', so LCA would be assigned for within each level of 
#' 'genus'.  This would mean that the first LCA would be assigned to the 
#' _Ceratodicyon_ genus, since there are multiple _Ceratodicyon_ species 
#' reported, and then a second LCA would be assigned to the _Haliclona_ 
#' _cymaeformis_ species.  Sorry it is so complicated.  Life is complicated.  
#' 
#' 
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
    min.taxid.table.length = 3,
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
        message(paste(source.number[i], " ",  source[i], '\n'))
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
  cid <- as.numeric(names(cid))
  # cid <- c(1, 2, 3, 5793, 1548943, 44254980, 19)
  # lca <-  vector(mode = 'integer', length(cid))
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
  
  ## create stable dopar function:
  `%dopar%` <- foreach::`%dopar%`
  
  
  cid.list <- as.list(cid)
  # .SD <- data.table::.SD
  
  # threads = 8
  cl <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  base::on.exit({
    parallel::stopCluster(cl)
    rm(cl)
  }) 
  results <- foreach::foreach(i = 1:(length(cid.list)), .errorhandling = "pass") %dopar% {
  # results <- foreach::foreach(i = 1:13000, .errorhandling = "pass") %dopar% {  
    
    requireNamespace('data.table')
    .datatable.aware = TRUE
    ..cols <- NULL
    
    taxids <- unique(cid.taxid$taxid[cid.taxid$cid == cid.list[[i]]])
    if(any(is.na(taxids))) {taxids <- taxids[!is.na(taxids)]}
    if(length(taxids) == 0) {
      out <- c(cid.list[[i]], NA)
      out <- matrix(out, ncol = 2, byrow = TRUE)
    }
    if(length(taxids) == 1) {
      out <- c(cid.list[[i]], taxids)
      out <- matrix(out, ncol = 2, byrow = TRUE)
    }
    if(length(taxids) > 1) {
      ## convert multicolumn data table to single column data table so single column matching
      ## then convert single column rows back to multicolumn row index
      mtch <- which(th.vec$taxid %in% taxids)
      if(length(mtch) == 0) {
        out <- rep(NA, 2)
        out <- matrix(out, ncol = 2, byrow = TRUE)
      } else {
        ## back convert to row numbers of original matrix
        tar.rows <- mtch - th.convert[mtch]
        ## subset the taxid for only those in containing our taxids of interest.
        sub.taxid.hierarchy <- taxid.hierarchy[sort(unique(tar.rows)),]
        
        ## remove any columsn with NA values
        has.na <- sapply(1:ncol(sub.taxid.hierarchy), FUN = function(x) any(is.na(sub.taxid.hierarchy[,x, with = FALSE])))
        sub.taxid.hierarchy <- sub.taxid.hierarchy[,!has.na, with = FALSE]
        
        ## determine the number of unique taxids at each taxonomic level. 
        unique.taxids.sub <- sapply(1:ncol(sub.taxid.hierarchy), FUN = function(x) {
          select_cols = names(sub.taxid.hierarchy)[x]
          length(unique(unlist(sub.taxid.hierarchy[ , select_cols, with = FALSE])))
        })
        na.taxids.sub <- as.integer(sapply(1:ncol(sub.taxid.hierarchy), FUN = function(x) {
          select_cols = names(sub.taxid.hierarchy)[x]
          any(is.na(sub.taxid.hierarchy[ , select_cols, with = FALSE]))
        }))
        
        unique.taxids.sub <- unique.taxids.sub - na.taxids.sub
        unique.taxids.sub.table <- table(unlist(unique.taxids.sub))
        # unique.taxids.sub.table <- unique.taxids.sub.table[!names(unique.taxids.sub.table) %in% c("0", "1")]
        unique.taxids.sub.table <- unique.taxids.sub.table[!names(unique.taxids.sub.table) %in% c("0")]
        # if(length(unique.taxids.sub.table) > 1 & any(names(unique.taxids.sub.table) == "1")) {
        #   unique.taxids.sub.table <- unique.taxids.sub.table[-which(names(unique.taxids.sub.table) == "1")]
        # }
        l.table <- length(unique.taxids.sub.table)
        if(l.table == 0) out <- c(cid.list[[i]], NA)
        if(is.na(l.table)) out <- c(cid.list[[i]], NA)
        if(l.table <= min.taxid.table.length) {
          unique.taxids.sub.table <- sort(unique.taxids.sub.table, decreasing = TRUE)
          tax.lev <- which(unique.taxids.sub == as.numeric(names(unique.taxids.sub.table)[1]))[1]
          sub.tax.ids <- unlist(sub.taxid.hierarchy[, tax.lev, with = FALSE])
          tax.lev.ids <- unique(sub.tax.ids)
          tax.lev.ids <- tax.lev.ids[which(!(is.na(tax.lev.ids)))]
          
          out <- matrix(nrow = length(tax.lev.ids), ncol = 2)
          out[,1] <- cid.list[[i]]
          for(j in 1:nrow(out)) {
            sub.tar.rows <- which(sub.tax.ids == tax.lev.ids[j])
            sub.sub.taxid.hierarchy <- sub.taxid.hierarchy[sub.tar.rows,]
            unique.taxids.sub.sub <- sapply(1:ncol(sub.sub.taxid.hierarchy), FUN = function(x) {
              select_cols = names(sub.sub.taxid.hierarchy)[x]
              length(unique(unlist(sub.sub.taxid.hierarchy[ , select_cols, with = FALSE])))
            })
            na.taxids.sub.sub <- as.integer(sapply(1:ncol(sub.sub.taxid.hierarchy), FUN = function(x) {
              select_cols = names(sub.sub.taxid.hierarchy)[x]
              any(is.na(sub.sub.taxid.hierarchy[ , select_cols, with = FALSE]))
            }))
            unique.taxids.sub.sub <- unique.taxids.sub.sub - na.taxids.sub.sub
            lca <- as.vector(unlist(sub.sub.taxid.hierarchy[1, which(unique.taxids.sub.sub == 1)[1], with = FALSE]))
            # for(k in 1:ncol(sub.taxid.hierarchy)) {
            #   sub.sub.taxids <- unique(unlist(sub.sub.taxid.hierarchy[, tax.lev, with = FALSE]))
            #   l.sub.sub.taxids <- length(sub.sub.taxids)
            #   if(is.na(l.sub.sub.taxids)) 
            #   if(l.sub.sub.taxids == 0) next
            #   if(length(sub.sub.taxids) == 1) {
            #     if(is.na(sub.sub.taxids[1])) {next}
            #     lca <- sub.sub.taxids
            #     break
            #   }
            # }
            if(length(lca)==1) {
              out[j, 2] <- as.numeric(lca)
            } else {
              out[j, 2] <- NA
            }
          }
        } else {
          # taxid.hierarchy.df <- data.frame(taxid.hierarchy)
          # percent.missing <- sapply(1:ncol(taxid.hierarchy), FUN = function(x) {length(which(is.na(taxid.hierarchy.df[,x])))})/nrow(taxid.hierarchy.df)
          out <- as.numeric(unlist(c(cid.list[[i]], sub.taxid.hierarchy[1, which(unique.taxids.sub == 1)[1], with = FALSE])))
          out <- matrix(out, ncol = 2, byrow = TRUE)
          out <- matrix(out, ncol = 2, byrow = TRUE)
        }
      }
    }
    return(out)
  }
  
  # results.old <- results
  # for(i in 1:length(results)) {
  #   results[[i]] <- as.numeric(results[[i]])
  # }
  # failed.at
  # cid.lengths <- sapply(1:length(results), FUN = function(x) length(results[[x]]))
  
  cid.lca <- do.call("rbind", results)
  if(any(is.na(cid.lca[,1]))) {
    cid.lca <- cid.lca[-which(is.na(cid.lca[,1])),]
  }
  
  # is.list(cid.lca[,1])
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
  # lca[is.na(lca)] <- max(th.vec, na.rm = TRUE)+max(lca, na.rm = TRUE) + 1
  taxid.dt <- data.table::data.table('taxid' = th.vec, 'level' = names(th.vec))
  # data.table::setkey(taxid.dt, "taxid")
  tmp <- match(cid.lca$lca, taxid.dt$taxid)
  tmp <- names(th.vec)[tmp]
  tmp <- gsub('[[:digit:]]+', '', tmp)
  
  ## at the column number of each lca.  will be used for taxon metabolomics
  lca.level <- match(tmp, names(taxid.hierarchy))
  rm(tmp); gc()
  cid.lca$lca.level <- lca.level
  
  ## add taxid.hierarchy to cid.lca, with values below the hierarchy column set to NA. will be used for taxon metabolome
  tmp <- match(cid.lca$lca, taxid.dt$taxid)
  tmp <- names(th.vec)[tmp]
  for(i in 1:length(names(taxid.hierarchy))) {
    tmp <- gsub(names(taxid.hierarchy)[i], '', tmp)
  }
  
  tmp <- as.numeric(tmp)
  if(any(is.na(tmp))) {
    tmp.rm <- which(is.na(tmp))
    cid.lca <- cid.lca[-tmp.rm,]
    tmp <- tmp[-tmp.rm]
  }

  cid.lca.h <- taxid.hierarchy[tmp,]
  for(i in (1:nrow(cid.lca))) {
    cid.lca.h[i, 1:max((cid.lca$lca.level[i]-1),1)] <- NA
  }
  
  cid.lca <- data.table::data.table(
    cid.lca,
    cid.lca.h
  )
  
  
  if(!is.null(out.dir)) {
    save(cid.lca, file = paste0(out.dir, "/cid.lca.Rdata"))
  }
  return(cid.lca)
}





