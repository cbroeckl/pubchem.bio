#' build.pubchem.bio
#'
#' utilizes downloaded and properly formatted local pubchem data created by 'get.pubchem.ftp' function
#' @details utilizes downloaded and properly formatted local pubchem data created by 'get.pubchem.ftp' function
#' @param pc.directory directory from which to load pubchem .Rdata files.  alternatively, provide  R data.tables for ALL cid._property_.object options defined below.  
#' @param use.bio.sources logical.  If TRUE (default) use the bio.source vector of sources, incorporating all CIDs from those bio databases.
#' @param bio.sources vector of source names from which to extract pubchem CIDs.  all can be found here: https://pubchem.ncbi.nlm.nih.gov/sources/.  deafults to c("Metabolomics Workbench", "Human Metabolome Database (HMDB)", "ChEBI", "LIPID MAPS",  "MassBank of North America (MoNA)")
#' @param use.pathways logical.  should all CIDs from any biological pathway data be incorporated into database? 
#' @param pathway.sources character. vector of sources to be used when adding metabolites to pubchem bio database. default = NULL, using all pathway sources.
#' @param use.taxid logical.  should all CIDs associated with a taxonomic identifier (taxid) be used? 
#' @param taxonomy.sources character. vector of sources to be used when adding taxonomically related metabolites to database.  Default = NULL, using all sources.
#' @param remove.salts logical.  should salts be removed from dataset?  default = TRUE.  salts recognized as '.' in smiles string.  performed after 'use.parent.cid'. 
#' @param use.parent.cid logical. should CIDs be replaced with parent CIDs?  default = TRUE.
#' @param get.properties logical. if TRUE, will return rcdk calculated properties:  XLogP, TPSA, HBondDonorCount and HBondAcceptorCount.
#' @param threads integer. how many threads to use when calculating rcdk properties.  parallel processing via DoParallel and foreach packages.  
#' @param rcdk.desc vector. character vector of valid rcdk descriptors.  default = rcdk.desc <- c("org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor", "org.openscience.cdk.qsar.descriptors.molecular.AcidicGroupCountDescriptor", "org.openscience.cdk.qsar.descriptors.molecular.BasicGroupCountDescriptor", "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor"). To see descriptor categories: 'dc <- rcdk::get.desc.categories(); dc' .  To see the descriptors within one category: 'dn <- rcdk::get.desc.names(dc\[4\]); dn'. Note that the four default parameters are relatively fast to calculate - some descriptors take a very long time to calculate.  you can calculate as many as you wish, but processing time will increase the more descriptors are added.   
#' @param cid.sid.object R data.table, generally produced by get.pubchem.ftp; alternatively, define pc.directory
#' @param cid.pwid.object R data.table, generally produced by get.pubchem.ftp; alternatively, define pc.directory
#' @param cid.parent.object R data.table, generally produced by get.pubchem.ftp; alternatively, define pc.directory
#' @param cid.taxid.object R data.table, generally produced by get.pubchem.ftp; alternatively, define pc.directory
#' @param cid.formula.object R data.table, generally produced by get.pubchem.ftp; alternatively, define pc.directory
#' @param cid.smiles.object R data.table, generally produced by get.pubchem.ftp; alternatively, define pc.directory
#' @param cid.inchikey.object R data.table, generally produced by get.pubchem.ftp; alternatively, define pc.directory
#' @param cid.monoisotopic.mass.object R data.table, generally produced by get.pubchem.ftp; alternatively, define pc.directory
#' @param cid.title.object R data.table, generally produced by get.pubchem.ftp; alternatively, define pc.directory
#' @param cid.cas.object R data.table, generally produced by get.pubchem.ftp; alternatively, define pc.directory
#' @param cid.pmid.ct.object R data.table, generally produced by get.pubchem.ftp; alternatively, define pc.directory
#' @return a data frame containing pubchem CID, title, formula, monoisotopic molecular weight, inchikey, smiles, cas, optionally rcdk properties
#' @author Corey Broeckling
#' 
#' @examples
#' data('cid.sid', package = "pubchem.bio")
#' data('cid.pwid', package = "pubchem.bio")
#' data('cid.parent', package = "pubchem.bio")
#' data('cid.taxid', package = "pubchem.bio")
#' data('cid.formula', package = "pubchem.bio")
#' data('cid.smiles, package = "pubchem.bio")
#' data('cid.inchikey', package = "pubchem.bio")
#' data('cid.monoisotopic.mass', package = "pubchem.bio")
#' data('cid.title', package = "pubchem.bio")
#' data('cid.cas', package = "pubchem.bio")
#' data(cid.pmid.ct', package = "pubchem.bio")
#' pc.bio.out <- build.pubchem.bio(use.pathways = FALSE, use.parent.cid = FALSE, get.properties = FALSE, threads = 1,
#' cid.sid.object = cid.sid, cid.pwid.object = cid.pwid, cid.parent.object = cid.parent, cid.taxid.object = cid.taxid,
#' cid.formula.object = cid.formula, cid.smiles.object = cid.smiles, cid.inchikey.object = cid.inchikey, cid.monoisotopic.mass.object = cid.monoisotopic.mass,
#' cid.title.object = cid.title, cid.cas.object = cid.cas, cid.pmid.ct.object = cid.pmid.ct)
#' head(pc.bio.out)
#' 
#' @export 
#' 
#' 
build.pubchem.bio <- function(
    pc.directory =  NULL,
    use.bio.sources = TRUE,
    bio.sources = c(
      "Metabolomics Workbench",
      "Human Metabolome Database (HMDB)",
      "ChEBI",
      "LIPID MAPS",
      "MassBank of North America (MoNA)"
    ),
    use.pathways = TRUE,
    pathway.sources = NULL, 
    use.taxid = TRUE,
    taxonomy.sources = NULL,
    use.parent.cid = TRUE,
    remove.salts = TRUE,
    get.properties = TRUE,
    threads = 8,
    rcdk.desc = c(
      "org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor",
      "org.openscience.cdk.qsar.descriptors.molecular.AcidicGroupCountDescriptor",
      "org.openscience.cdk.qsar.descriptors.molecular.BasicGroupCountDescriptor",
      "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor"
    ),
    cid.sid.object = NULL,
    cid.pwid.object = NULL,
    cid.parent.object = NULL,
    cid.taxid.object = NULL,
    cid.formula.object = NULL,
    cid.smiles.object = NULL,
    cid.inchikey.object = NULL,
    cid.monoisotopic.mass.object = NULL,
    cid.title.object = NULL,
    cid.cas.object = NULL,
    cid.pmid.ct.object = NULL
){  
  
  cid <- vector(length = 0, mode = 'integer')
  
  if(use.parent.cid) {
    if(is.null(cid.parent.object)) {
      load(paste0(pc.directory, "/cid.parent.Rdata"))
      cid.parent <- cid.parent
    } else {
      cid.parent <- cid.parent.object
    }
    
  }
  
  ## get CIDs by data source first
  if(use.bio.sources) {
    if(!is.null(bio.sources)) {
      if(is.null(cid.sid.object)) {
        load(paste0(pc.directory, "/cid.sid.Rdata"))
        cid.sid <- cid.sid
      } else {
        cid.sid <- cid.sid.object
      }

      data.table::setkey(cid.sid, "cid")
      keep <- cid.sid$source %in% bio.sources
      source.cid <- cid.sid$cid[keep]
      if(any(is.na(source.cid))) {
        source.cid <- source.cid[-which(is.na(source.cid))]
      }
      if(use.parent.cid) {
        data.table::setkey(cid.parent, "cid")
        m <- match(source.cid, cid.parent$cid)
        parent.cid <- cid.parent$parent.cid[m]
        parent.cid[which(is.na(parent.cid))] <- source.cid[which(is.na(parent.cid))]
        source.cid <- parent.cid
        # cid <- sort(unique(cid))
        # message(" - after replacing CID with parent CID, current unique cid count:" , length(cid), '\n')
        rm(parent.cid); rm(m); gc()
      }
      source.cid.table <- table(source.cid)
      cid <- sort(unique(c(cid, source.cid)))
      message(" - added", length(source.cid.table), "cids based on bio.sources", '\n')
      rm(cid.sid); rm(source.cid); rm(keep); gc()
    }  ## which(cid == 187)
  }
  
  if(use.pathways) {
    if(is.null(cid.pwid.object)) {
      load(paste0(pc.directory, "/cid.pwid.Rdata"))
      cid.pwid <- cid.pwid
    } else {
      cid.pwid <- cid.pwid.object
    }

    data.table::setkey(cid.pwid, "cid")
    if(!is.null(pathway.sources)) {
      keep <- cid.pwid$source %in% pathway.sources
    } else {
      keep <- 1:nrow(cid.pwid)
    }
    path.cid <- cid.pwid$cid[keep]
    if(any(is.na(path.cid))) {
      path.cid <- path.cid[-which(is.na(path.cid))]
    }
    if(use.parent.cid) {
      m <- match(path.cid, cid.parent$cid)
      parent.cid <- cid.parent$parent.cid[m]
      parent.cid[which(is.na(parent.cid))] <- path.cid[which(is.na(parent.cid))]
      path.cid <- parent.cid
      # cid <- sort(unique(cid))
      # message(" - after replacing CID with parent CID, current unique cid count:" , length(cid), '\n')
      rm(parent.cid); rm(m); gc()
    }
    path.cid.table <- table(path.cid)
    cid <- sort(unique(c(cid, path.cid)))
    message(" - added", length(path.cid.table), "cids based on pathways. current unique cid count:" , length(cid), '\n')
    rm(path.cid); rm(cid.pwid); rm(keep); gc()
  }
  
  if(use.taxid) {
    if(is.null(cid.taxid.object)) {
      load(paste0(pc.directory, "/cid.taxid.Rdata"))
      cid.taxid <- cid.taxid
    } else {
      cid.taxid <- cid.taxid.object
    }

    data.table::setkey(cid.taxid, "cid")
    if(!is.null(taxonomy.sources)) {
      keep <- cid.taxid$source %in% taxonomy.sources
    } else {
      keep <- 1:nrow(cid.taxid)
    }
    tax.cid <- cid.taxid$cid[keep]
    if(any(is.na(tax.cid))) {
      tax.cid <- tax.cid[-which(is.na(tax.cid))]
    }
    if(use.parent.cid) {
      data.table::setkey(cid.parent, "cid")
      m <- match(tax.cid, cid.parent$cid)
      parent.cid <- cid.parent$parent.cid[m]
      parent.cid[which(is.na(parent.cid))] <- tax.cid[which(is.na(parent.cid))]
      tax.cid <- parent.cid
      # cid <- sort(unique(cid))
      # message(" - after replacing CID with parent CID, current unique cid count:" , length(cid), '\n')
      rm(parent.cid); rm(m); gc()
    }
    tax.cid.table <- table(tax.cid)
    cid <- sort(unique(c(cid, tax.cid)))
    message(" - added", length(tax.cid.table), "cids based on taxonomy. current unique cid count:" , length(cid), '\n')
    rm(tax.cid); rm(cid.taxid); gc()
  }
  
  
  if(use.parent.cid) {
    rm(cid.parent); gc()
  }
  # create output data.frame
  # CID, name(title), formula, monoisotopic molecular weight, inchikey, smiles, cas, optionally pubchem properties
  
  message(" - extracting descriptors from files:" , '\n')
  
  ## formula
  message(" ---- formula" ,  format(Sys.time()), '\n')
  if(is.null(cid.formula.object)) {
    load(paste0(pc.directory, "/cid.formula.Rdata"))
    cid.formula <- cid.formula
  } else {
    cid.formula <- cid.formula.object
  }

  data.table::setkey(cid.formula, "cid")
  m <- match(cid, cid.formula$cid)
  formula <- cid.formula$formula[m]
  m <- match(cid, cid.formula$cid)
  parent.formula <- cid.formula$formula[m]
  rm(cid.formula); rm(m); gc()
  
  ## smiles
  message(" ---- smiles" ,  format(Sys.time()), '\n')
  if(is.null(cid.smiles)) {
    load(paste0(pc.directory, "/cid.smiles.Rdata"))
    cid.smiles <- cid.smiles
  } else {
    cid.smiles <- cid.smiles.object
  }

  data.table::setkey(cid.smiles, "cid")
  m <- match(cid, cid.smiles$cid)
  smiles <- cid.smiles$smiles[m]
  m <- match(cid, cid.smiles$cid)
  parent.smiles <- cid.smiles$smiles[m]
  rm(cid.smiles); rm(m); gc()
  
  ## choose whether to use original or parent CID.
  ## if neither is a salt, use parent
  ## if both are salts, use parent
  ## if only one is a salt, use the one that isn't
  # if(use.parent.cid) {
  #   bio.cid <- parent.cid
  #   # parent.charged <- grepl("-", parent.formula, fixed = TRUE) | grepl("+", parent.formula, fixed = TRUE)
  #   parent.salt <- grepl(".", parent.smiles, fixed = TRUE)
  #   # bio.charged <- grepl("-", formula, fixed = TRUE) | grepl("+", formula, fixed = TRUE)
  #   bio.salt <- grepl("-", smiles, fixed = TRUE)
  #   use.parent <- parent.salt == bio.salt
  #   disagreement <- which(!use.parent)    
  #   use.parent[disagreement] <- bio.salt[disagreement]
  #   cid[use.parent] <- parent.cid[use.parent]
  #   formula[use.parent] <- parent.formula[use.parent]
  #   smiles[use.parent] <- parent.smiles[use.parent]
  # }
  
  ## monoisotopic mass
  message(" ---- monoisotopic mass" ,  format(Sys.time()), '\n')
  if(is.null(cid.monoisotopic.mass.object)) {
    load(paste0(pc.directory, "/cid.monoisotopic.mass.Rdata"))
    cid.monoisotopic.mass <- cid.monoisotopic.mass
  } else {
    cid.monoisotopic.mass <- cid.monoisotopic.mass.object
  }

  data.table::setkey(cid.monoisotopic.mass, "cid")
  m <- match(cid, cid.monoisotopic.mass$cid)
  monoisotopic.mass <- cid.monoisotopic.mass$monoisotopic.mass[m]
  rm(cid.monoisotopic.mass); rm(m); gc()
  
  ## inchikey
  message(" ---- inchikey" ,  format(Sys.time()), '\n')
  if(is.null(cid.inchikey.object)) {
    load(paste0(pc.directory, "/cid.inchikey.Rdata"))
    cid.inchikey <- cid.inchikey
  } else {
    cid.inchikey <- cid.inchikey.object
  }

  data.table::setkey(cid.inchikey, "cid")
  m <- match(cid, cid.inchikey$cid)
  inchikey <- cid.inchikey$inchikey[m]
  rm(cid.inchikey); rm(m); gc()
  
  ## title
  message(" ---- compound name" ,  format(Sys.time()), '\n')
  if(is.null(cid.title.object)) {
    load(paste0(pc.directory, "/cid.title.Rdata"))
    cid.title <- cid.title
  } else {
    cid.title <- cid.title.object
  }

  data.table::setkey(cid.title, "cid")
  m <- match(cid, cid.title$cid)
  name <- cid.title$title[m]
  rm(cid.title); rm(m); gc()
  
  ## cas
  message(" ---- cas number" ,  format(Sys.time()), '\n')
  if(is.null(cid.cas.object)) {
    load(paste0(pc.directory, "/cid.cas.Rdata"))
    cid.cas <- cid.cas
  } else {
    cid.cas <- cid.cas.object
  }

  data.table::setkey(cid.cas, "cid")
  m <- match(cid, cid.cas$cid)
  cas <- cid.cas$cas[m]
  rm(cid.cas); rm(m); gc()
  
  ## pmid count
  message(" ---- pubmed count" ,  format(Sys.time()), '\n')
  if(is.null(cid.pmid.ct.object)) {
    load(paste0(pc.directory, "/cid.pmid.ct.Rdata"))
    cid.pmid.ct <- cid.pmid.ct
  } else {
    cid.pmid.ct <- cid.pmid.ct.object
  }

  data.table::setkey(cid.pmid.ct, "cid")
  m <- match(cid, cid.pmid.ct$cid)
  pmid.ct <- cid.pmid.ct$pmid.ct[m]
  if(any(is.na(pmid.ct))) {
    pmid.ct[is.na(pmid.ct)] <- 0
  }
  rm(cid.pmid.ct); rm(m); gc()
  
  ## first block of inchikey - same bonding
  message(" ---- first block of inchikey" , '\n')
  inchikey.first.block <- sapply(1:length(inchikey), FUN = function(x){unlist(strsplit(inchikey[x], "-"))[1]})
  
  out <- data.table::data.table(
    cid, 
    name, 
    formula,
    monoisotopic.mass,
    inchikey,
    inchikey.first.block,
    smiles,
    cas,
    pmid.ct
  )
  
  if(any(ls() == 'source.cid.table')) {
    source.cid.table <- data.table::data.table(
      'cid' = as.integer(as.numeric(names(source.cid.table))),
      'count' = as.vector(source.cid.table)
    )
    data.table::setkey(source.cid.table, "cid")
    m <- match(cid, source.cid.table$cid)
    source.ct <- as.integer(source.cid.table$count[m])
    source.ct[which(is.na(source.ct))] <- 0L
    out <- data.table::data.table(
      out,
      source.ct
    )
    rm(source.cid.table); rm(source.ct)
  }
  

  if(any(ls() == 'path.cid.table')) {
    path.cid.table <- data.table::data.table(
      'cid' = as.integer(as.numeric(names(path.cid.table))),
      'count' = as.vector(path.cid.table)
    )
    data.table::setkey(path.cid.table, "cid")
    m <- match(cid, path.cid.table$cid)
    pathway.ct <- as.integer(path.cid.table$count[m])
    pathway.ct[which(is.na(pathway.ct))] <- 0L
    out <- data.table::data.table(
      out,
      pathway.ct
    )
    rm(path.cid.table); rm(pathway.ct)
  }
  
  if(any(ls() == 'tax.cid.table')) {
    tax.cid.table <- data.table::data.table(
      'cid' = as.integer(as.numeric(names(tax.cid.table))),
      'count' = as.vector(tax.cid.table)
    )
    data.table::setkey(tax.cid.table, "cid")
    m <- match(cid, tax.cid.table$cid)
    taxonomy.ct <- as.integer(tax.cid.table$count[m])
    taxonomy.ct[which(is.na(taxonomy.ct))] <- 0L
    out <- data.table::data.table(
      out,
      taxonomy.ct
    )
    rm(tax.cid.table); rm(taxonomy.ct)
  }
  
  rm.rows <- unique(which(is.na(out$smiles)))
  if(length(rm.rows) > 0) {
    out <- out[-rm.rows,]
    message(" - removed", length(rm.rows), "with missing structures,", "current unique cid count:" , nrow(out), '\n')
  }
  
  ## remove duplicated rows
  out <- unique(out[,1:ncol(out)])
  
  ## remove salts
  if(remove.salts) {
    rm.rows <- grep(".", out$smiles, fixed = TRUE)
    if(length(rm.rows) > 0) {
      out <- out[-rm.rows,]
    }
  }
  
  ## create stable dopar function:
  `%dopar%` <- foreach::`%dopar%`
  
  ## get simple physical-chemical properties
  if(get.properties) {
    message(" - calclulating rcdk properties",  format(Sys.time()), '\n')
    cid.list <- as.list(out$cid)
    sm.list <- as.list(out$smiles)
    doParallel::registerDoParallel(cl <- parallel::makeCluster(threads))
    results <- foreach::foreach(i = 1:(length(cid.list))) %dopar% {
      i <- i
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
    doParallel::stopImplicitCluster()
  }
  message(" - rcdk properties completed", format(Sys.time()), '\n')
  pc.bio <- out
  rm(out)
  gc()
  
  
  save(pc.bio, file = paste0(pc.directory, "/pc.bio.Rdata"))
  return(pc.bio)
}
## pc.bio <- build.pubchem.bio(pc.directory = "R:/RSTOR-PMF/Software/db/met.db/20241216")


