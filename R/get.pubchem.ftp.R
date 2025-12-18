#' get.pubchem.ftp
#'
#' first step to building a local selective, biologically focused, pubchem data repository focused on metabolomics informatics 
#' @details this function downloads and unzips files from pubchem and NCBI taxonomy FTP sites as a first step in building a local metabolomics repository.
#' @param pc.directory character. directory to which data will be saved
#' @param timeout numeric.  timeout setting for FTP download.  setting options(timeout) value too small will generate errors for large files. default = 50000.
#' @param rm.tmp.files logical.  should temporary files be removed after completion of download and parsing?  Default = TRUE. 
#' @param threads integer. the number of parallel threads to be used by foreach %dopar% during processing of taxonomy hierarchy data.  
#' @return nothing.  all data are saved to disk for later loading
#' @author Corey Broeckling
#' 
#' @examples
#' \dontrun{
#' my.dir <- "C:/Temp/20250725"
#' # or some other valid directory.
#' # this will be created assuming 'C:/Temp' exists.
#' get.pubchem.ftp(
#'     pc.directory = my.dir,
#'     timeout = 50000,
#'     rm.tmp.files = TRUE
#' )
#' }
#' 
#' @export 
#' 

get.pubchem.ftp <- function(
    pc.directory = NULL,
    timeout = 50000,
    rm.tmp.files = TRUE,
    threads = 2
) {
  
  ## setup directory structure
  if(is.null(pc.directory)) stop("you must define the directory that you wish to save PubChem data to.", '\n')
  pc.directory <- suppressWarnings(normalizePath(pc.directory))
  if(!dir.exists(pc.directory)) dir.create(pc.directory)
  tmp.dir <- suppressWarnings(paste0(pc.directory, "/tmp/"))
  dir.create(tmp.dir)
  
  loadNamespace('data.table')
  .datatable.aware = TRUE
  ..pc.bio.keep.cols <- NULL
  
  message(" -- writing data to ", pc.directory, '\n')
  
  readme <- c(
    "Data collection derived from pubchem, primarily using the 'compound' FTP data.", '\n',
    "generation date:", date(), '\n')
  
  ## set download timeout to prevent errors
  old <- base::options()
  on.exit(base::options(old))
  options(timeout=timeout)
  
  
  ## data are imported, generally, as a data.table using fread
  ## these data can also be indexed by cid using the data.table::setkey function
  ## i.e. data.table::setkey(d, "cid")
  ## making "d$cid %in% my.cid.vector" really fast for later searching.
  ## store dataset as a data.table, indexed by 'cid'. 
  ## if when loading one needs to search by anything other than cid,
  ## best to reindex by search variable first, to dramatically reduce search time
  
  ########################
  ## CID to Preferred CID
  ########################
  
  if(file.exists(paste0(pc.directory, "/cid.preferred.Rdata"))) {load(paste0(pc.directory, "/cid.preferred.Rdata"))} else {
    download.url <- "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Preferred.gz"
    tmp.file <- paste0(tmp.dir, basename(download.url))
    utils::download.file(url = download.url, destfile = tmp.file)
    R.utils::gunzip(tmp.file, remove = FALSE, destname = paste0(tmp.dir, gsub(".gz", "", basename(tmp.file))))
    # untar(paste0(gsub(".gz", "", tmp.file)), exdir = tmp.dir)
    
    #### cid.to.preferred.cid.relationship
    d <- data.table::fread(paste0(tmp.dir, basename(tmp.file)))
    names(d) <- c("cid", "preferred.cid")
    data.table::setkey(d, "cid")
    cid.preferred <- d
    rm(d)
    data.table::setkey(cid.preferred, "cid")
    save(cid.preferred, file = paste0(pc.directory, "/cid.preferred.Rdata"))
    ### this is the only dataset we leave open, all others are immediately removed after saving
    ### replace CID with preferred CID for all datasets
    
    readme <- c(
      readme,
      " - cid.preferred.Rdata, is a data.table of two columns, headers 'cid' and 'preferred.cid', each integer values. data.table is indexed by 'cid'. ", '\n',
      "   NOTE THAT: for all values in all downloaded datasets, FTP available CID values for all associations have been replaced with their preferred CIDs, based on this table. ", '\n', '\n'
    )
  }
  
  ########################
  ## CID to parent CID
  ########################
  if(file.exists(paste0(pc.directory, "/cid.parent.Rdata"))) {
    
  } else {
    download.url <- "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Parent.gz"
    tmp.file <- paste0(tmp.dir, basename(download.url))
    utils::download.file(url = download.url, destfile = tmp.file)
    R.utils::gunzip(tmp.file, remove = FALSE, destname = paste0(tmp.dir, gsub(".gz", "", basename(tmp.file))))
    # untar(paste0(gsub(".gz", "", tmp.file)), exdir = tmp.dir)
    
    #### cid.to.parent.relationship
    d <- data.table::fread(paste0(tmp.dir, basename(tmp.file)))
    names(d) <- c("cid", "parent.cid")
    data.table::setkey(d, "cid")
    tmp <- match(d$cid, cid.preferred$cid)
    use <- which(!is.na(tmp))
    if(length(use)>0) {
      d$cid[use] <- cid.preferred$preferred.cid[tmp[use]]
    }
    cid.parent <- d
    rm(d)
    data.table::setkey(cid.parent, "cid")
    save(cid.parent, file = paste0(pc.directory, "/cid.parent.Rdata"))
    rm(cid.parent)
    gc()
    
    readme <- c(
      readme,
      " - cid.parent.Rdata, is a data.table of two columns, headers 'cid' and 'parent.cid', each integer values. data.table is indexed by 'cid'. ", '\n', '\n'
    )
  }
  ########################
  ## CID to InChIKey
  ########################
  if(!file.exists(paste0(pc.directory, "/cid.inchi.Rdata"))) {
    download.url <- "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-InChI-Key.gz"
    tmp.file <- paste0(tmp.dir, basename(download.url))
    utils::download.file(url = download.url, destfile = tmp.file)
    R.utils::gunzip(tmp.file, remove = FALSE, destname = paste0(tmp.dir, gsub(".gz", "", basename(tmp.file))))
    # untar(paste0(gsub(".gz", "", tmp.file)), exdir = tmp.dir)
    
    #### cid.to.inchikey.relationship
    d <- data.table::fread(paste0(tmp.dir, basename(tmp.file)))
    names(d) <- c("cid", "inchi", "inchikey")
    data.table::setkey(d, "cid")
    tmp <- match(d$cid, cid.preferred$cid)
    use <- which(!is.na(tmp))
    if(length(use)>0) {
      d$cid[use] <- cid.preferred$preferred.cid[tmp[use]]
    }
    pc.bio.keep.cols <- c("cid", "inchikey")
    cid.inchikey <- d[, ..pc.bio.keep.cols]
    data.table::setkey(cid.inchikey, "cid")
    save(cid.inchikey, file = paste0(pc.directory, "/cid.inchikey.Rdata"))
    rm(cid.inchikey);   gc()
    pc.bio.keep.cols <- c("cid", "inchi")
    cid.inchi <- d[, ..pc.bio.keep.cols]
    data.table::setkey(cid.inchi, "cid")
    save(cid.inchi, file = paste0(pc.directory, "/cid.inchi.Rdata"))
    rm(cid.inchi); gc()
    rm(d); gc()
  }
  readme <- c(
    readme,
    " - cid.inchikey.Rdata, is a data.table of two columns, headers 'cid'and 'inchikey', with integer and character values, respectively. data.table is indexed by 'cid'. ", '\n', '\n',
    " - cid.inchi.Rdata, is a data.table of two columns, headers 'cid'and 'inchi', with integer and character values, respectively. data.table is indexed by 'cid'. ", '\n', '\n'
  )
  
  
  ########################
  ## CID to SMILES
  ########################
  if(!file.exists(paste0(pc.directory, "/cid.smiles.Rdata"))) {
    download.url <- "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz"
    tmp.file <- paste0(tmp.dir, basename(download.url))
    utils::download.file(url = download.url, destfile = tmp.file)
    R.utils::gunzip(tmp.file, remove = FALSE, destname = paste0(tmp.dir, gsub(".gz", "", basename(tmp.file))))
    # untar(paste0(gsub(".gz", "", tmp.file)), exdir = tmp.dir)
    
    #### cid.to.smiles.relationship
    d <- data.table::fread(paste0(tmp.dir, basename(tmp.file)))
    names(d) <- c("cid", "smiles")
    data.table::setkey(d, "cid")
    tmp <- match(d$cid, cid.preferred$cid)
    use <- which(!is.na(tmp))
    if(length(use)>0) {
      d$cid[use] <- cid.preferred$preferred.cid[tmp[use]]
    }
    cid.smiles <- d
    rm(d); gc()
    data.table::setkey(cid.smiles, "cid")
    save(cid.smiles, file = paste0(pc.directory, "/cid.smiles.Rdata"))
    rm(cid.smiles); gc()
  }
  readme <- c(
    readme,
    " - cid.smiles.Rdata, is a data.table of two columns, headers 'cid'and 'smiles', with integer and character values, respectively. data.table is indexed by 'cid'. ", '\n', '\n'
  )
  
  
  ########################
  ## CID to Title
  ########################
  if(!file.exists(paste0(pc.directory, "/cid.title.Rdata"))) {
    download.url <- "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Title.gz"
    tmp.file <- paste0(tmp.dir, basename(download.url))
    utils::download.file(url = download.url, destfile = tmp.file)
    R.utils::gunzip(tmp.file, remove = FALSE, destname = paste0(tmp.dir, gsub(".gz", "", basename(tmp.file))))
    # untar(paste0(gsub(".gz", "", tmp.file)), exdir = tmp.dir)
    
    #### cid.to.title.relationship
    d <- data.table::fread(paste0(tmp.dir, basename(tmp.file)), quote="")
    names(d) <- c("cid", "title")
    data.table::setkey(d, "cid")
    tmp <- match(d$cid, cid.preferred$cid)
    use <- which(!is.na(tmp))
    if(length(use)>0) {
      d$cid[use] <- cid.preferred$preferred.cid[tmp[use]]
    }
    cid.title <- d
    rm(d); gc()
    data.table::setkey(cid.title, "cid")
    save(cid.title, file = paste0(pc.directory, "/cid.title.Rdata"))
    rm(cid.title); gc
    
  }
  readme <- c(
    readme,
    " - cid.title.Rdata, is a data.table of two columns, headers 'cid'and 'title', with integer and character values, respectively. data.table is indexed by 'cid'. ", '\n', '\n'
  )
  
  ########################
  ## CID to PMID
  ########################
  if(!file.exists(paste0(pc.directory, "/cid.pmid.ct.Rdata"))) {
    
    download.url <- "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-PMID.gz"
    tmp.file <- paste0(tmp.dir, basename(download.url))
    utils::download.file(url = download.url, destfile = tmp.file)
    R.utils::gunzip(tmp.file, remove = FALSE, destname = paste0(tmp.dir, gsub(".gz", "", basename(tmp.file))))
    # untar(paste0(gsub(".gz", "", tmp.file)), exdir = tmp.dir)
    
    #### cid.to.pmid.relationship
    d <- data.table::fread(paste0(tmp.dir, basename(tmp.file)))
    if(ncol(d) > 2) {
      names(d) <- c(c("cid", "pmid"), paste0("unk", 1:(ncol(d)-2)))
    } else {
      names(d) <- c("cid", "pmid")
    }
    pc.bio.keep.cols <- c("cid", "pmid")
    d <- d[, ..pc.bio.keep.cols]
    
    data.table::setkey(d, "cid")
    tmp <- match(d$cid, cid.preferred$cid)
    use <- which(!is.na(tmp))
    if(length(use)>0) {
      d$cid[use] <- cid.preferred$preferred.cid[tmp[use]]
    }
    cid.pmid <- d
    rm(d); gc()
    data.table::setkey(cid.pmid, "cid")
    save(cid.pmid, file = paste0(pc.directory, "/cid.pmid.Rdata"))
    
    cid.pmid.ct <- table(cid.pmid$cid)
    cid.pmid.ct <- data.table::data.table(cid.pmid.ct)
    names(cid.pmid.ct) <- c('cid', 'pmid.ct')
    cid.pmid.ct$cid <- as.integer(as.numeric(cid.pmid.ct$cid))
    data.table::setkey(cid.pmid.ct, "cid")
    save(cid.pmid.ct, file = paste0(pc.directory, "/cid.pmid.ct.Rdata"))
    
    rm(cid.pmid); gc()
    
  }
  readme <- c(
    readme,
    " - cid.pmid.Rdata, is a data.table of two columns, headers 'cid'and 'pmid', with integer and character values, respectively. data.table is indexed by 'cid'. ", '\n', '\n',
    " - cid.pmid.ct.Rdata, is a data.table of two columns, headers 'cid'and 'pmid.ct', each with integer values.  'pmid.ct' represents the number of pubmed references for that compound. data.table is indexed by 'cid'. ", '\n', '\n'
  )
  
  ########################
  ## CID to Mass
  ########################
  if(!file.exists(paste0(pc.directory, "/cid.accurate.mass.Rdata"))) {
    download.url <- "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Mass.gz"
    tmp.file <- paste0(tmp.dir, basename(download.url))
    if(!file.exists(paste0(tmp.dir, gsub(".gz", "", basename(tmp.file))))) {
      utils::download.file(url = download.url, destfile = tmp.file)
      R.utils::gunzip(tmp.file, remove = FALSE, destname = paste0(tmp.dir, gsub(".gz", "", basename(tmp.file))))
    }

    # untar(paste0(gsub(".gz", "", tmp.file)), exdir = tmp.dir)
    
    #### cid.to.mass.relationship
    d <- data.table::fread(paste0(tmp.dir, basename(tmp.file)))
    names(d) <- c("cid", "formula", "monoisotopic.mass", "accurate.mass")
    # data.table::setkey(d, "cid")
    # d <- as.data.frame(d)
    is.750 <- which(d$cid == 750)
    check.750 <- round(d$monoisotopic.mass[is.750], 2) == 75.03
    if(!check.750) {
      stop("accurate mass of cid 750 in d is reported as", round(d$monoisotopic.mass[is.750], 2), "and should be 75.03" )
    } else {
      warning("accurate mass of cid 750 in d is as expected")
    }
    tmp <- match(d$cid, cid.preferred$cid)
    use <- which(!is.na(tmp))
    if(length(use)>0) {
      d$cid[use] <- cid.preferred$preferred.cid[tmp[use]]
    }
    pc.bio.keep.cols <- c("cid", "formula")
    cid.formula <- d[, ..pc.bio.keep.cols]
    data.table::setkey(cid.formula, "cid")
    save(cid.formula, file = paste0(pc.directory, "/cid.formula.Rdata"))
    rm(cid.formula); gc()
    pc.bio.keep.cols <- c("cid", "monoisotopic.mass")
    cid.monoisotopic.mass <- d[, ..pc.bio.keep.cols]
    is.750 <- which(cid.monoisotopic.mass$cid == 750)
    check.750 <- round(cid.monoisotopic.mass$monoisotopic.mass[is.750], 2) == 75.03
    if(!check.750) {
      stop("accurate mass of cid 750 reported as", round(cid.monoisotopic.mass$monoisotopic.mass[is.750], 2), "and should be 75.03" )
    } else {
      warning("accurate mass of cid 750 is as expected")
      
    }
    data.table::setkey(cid.monoisotopic.mass, "cid")
    save(cid.monoisotopic.mass, file = paste0(pc.directory, "/cid.monoisotopic.mass.Rdata"))
    rm(cid.monoisotopic.mass); gc()
    pc.bio.keep.cols <- c("cid", "accurate.mass")
    cid.accurate.mass <- d[, ..pc.bio.keep.cols]
    data.table::setkey(cid.accurate.mass, "cid")
    save(cid.accurate.mass, file = paste0(pc.directory, "/cid.accurate.mass.Rdata"))
    rm(cid.accurate.mass); gc()
    rm(d); gc()
  }
  readme <- c(
    readme,
    " - cid.formula.Rdata, is a data.table of two columns, headers 'cid'and 'formula', with integer and character values, respectively. data.table is indexed by 'cid'. ", '\n', '\n',
    " - cid.monoisotopic.mass.Rdata, is a data.table of two columns, headers 'cid'and 'monoisotopic.mass', with integer and character values, respectively. data.table is indexed by 'cid'. ", '\n', '\n',
    " - cid.accurate.mass.Rdata, is a data.table of two columns, headers 'cid'and 'accurate.mass', with integer and character values, respectively. data.table is indexed by 'cid'. ", '\n', '\n'
  )
  
  ########################
  ## CID to MeSH name
  ########################
  if(!file.exists(paste0(pc.directory, "/cid.mesh.name.Rdata"))) {
    download.url <- "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-MeSH"
    tmp.file <- paste0(tmp.dir, basename(download.url))
    utils::download.file(url = download.url, destfile = tmp.file)
    # R.utils::gunzip(tmp.file, remove = FALSE)
    # untar(paste0(gsub(".gz", "", tmp.file)), exdir = tmp.dir)
    
    #### cid.to.mesh.name relationship
    d <- readLines(paste0(tmp.dir, basename(tmp.file)), encoding = "UTF-8")
    d <- lapply(1:length(d), FUN = function(x) {
      # cat(d[x], '\n')
      tmp <- unlist(strsplit(d[x], '\t', fixed = TRUE, useBytes = TRUE))
      tmp <- data.frame(
        "cid" = rep(tmp[1], length(tmp)-1),
        "mesh.name" = tmp[2:length(tmp)]
      )
      tmp
    })
    d <- as.data.frame(do.call(rbind, d))
    d[,1] <- as.numeric(d[,1])
    d <- data.table::as.data.table(d)
    data.table::setkey(d, "cid")
    tmp <- match(d$cid, cid.preferred$cid)
    use <- which(!is.na(tmp))
    if(length(use)>0) {
      d$cid[use] <- cid.preferred$preferred.cid[tmp[use]]
    }
    cid.mesh.name <- d
    rm(d)
    data.table::setkey(cid.mesh.name, "cid")
    save(cid.mesh.name, file = paste0(pc.directory, "/cid.mesh.name.Rdata"))
    
  }
  
  readme <- c(
    readme,
    " - cid.mesh.name.Rdata, is a data.table of two columns, headers 'cid'and 'mesh.name', with integer and character values, respectively. data.table is indexed by 'cid'. ", '\n', '\n'
  )
  
  
  #### must use separate file for mesh name to mesh function
  if(!file.exists(paste0(pc.directory, "/cid.mesh.function.Rdata"))) {
    download.url <- "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/MeSH-Pharm"
    tmp.file <- paste0(tmp.dir, basename(download.url))
    utils::download.file(url = download.url, destfile = tmp.file)
    # R.utils::gunzip(tmp.file, remove = FALSE)
    # untar(paste0(gsub(".gz", "", tmp.file)), exdir = tmp.dir)
    
    #### cid.to.mesh.relationship
    d <- readLines(paste0(tmp.dir, basename(tmp.file)))
    d <- lapply(1:length(d), FUN = function(x) {
      tmp <- unlist(strsplit(d[x], '\t', fixed = TRUE))
      tmp <- data.frame(
        "mesh.name" = rep(tmp[1], length(tmp)-1),
        "mesh.function" = tmp[2:length(tmp)]
      )
      tmp
    })
    d <- as.data.frame(do.call(rbind, d))
    d <- data.table::as.data.table(d)
    data.table::setkey(d, "mesh.name")
    tmp <- match(d$cid, cid.preferred$cid)
    use <- which(!is.na(tmp))
    if(length(use)>0) {
      d$cid[use] <- cid.preferred$preferred.cid[tmp[use]]
    }
    mesh.name.mesh.function <- d
    rm(d)
    
    
    ## merge two files into one
    d <- merge(cid.mesh.name, mesh.name.mesh.function, by = "mesh.name", all.y = TRUE, all.x = FALSE)
    tmp <- match(d$cid, cid.preferred$cid)
    use <- which(!is.na(tmp))
    if(length(use)>0) {
      d$cid[use] <- cid.preferred$preferred.cid[tmp[use]]
    }
    cid.mesh.function <- d
    cid.mesh.function <- data.table::as.data.table((cid.mesh.function))
    cid.mesh.function <- cid.mesh.function[,c(2:3)]
    data.table::setkey(d, "cid")
    rm(d)
    data.table::setkey(cid.mesh.function, "cid")
    save(cid.mesh.function, file = paste0(pc.directory, "/cid.mesh.function.Rdata"))
    rm(mesh.name.mesh.function); gc()
    rm(cid.mesh.name); gc()
    
  }
  
  readme <- c(
    readme,
    " - cid.mesh.function.Rdata, is a data.table of two columns, headers 'cid'and 'mesh.function', with integer and character values, respectively. data.table is indexed by 'cid'. ", '\n', '\n'
  )
  
  
  ########################
  ## TaxID to taxonomy
  ########################
  download.url <- "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
  tmp.file <- paste0(tmp.dir, download.url)
  if(!file.exists(paste0(tmp.dir, basename(tmp.file)))) {
    utils::download.file(url = download.url, destfile = paste0(tmp.dir, basename(tmp.file)))
    R.utils::gunzip(paste0(tmp.dir, basename(tmp.file)), remove = FALSE, 
                    destname = gsub(".gz", "", paste0(tmp.dir, basename(tmp.file))))
    utils::untar(gsub(".gz", "", paste0(tmp.dir, basename(tmp.file))), exdir = tmp.dir)
  }
  
  #### taxid.to.rank.relationship
  ## get taxid taxname first
  suppressWarnings(d <- readLines(
    paste0(tmp.dir, "names.dmp"))
  )
  d <- lapply(1:length(d), FUN = function(x) {unlist(strsplit(d[x], "\t|\t", fixed = TRUE))})
  d <- as.data.frame(do.call(rbind, d))
  d <- d[grepl("scientific name", d$V4),]
  d <- d[,c(1, 2),]
  names(d) <- c("taxid", "taxname")
  d <- data.table::as.data.table(d)
  taxid.taxname <- d
  row.names(taxid.taxname) <- NULL
  rm(d); gc()
  
  suppressWarnings(d <- readLines(
    paste0(tmp.dir, "nodes.dmp"))
  )
  d <- lapply(1:length(d), FUN = function(x) {unlist(strsplit(d[x], "\t|\t", fixed = TRUE))})
  d <- as.data.frame(do.call(rbind, d))
  names(d) <- c("taxid", "parent.taxid", "rank", "embl.code", "division.id", 
                "inherited.division.flag", "genetic.code.id",
                "inherited.gc.flag", "mitocondrial.genetic.code.id", 
                "inherited.mgc.flag", "genbank.hidden.flag", "hidden.subtree.root.flag", 
                "comments")
  d <- d[, c("taxid", "parent.taxid", "rank", "division.id"),]
  gc()
  
  d <- data.table::as.data.table(d)
  data.table::setkey(d, "taxid")
  d <- merge(d, taxid.taxname, by = 'taxid', all.x = TRUE, all.y = FALSE)
  d <- data.table::as.data.table(d)
  data.table::setkey(d, "taxid")
  taxid.taxname <- d
  save(taxid.taxname, file = paste0(pc.directory, "/taxid.taxname"))
  rm(taxid.taxname); gc()
  
  ## determine the rank of the parent taxid
  d$rank.parent.taxid <- d$rank[match(d$parent.taxid, d$taxid)]
  
  ## define '1' as root
  d[which(d$taxid == "1" & d$parent.taxid == "1"), "rank.parent.taxid"] <- "root"
  d[which(d$taxid == "1" & d$parent.taxid == "1"), "rank"] <- "root"
  
  
  ## start from child as species
  sp.id <- d$taxid[which(d$rank == "species")]
  names(sp.id) <- rep("species", length(sp.id))
  
  ## parent lineage
  par.id <- list()
  par.id[[1]] <- d$taxid[which(d$rank == "species")]
  names(par.id[[1]]) <- rep("species", length(par.id[[1]]))
  
  keep.going <- TRUE
  i <- 1
  unique.levs <- length(sp.id) + 1
  while(keep.going) {
    mtch <- match(par.id[[i]], d$taxid)
    par.id[[i+1]] <- d$parent.taxid[mtch]
    names(par.id[[i+1]]) <- d$rank.parent.taxid[mtch]
    old.unique.levs <- unique.levs
    unique.levs <- length(unique(par.id[[i+1]]))
    i <- i + 1
    if(old.unique.levs == unique.levs) {keep.going = FALSE}
  }
  
  # table(par.id[[36]])
  # threads = 8
  cl <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  base::on.exit({
    parallel::stopCluster(cl)
    rm(cl)
  }) 
  
  `%dopar%` <- foreach::`%dopar%`
  j <- NULL
  x <- NULL
  i <- NULL
  k <- NULL
  
  ### turn par.id into two data frames, 
  ### one contains values, the second names
  par.id.names <- par.id
  for(i in 1:length(par.id)){
    par.id.names[[i]] <- names(par.id[[i]])
  }
  par.id.names <- data.frame(par.id.names)
  par.id.df <- data.frame(par.id)
  ranks <- c( "sub.subspecies", "subspecies",
              "species", "genus", "subfamily", "family",  
              "suborder", "order", "subclass", "class", 
              "subphylum", "phylum", "kingdom", "domain")
  
  ## create species lineage for each species, from par.id
  ## output is a list if length n.species, which will be coerced into a data.frame
  
  sp.lineage <- foreach::foreach(j = 1:length(sp.id)) %dopar% {
    lineage <- c(sp.id[j], sapply(1:length(par.id), FUN = function(x) par.id[[x]][j]))
    lineage <- lineage[1:which(lineage == "1")[1]]
    lineage <- lineage[names(lineage) %in% ranks]
    lineage <- lineage[ranks]
    lineage
  }
  
  # sp.lineage <- lapply(1:length(sp.id), FUN = function(x, ranks) {
  #   use <- which(!is.na(par.id.df[x,]))
  #   tmp.v <- as.numeric(par.id.df[x, use])
  #   names(tmp.v) <- par.id.names[x, use]
  #   tmp.v <- tmp.v[ranks]
  # })
  
  par.lineage <- data.frame(t(data.frame(sp.lineage)))
  names(par.lineage) <- ranks
  row.names(par.lineage) <- NULL
  par.lineage <- par.lineage[,which(ranks == "species"):ncol(par.lineage)]
  
  ## subspecies lineage
  sp.id <- d$parent.taxid[which(d$rank.parent.taxid == "species")]
  u.sp.id <- unique(sp.id)
  # u.sp.id.l <- split(u.sp.id, ceiling(seq_along(u.sp.id)/(length(u.sp.id)/threads)))
  sub.sp <- foreach::foreach(x = 1:length(u.sp.id)) %dopar% {
    children <- d$taxid[which(d$parent.taxid == u.sp.id[x])]
    out <- data.frame(
      "species" = rep(u.sp.id[x], length(children)),
      "subspecies" = children
    )
    out
  }
  
  sub.sp <- do.call("rbind", sub.sp)
  
  ## and now, with subspecies as parent, where all clades immediately below species were assigned as 'subspecies'
  u.sp.id <- unique(sub.sp$subspecies)
  # u.sp.id.l <- split(u.sp.id, ceiling(seq_along(u.sp.id)/(length(u.sp.id)/threads)))
  sub.sub.sp <- foreach::foreach(x = 1:length(u.sp.id)) %dopar% {
    children <- d$taxid[which(d$parent.taxid == sub.sp$subspecies[x])]
    out <- data.frame(
      "subspecies" = rep(u.sp.id[x], length(children)),
      "sub.subspecies" = children
    )
    out
  }
  sub.sub.sp <- do.call("rbind", sub.sub.sp)
  
  tmp <- merge(sub.sp, sub.sub.sp, by = "subspecies", all = TRUE)
  
  tmp <- merge(tmp, par.lineage, by = "species", all = TRUE)
  
  # head(tmp)
  
  taxid.hierarchy <- tmp[,ranks]
  taxid.hierarchy$root <- "1"
  for(i in 1:ncol(taxid.hierarchy)) {
    taxid.hierarchy[,i] <- as.numeric(taxid.hierarchy[,i])
  }
  
  taxid.hierarchy <- data.table::as.data.table(taxid.hierarchy)
  data.table::setkey(taxid.hierarchy, "species")
  rm(d); gc() 
  rm(i, keep.going, mtch, old.unique.levs, sp.id, u.sp.id, unique.levs, tmp, par.id, par.lineage, sp.lineage, sub.sp, sub.sub.sp)
  rm(par.id.df, par.id.names, j, k, x)
  gc() 
  save(taxid.hierarchy, file = paste0(pc.directory, "/taxid.hierarchy.Rdata"))
  
  readme <- c(
    readme,
    " - taxid.taxname.Rdata, is a data.table of five columns, including taxid, the parent.taxid, the taxid rank, division.id, and the taxid's taxname", '\n', '\n',
    " - taxid.hierarchy.Rdata, is a data.table of several columns, each representing a taxonomic level, ordered from lowest to highest.  Each cell contains an integer taxid value from the NCBI taxonomy database. If hierarchy isn't defined, cell will contain NA. ", '\n',
    paste0(ranks, collapse = '\n', sep = ""),
    '\n', '\n'
  )
  rm(ranks)
  
  ########################
  ## CID to TaxID
  ########################
  download.url <- "https://ftp.ncbi.nlm.nih.gov/pubchem/Target/taxonomy2cid.gz"
  tmp.file <- paste0(tmp.dir, basename(download.url))
  utils::download.file(url = download.url, destfile = tmp.file)
  R.utils::gunzip(tmp.file, remove = FALSE)
  # untar(paste0(gsub(".gz", "", tmp.file)), exdir = tmp.dir)
  
  #### cid.to.taxid.relationship
  d <- data.table::fread(paste0(tmp.dir, basename(tmp.file)), sep="\t", fill=TRUE)
  names(d) <- c("taxid", "cid", "data.source", "source.compound", "source.compound.id", 
                "source.taxonomy", "source.taxonomy.id")
  
  
  
  ## consider whether keep only the natural products type databases
  ## HMDB has too many exogenous compounds, which reduces the value of the lca estimate, and creates associations which are not biologically real (i.e. pesticides are not produced by humans)
  ## FooDB seems to have many genome-predicted records, which seem at best speculative (daidzein is everywhere)
  ## though it seems that many of the records are not visible on the HMDB or FooDB websites any more.  
  ## can use yeast/ecoli pathway data if we want, no need for it to be here. 
  ## individual lab records are too hard to validate.  
  
  # tax.sources <- c("LOTUS - the natural products occurrence database", "The Natural Products Atlas", 
  #                  "KNApSAcK Species-Metabolite Database", "Natural Product Activity and Species Source (NPASS)")
  # d <- d[d$data.source %in% tax.sources,]
  pc.bio.keep.cols <- c("taxid", "cid", "data.source")
  d <- d[, ..pc.bio.keep.cols]
  tmp <- match(d$cid, cid.preferred$cid)
  use <- which(!is.na(tmp))
  if(length(use)>0) {
    d$cid[use] <- cid.preferred$preferred.cid[tmp[use]]
  } 
  cid.taxid <- d
  cid.taxid <- data.table::as.data.table(cid.taxid)
  cid.taxid <- cid.taxid[!duplicated(cid.taxid), ]
  data.table::setkey(cid.taxid, "cid")
  rm(d); gc()
  save(cid.taxid, file = paste0(pc.directory, "/cid.taxid.Rdata"))
  
  ## load("C:/Temp/20250718/cid.taxid.Rdata")
  tax.source.table <- table(cid.taxid$data.source)
  tax.source.table <- data.frame(tax.source.table)
  names(tax.source.table) <- c("data.source", "count")
  utils::write.csv(tax.source.table, file = paste0(pc.directory, "/taxid.sources.csv"), row.names = FALSE)
  rm(tax.source.table)
  
  rm(cid.taxid); gc() # wait to remove cid.taxid, as we will add to it for pathway data
  readme <- c(
    readme,
    " - cid.taxid.Rdata, is a data.table of two columns, headers 'cid'and 'taxid', each integer values. data.table is indexed by 'cid'. ", '\n', '\n'
  )
  
  ########################
  ## cid to pathway
  ########################
  
  pathway.sources <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/sourcetable/all/CSV/?response_type=save&response_basename=PubChemDataSources_all"
  d <- utils::read.csv(pathway.sources)
  d <- d[which(d$Pathway.Count > 0),]
  
  source.name <- d$Source.Name
  
  download.urls <- 
    paste0("https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22pathway%22,%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_pathway_text_", 
           gsub(" ", "%20", source.name),  #  "Plant%20Reactome", 
           "%22,%22where%22:{%22ands%22:[{%22*%22:%22", 
           gsub(" ", "%22},{%22*%22:%22", source.name),   #  "Plant%22},{%22*%22:%22Reactome", 
           "%22},{%22source%22:%22", 
           gsub(" ", "%22},{%22source%22:%22", source.name),   #  "Plant%22},{%22source%22:%22Reactome", 
           "%22}]}}"
    )
  out.files <- paste0(tmp.dir, "/", "pw", 1:length(source.name),".tsv")
  
  # download.urls <- c(
  #   "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22pathway%22,%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_pathway_text_Reactome%22,%22where%22:{%22ands%22:[{%22*%22:%22Reactome%22},{%22source%22:%22Reactome%22}]}}",
  #   "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22pathway%22,%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_pathway_text_PlantCyc%22,%22where%22:{%22ands%22:[{%22*%22:%22PlantCyc%22},{%22source%22:%22PlantCyc%22}]}}",
  #   "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22pathway%22,%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_pathway_text_Plant%20Reactome%22,%22where%22:{%22ands%22:[{%22*%22:%22Plant%22},{%22*%22:%22Reactome%22},{%22source%22:%22Plant%22},{%22source%22:%22Reactome%22}]}}",
  #   "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22pathway%22,%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_pathway_text_PathBank%22,%22where%22:{%22ands%22:[{%22*%22:%22PathBank%22},{%22source%22:%22PathBank%22}]}}",
  #   "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22pathway%22,%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_pathway_text_WikiPathways%22,%22where%22:{%22ands%22:[{%22*%22:%22WikiPathways%22},{%22source%22:%22WikiPathways%22}]}}"
  # )  
  
  for(i in 1:length(out.files)) {
    curl::curl_download(url = download.urls[i], destfile = out.files[i])
  }
  
  #### cid.to.pathway.relationship - read downloaded files
  for(i in 1:length(download.urls)) {
    if(i == 1) {
      d <- data.table::fread(out.files[i])
    } else {
      d <- rbind(d, data.table::fread(out.files[i]))
    }
  }
  
  ## 'there something happenin' here.... 
  ## 'what it is ain't exactly clear....
  d.2 <- lapply(1:nrow(d), FUN = function(i) {
    cids <- d$cids[i]
    cids <- as.numeric(unlist(strsplit(cids, "|", fixed = TRUE)))
    tmp <- do.call("rbind", replicate(length(cids), d[i,], simplify = FALSE))
    tmp$cids <- cids
    tmp
  })
  gc()
  d.2 <- dplyr::bind_rows(d.2)
  d <- d.2
  rm(d.2)
  gc()
  pc.bio.keep.cols <- c("pwacc", "name", "pwtype", "source", 
            "taxid", "taxname", "cids")
  d <- d[, ..pc.bio.keep.cols]
  names(d) <- c("pwid", "name", "pwtype", "source", 
                "taxid", "taxname", "cid")
  d <- d[-which(is.na(d$cid)),]
  tmp <- match(d$cid, cid.preferred$cid)
  use <- which(!is.na(tmp))
  if(length(tmp)>0) {
    d$cid[use] <- cid.preferred$preferred.cid[tmp[use]]
  }
  cid.pwid <- d
  rm(d); gc()
  gc()
  data.table::setkey(cid.pwid, "cid")
  gc()
  save(cid.pwid, file = paste0(pc.directory, "/cid.pwid.Rdata"))
  
  ## load("C:/Temp/20250718/cid.pwid.Rdata")
  pw.source.table <- table(cid.pwid$source)
  pw.source.table <- data.frame(pw.source.table)
  names(pw.source.table) <- c("data.source", "count")
  utils::write.csv(pw.source.table, file = paste0(pc.directory, "/pw.sources.csv"), row.names = FALSE)
  rm(pw.source.table)
  
  rm(cid.pwid)
  gc()
  readme <- c(
    readme,
    " - cid.pwid.Rdata, is a data.table of seven columns, 'pwid' (integer), 'name' (character), 'pwtype' (character), 'source' (character), 'taxid' (ingeter), 'taxname' (character), and 'cid' (integer). data.table is indexed by 'cid'. ", '\n', '\n'
  )
  
  # sp.spec <- which(!is.na(cid.pwid$taxid))
  # cid.taxid.2 <- data.frame(
  #   taxid = cid.pwid$taxid[sp.spec],
  #   cid = cid.pwid$cid[sp.spec]
  # )
  # cid.taxid.2 <- data.table::as.data.table(cid.taxid.2)
  # cid.taxid.2 <- cid.taxid.2[!duplicated(cid.taxid.2), ]
  # 
  # ## conserved plantcyc pathways
  # con.path <- which(cid.pwid$pwtype == 'conserved' & cid.pwid$source == "PlantCyc")
  # cids <- unique(cid.pwid$cid[con.path])
  # taxids <- 33090
  # cid.taxid.3 <- expand.grid(taxid = taxids, cid = cids)
  # 
  # cid.taxid <- rbind(
  #   cid.taxid,
  #   cid.taxid.2
  # )
  
  
  ########################
  ## SID to Data Source
  ########################
  download.url <- "https://ftp.ncbi.nlm.nih.gov/pubchem/Substance/Extras/SID-Map.gz"
  tmp.file <- paste0(tmp.dir, basename(download.url))
  utils::download.file(url = download.url, destfile = tmp.file)
  R.utils::gunzip(tmp.file, remove = FALSE)
  gc()
  
  #### sid.to.data.source.relationship
  d <- data.table::fread(paste0(tmp.dir, basename(tmp.file)), sep="\t", fill=TRUE) 
  # dim(d)
  
  names(d) <- c("sid", "source", "source.id", "cid")
  gc()
  
  # bio.cid <- bio.sources %in% d[,2]
  
  tmp <- match(d$cid, cid.preferred$cid)
  use <- which(!is.na(tmp))
  if(length(tmp)>0) {
    d$cid[use] <- cid.preferred$preferred.cid[tmp[use]]
  }
  cid.sid <- d
  rm(d); gc()
  data.table::setkey(cid.sid, "cid")
  save(cid.sid, file = paste0(pc.directory, "/cid.sid.Rdata"))
  
  ## load("C:/Temp/20250718/cid.sid.Rdata")
  all.source.table <- table(cid.sid$source)
  all.source.table <- data.frame(all.source.table)
  names(all.source.table) <- c("data.source", "count")
  utils::write.csv(all.source.table, file = paste0(pc.directory, "/all.sources.csv"), row.names = FALSE)
  rm(all.source.table)
  
  rm(cid.sid); gc()
  readme <- c(
    readme,
    " - cid.sid.Rdata, is a data.table of four columns, 'sid' (integer), 'source' (character), 'source.id' (character), and 'cid' (integer). data.table is indexed by 'cid'. ", '\n', '\n'
  )
  
  ########################
  ## CID to Synonym
  ########################
  download.url <- "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Synonym-unfiltered.gz"
  tmp.file <- paste0(tmp.dir, basename(download.url))
  utils::download.file(url = download.url, destfile = tmp.file)
  R.utils::gunzip(tmp.file, remove = FALSE)
  # untar(paste0(gsub(".gz", "", tmp.file)), exdir = tmp.dir)
  
  #### cid.to.synonyms.relationship
  d <- data.table::fread(gsub(".gz", "", paste0(tmp.dir, basename(tmp.file))), quote = "")
  names(d) <- c("cid", "synonym")
  cid.synonym <- as.vector(d$cid)
  cid.names <- gsub(" ", ".", d$synonym)
  names(cid.synonym) <- cid.names
  
  cas.like <- which(stringr::str_detect(d$synonym, "[:digit:]{2,7}[-][:digit:]{2}[-][:digit:]{1}"))
  cid.cas <- d[cas.like,]
  names(cid.cas) <- c("cid", "cas")
  data.table::setkey(cid.cas, "cid")
  save(cid.cas, file = paste0(pc.directory, "/cid.cas.Rdata"))
  rm(cid.cas); rm(cas.like)
  
  cid.synonym <- d
  rm(d)
  data.table::setkey(cid.synonym, "cid")
  save(cid.synonym, file = paste0(pc.directory, "/cid.synonym.Rdata"))
  cid.synonyms <- cid.synonym[, list(synonym = list(synonym)),
                              by = list(cid)]
  rm(cid.synonym); gc()
  data.table::setkey(cid.synonyms, "cid")
  save(cid.synonyms, file =  paste0(pc.directory, "/cid.synonyms.Rdata"))
  rm(cid.synonyms); gc()
  readme <- c(
    readme,
    " - cid.cas.Rdata, is a data.table of two columns, 'cid' (integer), 'cas' (character).  Note that CAS is extracted from 'synonymns' using string pattern matching.  It is possible there are errors. data.table is indexed by 'cid'. ", '\n', '\n',
    " - cid.synonym.Rdata, is a data.table of two columns, 'cid' (integer), 'synonym' (character). data.table is indexed by 'cid'. ", '\n', '\n',
    " - cid.synonyms.Rdata, is a data.table of two columns, 'cid' (integer), 'synonym' (character), where synonym is a list of synonyms for each cid. data.table is indexed by 'cid'. ", '\n', '\n'
  )
  
  ## record data source table
  all.sources <- utils::read.csv("https://pubchem.ncbi.nlm.nih.gov/rest/pug/sourcetable/substance/CSV")
  all.sources$source.pubchem.url <- paste0("https://pubchem.ncbi.nlm.nih.gov/source/", all.sources$DSN)
  save(all.sources, file =  paste0(pc.directory, "/all.sources.Rdata"))
  readme <- c(
    readme,
    " - all.sources.Rdata, is a data.table several columns, reported as directly downloaded from Pubchem via https://pubchem.ncbi.nlm.nih.gov/rest/pug/sourcetable/substance/CSV. ", '\n', '\n'
  )
  
  
  fileConn<-file(paste0(pc.directory, '/readme.txt'))
  writeLines(readme, fileConn)
  close(fileConn)
  
  if(rm.tmp.files) {
    unlink(paste0(pc.directory, "/tmp/"), recursive=TRUE)
  }
  

  
  
  message(' -- finished', '\n')
  
  
  return(readme)
}


