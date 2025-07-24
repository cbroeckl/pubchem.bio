#' get.pubchem.ftp
#'
#' first step to building a local selective, biologically focused, pubchem data repository focused on metabolomics informatics 
#' @details this function downloads and unzips files from pubchem and NCBI taxonomy FTP sites as a first step in building a local metabolomics repository.
#' @param pc.directory character. directory to which data will be saved
#' @param timeout numeric.  timeout setting for FTP download.  setting options(timeout) value too small will generate errors for large files. default = 50000.
#' @param rm.tmp.files logical.  should temporary files be removed after completion of download and parsing?  Default = TRUE. 
#' @return nothing.  all data are saved to disk for later loading
#' @author Corey Broeckling
#' 
#' @examples
#' /donttest{
#' my.dir <- paste0(getwd(), "/", format(Sys.Date(), "%Y%m%d"))
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
    rm.tmp.files = TRUE
) {
  
  ## setup directory structure
  if(is.null(pc.directory)) pc.directory <- paste0(getwd(), "/", format(Sys.time(), "%Y%m%d"), "/")
  pc.directory <- suppressWarnings(normalizePath(pc.directory))
  if(!dir.exists(pc.directory)) dir.create(pc.directory)
  tmp.dir <- suppressWarnings(paste0(pc.directory, "/tmp/"))
  dir.create(tmp.dir)
  
  message(" -- writing data to", pc.directory, '\n')
  
  readme <- c(
    "Data collection derived from pubchem, primarily using the 'compound' FTP data.", '\n',
    "generation date:", date(), '\n')
  
  ## set download timeout to prevent errors
  oldpar <- par(no.readonly = TRUE) 
  on.exit(par(oldpar)) 
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
  
  
  ########################
  ## CID to parent CID
  ########################
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
  Sys.sleep(2)
  
  readme <- c(
    readme,
    " - cid.parent.Rdata, is a data.table of two columns, headers 'cid' and 'parent.cid', each integer values. data.table is indexed by 'cid'. ", '\n', '\n'
  )
  
  ########################
  ## CID to InChIKey
  ########################
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
  cid.inchikey <- d[,c(1,3)]
  data.table::setkey(cid.inchikey, "cid")
  save(cid.inchikey, file = paste0(pc.directory, "/cid.inchikey.Rdata"))
  rm(cid.inchikey);   gc()
  Sys.sleep(2)
  
  cid.inchi <- d[,c(1,2)]
  data.table::setkey(cid.inchi, "cid")
  save(cid.inchi, file = paste0(pc.directory, "/cid.inchi.Rdata"))
  rm(cid.inchi); gc()
  Sys.sleep(2)
  rm(d); gc()
  Sys.sleep(2)
  
  readme <- c(
    readme,
    " - cid.inchikey.Rdata, is a data.table of two columns, headers 'cid'and 'inchikey', with integer and character values, respectively. data.table is indexed by 'cid'. ", '\n', '\n',
    " - cid.inchi.Rdata, is a data.table of two columns, headers 'cid'and 'inchi', with integer and character values, respectively. data.table is indexed by 'cid'. ", '\n', '\n'
  )
  
  ########################
  ## CID to SMILES
  ########################
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
  Sys.sleep(2)
  data.table::setkey(cid.smiles, "cid")
  save(cid.smiles, file = paste0(pc.directory, "/cid.smiles.Rdata"))
  rm(cid.smiles); gc()
  Sys.sleep(2)
  
  readme <- c(
    readme,
    " - cid.smiles.Rdata, is a data.table of two columns, headers 'cid'and 'smiles', with integer and character values, respectively. data.table is indexed by 'cid'. ", '\n', '\n'
  )
  
  ########################
  ## CID to Title
  ########################
  download.url <- "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Title.gz"
  tmp.file <- paste0(tmp.dir, basename(download.url))
  utils::download.file(url = download.url, destfile = tmp.file)
  R.utils::gunzip(tmp.file, remove = FALSE, destname = paste0(tmp.dir, gsub(".gz", "", basename(tmp.file))))
  # untar(paste0(gsub(".gz", "", tmp.file)), exdir = tmp.dir)
  
  #### cid.to.title.relationship
  d <- data.table::fread(paste0(tmp.dir, basename(tmp.file)))
  names(d) <- c("cid", "title")
  data.table::setkey(d, "cid")
  tmp <- match(d$cid, cid.preferred$cid)
  use <- which(!is.na(tmp))
  if(length(use)>0) {
    d$cid[use] <- cid.preferred$preferred.cid[tmp[use]]
  }  
  cid.title <- d
  rm(d); gc()
  Sys.sleep(2)
  data.table::setkey(cid.title, "cid")
  save(cid.title, file = paste0(pc.directory, "/cid.title.Rdata"))
  rm(cid.title); gc
  Sys.sleep(2)
  
  readme <- c(
    readme,
    " - cid.title.Rdata, is a data.table of two columns, headers 'cid'and 'title', with integer and character values, respectively. data.table is indexed by 'cid'. ", '\n', '\n'
  )
  
  ########################
  ## CID to PMID
  ########################
  download.url <- "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-PMID.gz"
  tmp.file <- paste0(tmp.dir, basename(download.url))
  utils::download.file(url = download.url, destfile = tmp.file)
  R.utils::gunzip(tmp.file, remove = FALSE, destname = paste0(tmp.dir, gsub(".gz", "", basename(tmp.file))))
  # untar(paste0(gsub(".gz", "", tmp.file)), exdir = tmp.dir)
  
  #### cid.to.pmid.relationship
  d <- data.table::fread(paste0(tmp.dir, basename(tmp.file)))
  d <- d[,c(1,2)]
  names(d) <- c("cid", "pmid")
  data.table::setkey(d, "cid")
  tmp <- match(d$cid, cid.preferred$cid)
  use <- which(!is.na(tmp))
  if(length(use)>0) {
    d$cid[use] <- cid.preferred$preferred.cid[tmp[use]]
  } 
  cid.pmid <- d
  rm(d); gc()
  Sys.sleep(2)
  data.table::setkey(cid.pmid, "cid")
  save(cid.pmid, file = paste0(pc.directory, "/cid.pmid.Rdata"))
  
  cid.pmid.ct <- table(cid.pmid$cid)
  cid.pmid.ct <- data.table::data.table(cid.pmid.ct)
  names(cid.pmid.ct) <- c('cid', 'pmid.ct')
  cid.pmid.ct$cid <- as.integer(as.numeric(cid.pmid.ct$cid))
  data.table::setkey(cid.pmid.ct, "cid")
  save(cid.pmid.ct, file = paste0(pc.directory, "/cid.pmid.ct.Rdata"))
  
  rm(cid.pmid); gc()
  Sys.sleep(2)
  
  readme <- c(
    readme,
    " - cid.pmid.Rdata, is a data.table of two columns, headers 'cid'and 'pmid', with integer and character values, respectively. data.table is indexed by 'cid'. ", '\n', '\n',
    " - cid.pmid.ct.Rdata, is a data.table of two columns, headers 'cid'and 'pmid.ct', each with integer values.  'pmid.ct' represents the number of pubmed references for that compound. data.table is indexed by 'cid'. ", '\n', '\n'
  )
  
  ########################
  ## CID to Mass
  ########################
  download.url <- "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Mass.gz"
  tmp.file <- paste0(tmp.dir, basename(download.url))
  utils::download.file(url = download.url, destfile = tmp.file)
  R.utils::gunzip(tmp.file, remove = FALSE, destname = paste0(tmp.dir, gsub(".gz", "", basename(tmp.file))))
  # untar(paste0(gsub(".gz", "", tmp.file)), exdir = tmp.dir)
  
  #### cid.to.mass.relationship
  d <- data.table::fread(paste0(tmp.dir, basename(tmp.file)))
  names(d) <- c("cid", "formula", "monoisotopic.mass", "exact.mass")
  data.table::setkey(d, "cid")
  tmp <- match(d$cid, cid.preferred$cid)
  use <- which(!is.na(tmp))
  if(length(use)>0) {
    d$cid[use] <- cid.preferred$preferred.cid[tmp[use]]
  } 
  
  cid.formula <- d[,c(1,2)]
  data.table::setkey(cid.formula, "cid")
  save(cid.formula, file = paste0(pc.directory, "/cid.formula.Rdata"))
  rm(cid.formula); gc()
  Sys.sleep(2)
  cid.monoisotopic.mass <- d[,c(1,3)]
  data.table::setkey(cid.monoisotopic.mass, "cid")
  save(cid.monoisotopic.mass, file = paste0(pc.directory, "/cid.monoisotopic.mass.Rdata"))
  rm(cid.monoisotopic.mass); gc()
  Sys.sleep(2)
  cid.accurate.mass <- d[,c(1,4)]
  data.table::setkey(cid.accurate.mass, "cid")
  save(cid.accurate.mass, file = paste0(pc.directory, "/cid.accurate.mass.Rdata"))
  rm(cid.accurate.mass); gc()
  Sys.sleep(2)
  rm(d); gc()
  Sys.sleep(2)
  
  readme <- c(
    readme,
    " - cid.formula.Rdata, is a data.table of two columns, headers 'cid'and 'formula', with integer and character values, respectively. data.table is indexed by 'cid'. ", '\n', '\n',
    " - cid.monoisotopic.mass.Rdata, is a data.table of two columns, headers 'cid'and 'monoisotopic.mass', with integer and character values, respectively. data.table is indexed by 'cid'. ", '\n', '\n'
  )
  
  ########################
  ## CID to MeSH name
  ########################
  download.url <- "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-MeSH"
  tmp.file <- paste0(tmp.dir, basename(download.url))
  utils::download.file(url = download.url, destfile = tmp.file)
  # R.utils::gunzip(tmp.file, remove = FALSE)
  # untar(paste0(gsub(".gz", "", tmp.file)), exdir = tmp.dir)
  
  #### cid.to.mesh.name relationship
  d <- readLines(paste0(tmp.dir, basename(tmp.file)))
  d <- lapply(1:length(d), FUN = function(x) {
    tmp <- unlist(strsplit(d[x], '\t', fixed = TRUE))
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
  readme <- c(
    readme,
    " - cid.mesh.name.Rdata, is a data.table of two columns, headers 'cid'and 'mesh.name', with integer and character values, respectively. data.table is indexed by 'cid'. ", '\n', '\n'
  )
  
  
  #### separate file for mesh name to mesh function
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
  Sys.sleep(2)
  rm(cid.mesh.name); gc()
  Sys.sleep(2)
  readme <- c(
    readme,
    " - cid.mesh.function.Rdata, is a data.table of two columns, headers 'cid'and 'mesh.function', with integer and character values, respectively. data.table is indexed by 'cid'. ", '\n', '\n'
  )
  Sys.sleep(2)
  
  
  ########################
  ## TaxID to taxonomy
  ########################
  download.url <- "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
  tmp.file <- paste0(tmp.dir, download.url)
  utils::download.file(url = download.url, destfile = paste0(tmp.dir, basename(tmp.file)))
  R.utils::gunzip(paste0(tmp.dir, basename(tmp.file)), remove = FALSE, 
                  destname = gsub(".gz", "", paste0(tmp.dir, basename(tmp.file))))
  utils::untar(gsub(".gz", "", paste0(tmp.dir, basename(tmp.file))), exdir = tmp.dir)
  
  #### taxid.to.rank.relationship
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
  d <- d[,c("taxid", "parent.taxid", "rank"),]
  gc()
  Sys.sleep(2)
  
  
  d <- data.table::as.data.table(d)
  data.table::setkey(d, "taxid")
  
  ## check.taxid <- 1280
  ## check.lineage <- c(1280, 1279, 90964, 1385, 91061, 1239, 2, 1)
  ## d[d$taxid %in% check.lineage,]
  ##  1280  1279 90964  1385 91061  1239     2     1
  ### try replacing all parent ids that are not in the rank list
  ## with their increasing parentage until you hit a rank list level
  #  consider adding 'strain', 'isolate', 'subspecies, - 
  ranks <- c("species", "genus", "subtribe", "tribe", "subfamily", "family", "superfamily", "infraorder", "suborder", "order", "superorder", "infraclass", "subclass", "class", "superclass", "subphylum", "phylum", "subkingdom", "kingdom", "superkingdom")
  d$rank.parent.taxid <- d$rank[match(d$parent.taxid, d$taxid)]
  missing.old <- length(which(!d$rank.parent.taxid %in% ranks))
  missing.new <- 1 + missing.old
  while(missing.new != missing.old) {
    missing.old <- missing.new
    do <- which(!(d$rank.parent.taxid %in% ranks))
    old.parent <- d$parent.taxid[do]
    repl.index <- match(old.parent, d$taxid)
    new.parent <- d$parent.taxid[repl.index]
    new.parent.rank <- d$rank.parent.taxid[repl.index]
    d$parent.taxid[do] <- new.parent
    d$rank.parent.taxid[do] <- new.parent.rank
    parent.rank <- d$rank[match(d$parent.taxid, d$taxid)]
    missing.new <- length(which(!parent.rank %in% ranks))
    # message(missing.new, '\n')
    # d[d$taxid %in% check.lineage,]
  }
  # d.orig <- d
  ## now each taxid to parent relationship is defined with only specified ranks allowed
  ## need to reorganize into full hierarchy table now, or create network
  
  ## data.frame
  ## for each level of rank, find parent id and parent taxid rank
  #### for each level of parent taxid rank, 
  #### insert the parent taxid in appropriate column
  tar.rank <- which(d$rank == "species")
  taxid.hierarchy <- matrix(nrow = length(tar.rank), ncol = length(ranks)) 
  dimnames(taxid.hierarchy)[[2]] <- ranks
  taxid.hierarchy[,"species"] <- as.numeric(as.vector(d$taxid[tar.rank]))
  child.taxid <- as.numeric(d$taxid)
  child.taxid.rank <- d$rank
  child.taxid.ind <- match(child.taxid.rank, ranks)
  parent.taxid <- as.numeric(d$parent.taxid)
  parent.taxid.rank <- d$rank.parent.taxid
  parent.taxid.ind <- match(parent.taxid.rank, ranks)
  tab <- data.frame(child.taxid, child.taxid.rank, child.taxid.ind, parent.taxid, parent.taxid.rank, parent.taxid.ind)
  # head(tab[do.d,], 10)
  for(i in 1:ncol(taxid.hierarchy)) {
    for(j in 2:ncol(taxid.hierarchy)) {
      do.d <- which(child.taxid.ind == i & parent.taxid.ind == j)
      do.th.row <- which(taxid.hierarchy[,i] %in% child.taxid[do.d])
      child.vals <- match(taxid.hierarchy[do.th.row,i], child.taxid)
      parent.vals <- parent.taxid[child.vals]
      # length(do.th.row); length(child.vals); length(parent.vals); length(child.taxid); length(parent.taxid)
      # if(length(rm.na) > 0) {
      #   do.d <- do.d[-rm.na]
      #   do.th.row <- do.th.row[-rm.na]
      # }
      taxid.hierarchy[do.th.row,j] <- parent.vals  # match(as.numeric(parent.taxid[do.d]), taxid.hierarchy[,i])
    }
  }
  
  ## tried this option for display of hierarchy as a network, where NA values were causing trouble
  ## but abandoned that - didn't seem worth the extra effort.  
  ## now replace all empty level values (NA) with a placeholder value to build a network in which 
  ## there is a constant number of levels. 
  ## and add a 'root' node to the end column to ensure all organisms share a common ancestor
  
  ## need to replace NAs one column at a time
  ## for each column, NA value should be identical if preceding value is identical
  # for(i in 2:ncol(taxid.hierarchy)) {
  #   nas <- which(is.na(taxid.hierarchy[, i]))
  #   if(length(nas) == 0) next
  #   to.replace <- taxid.hierarchy[, (i-1)][nas]
  #   replace.with <- max(max(taxid.hierarchy, na.rm = TRUE), 900000000) + as.integer(as.factor(as.character(to.replace)))
  #   taxid.hierarchy[nas, i] <- replace.with
  # }
  
  taxid.hierarchy <- data.table::as.data.table(taxid.hierarchy)
  taxid.hierarchy$root <- rep(1, nrow(taxid.hierarchy))
  data.table::setkey(taxid.hierarchy, "species")
  rm(d); gc() 
  Sys.sleep(2)
  save(taxid.hierarchy, file = paste0(pc.directory, "/taxid.hierarchy.Rdata"))
  
  readme <- c(
    readme,
    " - taxid.hierarchy.Rdata, is a data.table of 21 columns, each representing a taxonomic level.  Each cell contains an integer taxid value from the NCBI taxonomy database. If hierarchy isn't defined, cell will contain NA. ",
    '\n', "  species",
    '\n', "  genus",
    '\n', "  subtribe",
    '\n', "  tribe",
    '\n', "  subfamily",
    '\n', "  family",
    '\n', "  superfamily",
    '\n', "  infraorder",
    '\n', "  suborder",
    '\n', "  order",
    '\n', "  superorder",
    '\n', "  infraclass",
    '\n', "  subclass",
    '\n', "  class",
    '\n', "  superclass",
    '\n', "  subphylum",
    '\n', "  phylum",
    '\n', "  subkingdom",
    '\n', "  kingdom",
    '\n', "  superkingdom",
    '\n', "  root",
    '\n', '\n'
  )
  
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
  
  d <- d[,c(1:3)]
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
  Sys.sleep(2)
  save(cid.taxid, file = paste0(pc.directory, "/cid.taxid.Rdata"))
  rm(cid.taxid); gc() # wait to remove cid.taxid, as we will add to it for pathway data
  Sys.sleep(2)
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
  
  # download.urls <- c(
  #   "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22pathway%22,%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_pathway_text_Reactome%22,%22where%22:{%22ands%22:[{%22*%22:%22Reactome%22},{%22source%22:%22Reactome%22}]}}",
  #   "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22pathway%22,%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_pathway_text_PlantCyc%22,%22where%22:{%22ands%22:[{%22*%22:%22PlantCyc%22},{%22source%22:%22PlantCyc%22}]}}",
  #   "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22pathway%22,%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_pathway_text_Plant%20Reactome%22,%22where%22:{%22ands%22:[{%22*%22:%22Plant%22},{%22*%22:%22Reactome%22},{%22source%22:%22Plant%22},{%22source%22:%22Reactome%22}]}}",
  #   "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22pathway%22,%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_pathway_text_PathBank%22,%22where%22:{%22ands%22:[{%22*%22:%22PathBank%22},{%22source%22:%22PathBank%22}]}}",
  #   "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22pathway%22,%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_pathway_text_WikiPathways%22,%22where%22:{%22ands%22:[{%22*%22:%22WikiPathways%22},{%22source%22:%22WikiPathways%22}]}}"
  # )  
  
  #### cid.to.pathway.relationship
  for(i in 1:length(download.urls)) {
    if(i == 1) {
      d <- data.table::fread(download.urls[i])
    } else {
      d <- rbind(d, data.table::fread(download.urls[i]))
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
  d <- d[,c("pwacc", "name", "pwtype", "source", 
            "taxid", "taxname", "cids")]
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
  Sys.sleep(2)
  gc()
  data.table::setkey(cid.pwid, "cid")
  gc()
  save(cid.pwid, file = paste0(pc.directory, "/cid.pwid.Rdata"))
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
  Sys.sleep(2)
  gc()
  
  #### sid.to.data.source.relationship
  d <- data.table::fread(paste0(tmp.dir, basename(tmp.file)), sep="\t", fill=TRUE) 
  # dim(d)
  
  names(d) <- c("sid", "source", "source.id", "cid")
  gc()
  Sys.sleep(2)
  
  # bio.cid <- bio.sources %in% d[,2]
  
  tmp <- match(d$cid, cid.preferred$cid)
  use <- which(!is.na(tmp))
  if(length(tmp)>0) {
    d$cid[use] <- cid.preferred$preferred.cid[tmp[use]]
  }
  cid.sid <- d
  rm(d); gc()
  Sys.sleep(2)
  data.table::setkey(cid.sid, "cid")
  save(cid.sid, file = paste0(pc.directory, "/cid.sid.Rdata"))
  rm(cid.sid); gc()
  Sys.sleep(2)
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
  readme <- c(
    readme,
    " - cid.cas.Rdata, is a data.table of two columns, 'cid' (integer), 'cas' (character).  Note that CAS is extracted from 'synonymns' using string pattern matching.  It is possible there are errors. data.table is indexed by 'cid'. ", '\n', '\n',
    " - cid.synonym.Rdata, is a data.table of two columns, 'cid' (integer), 'synonym' (character). data.table is indexed by 'cid'. ", '\n', '\n'
  )
  
  sink(paste0(pc.directory, '/readme.txt'))
  message(readme)
  sink()
  
  if(rm.tmp.files) {
    unlink(paste0(pc.directory, "tmp/"), recursive=TRUE)
  }
  
  message(' -- finished', '\n')
  
  return(readme)
}


