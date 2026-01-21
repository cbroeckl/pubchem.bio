#' export.ComboundDb
#'
#' export pubchem.bio pc.bio syle data.table to a CompoundDb database
#' @details takes output from 'build.pubchem.bio' or 'build.taxon.metabolome' functions, reformatting, and converting to a CompoundDb object and database. see see vignette('create-compounddb', package = 'CompoundDb') for more details on the CompoundDb package
#' @param pc.bio.object input data.table, generated from 'build.pubchem.bio' or 'build.taxon.metabolome' functions
#' @param pc.directory directory from which to load pubchem .Rdata files.  alternatively, provide  R data.tables for ALL cid._property_.object options defined below.
#' @param path valid file path specifying the directory to which the database will be saved.
#' @param dbFile valid file name - use .sqlite extension. if NULL, will set to 'pubchem.bio.sqlite'
#' @param source default = NULL.  if NULL, will set to 'pubchem.bio'
#' @param source_version default = NULL.  if NULL, will set to pubchem.bio version.
#' @param url default = NULL. if NULL, will set to "https://cran.r-project.org/web/packages/pubchem.bio/index.html"
#' @param organism character.  default = NA_character_.  if this is a species specific database, set to species name or NCBI taxonomy ID number.
#' @param source_date default = NULL. if NULL, will set to current date via format(Sys.Date(), format = "%Y-%m-%d")
#' @param get.synonyms logical. default = FALSE.  When TRUE, will load the 'cid.synonyms.Rdata' file from the pc.directory, which must be defined. database product will have all synonyms included. If FALSE, synonyms will be an empty list.
#' @return 'CompDb' from CompoundDb package.  see vignette('create-compounddb', package = 'CompoundDb')
#' @author Corey Broeckling
#'
#' @export
#'

export.ComboundDb <- function(
    pc.bio.object = NULL,
    pc.directory = NULL,
    path = NULL,
    dbFile = NULL,
    source = NULL,
    source_version = NULL,
    url = NULL,
    source_date = NULL,
    organism = NA_character_,
    get.synonyms = FALSE
) {

  ## check that pc.bio.object has all the correct headers
  ## change names a bit for compatibility
  names(pc.bio.object)[which(names(pc.bio.object) == "monoisotopic.mass")] <- 'exactmass'

  ## add 'compound_id' column as pubchem cid
  pc.bio.object <- data.table::data.table(
    pc.bio.object,
    compound_id = pc.bio.object$cid
  )

  names(pc.bio.object) <- gsub(".", "_", names(pc.bio.object), fixed = TRUE)

  ## check names for completeness of required
  if(!all(c("compound_id", "name", "inchi", "inchikey", "formula", "exactmass") %in% names(pc.bio.object))) {
    stop('the column names "compound_id", "name", "inchi", "inchikey", "formula", and "exactmass" are required', '\n')
  }

  if(is.null(dbFile)) dbFile = "pubchem.bio.sqlite"
  if(is.null(source)) source = "pubchem.bio"
  if(is.null(source_version)) source_version = as.character(utils::packageVersion('pubchem.bio'))
  if(is.null(url)) url = "https://cran.r-project.org/web/packages/pubchem.bio/index.html"
  if(is.null(source_date)) source_date = format(Sys.Date(), format = "%Y-%m-%d")
  if(is.null(organism)) organism = NA_character_

  ## either get the pubchem.bio synonyms or make empty list of synonyms
  if(get.synonyms) {

    synonym <- NA
    load(paste0(pc.directory, "/cid.synonym.Rdata"))
    cid.synonym <- cid.synonym
    data.table::setkey(cid.synonym, "cid")
    .datatable.aware = TRUE
    cid.synonyms <- cid.synonym[, list('synonym' = list(synonym)),
                                by = list(cid)]
    rm(cid.synonym); gc()
    m <- match(pc.bio.object$cid, cid.synonyms$cid)
    synonyms <- cid.synonyms$synonym[m]
    rm(m); rm(cid.synonyms); gc()
  } else {
    synonyms <- as.list(rep(NA_character_, nrow(pc.bio.object)))
  }

  ##
  pc.bio.object <- data.table::data.table(
    pc.bio.object,
    synonyms
  )

  db.con <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = paste0(path, "/", dbFile))
  on.exit(
    if(RSQLite::dbIsValid(db.con)) {
      RSQLite::dbDisconnect(db.con)
    }
  )

  metad <- CompoundDb::make_metadata(source = source,
                                     source_version = source_version,
                                     url = url,
                                     source_date = source_date,
                                     organism = organism)

  db_file <- CompoundDb::createCompDb(pc.bio.object, metadata = metad, path = path,
                                      dbFile = dbFile)
  cdb <- CompoundDb::CompDb(db_file, flags = RSQLite::SQLITE_RW)
  cdb

}

# cpdb <- export.ComboundDb(pc.bio.object = pubchem.bio)
