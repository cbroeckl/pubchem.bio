#' build.element.count
#'
#' takes as input an pubchem.bio data.table (generally produced by 'build.pubchem.bio' or 'build.taxon.metabolome') and removes inorganic compounds (those without any carbon).  
#' @details utilizes downloaded and properly formatted local pubchem data created by 'build.pubchem.bio' function, removes inorganic compounds. note that for any formula which have charge noted '+' or '-', the charge is removed for tabulation of the element count. 
#' @param pc.directory directory from which to load pubchem .Rdata files
#' @param pubchem.bio.object R data.table, generally produced by build.pubchem.bio; preferably, define pc.directory
#' @param remove.homonuclear.molecules logical. default = TRUE.  should molecules with one element (or fewer) be removed?
#' @param remove.inorganics logical. default = TRUE.  Should inorganic molecules (no carbon) be removed?  
#' @return a data.table, same as input, except inorganics will have been removed.  
#' @author Corey Broeckling
#' 
#' @export
#' 
#'
#'
build.element.count <- function(
    pc.directory = NULL,
    pubchem.bio.object = NULL,
    remove.homonuclear.molecules = TRUE, 
    remove.inorganics = TRUE
    ) {
  
  loadNamespace("data.table")
  .datatable.aware = TRUE
  
  if(is.null(pubchem.bio.object)) {
    load(paste0(pc.directory, "/pc.bio.Rdata"))
    pc.bio <- pc.bio
  } else {
    pc.bio <- pubchem.bio.object
  }
  
  forms <- pc.bio$formula
  if(any(grepl("+", forms, fixed = TRUE))) {
    forms <- sapply(1:length(forms), FUN = function(x) {unlist(strsplit(forms[[x]], "+", fixed = TRUE))[1]})
  }
  
  if(any(grepl("-", forms, fixed = TRUE))) {
    forms <- sapply(1:length(forms), FUN = function(x) {unlist(strsplit(forms[[x]], "-", fixed = TRUE))[1]})
  }
  
  suppressWarnings(elems <- MetaboCoreUtils::countElements(forms))
  elem.l <- sapply(1:length(elems), FUN = function(x) {length(elems[[x]])})
  if(any(elem.l == 0)) {
    elem.rep <- which(elem.l == 0)
    for(i in elem.rep) {
      elems[[i]] <- c("C" = 0)
    }
  }
  
  if(remove.homonuclear.molecules) {
    pc.bio <- pc.bio[which(elem.l >= 2),]
    forms <- forms[which(elem.l >= 2)]
    elems <- elems[which(elem.l >= 2)]
    elem.l <- elem.l[which(elem.l >= 2)]
  }
  
  elems <- data.table::rbindlist(lapply(elems, as.data.frame.list), fill=TRUE)
  no.name <- grep("integer", names(elems))
  if(length(no.name) > 0) {
    elems <- elems[, grep("integer", names(elems)):=NULL]
  }
  data.table::setnafill(elems, fill = 0)
  if(remove.inorganics) {
    inorg <- which(elems$C == 0)
    if(length(inorg) > 0) {
      pc.bio <- pc.bio[-inorg,]
      elems <- elems[-inorg,]
      col.sums <- colSums(elems)
      keep.cols <- which(col.sums > 0)
      elems <- elems[,..keep.cols]
    }
  }
  
  out <- cbind(pc.bio, elems)
  return(out)
}

