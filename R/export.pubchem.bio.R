#' export.pubchem.bio
#'
#' export pubchem.bio pc.bio syle data.table to tab delimited text file for import into other programs.  all columns exported.    
#' @details takes output from 'build.pubchem.bio' or 'build.taxon.metabolome' functions, reformatting, and exporting to input format suitable for MSFinder.  
#' @param pc.bio.object input data.table, generated from 'build.pubchem.bio' or 'build.taxon.metabolome' functions
#' @param export.file.name valid file path and name.  Extension should be listed as '.tsv'.
#' @return nothing - file written to disk. 
#' @author Corey Broeckling
#' 
#' @export 
#' 

export.pubchem.bio <- function(
    pc.bio.object = NULL,
    export.file.name = NULL
) {
  
  utils::write.table(pc.bio.object, file = export.file.name, sep = "\t", col.names=TRUE, row.names = FALSE)
  
}

# export.pubchem.bio(pc.bio.object = pc.bio.sub, export.file.name = "C:/Temp/20250703/pc.bio.full.pubchem.bio.tsv")
