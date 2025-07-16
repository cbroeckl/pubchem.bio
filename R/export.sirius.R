#' export.sirius
#'
#' export pubchem.bio pc.bio syle data.table to format suitable for sirius input.  
#' @details takes output from 'build.pubchem.bio' or 'build.taxon.metabolome' functions, reformatting, and exporting to input format suitable for Sirius.  
#' @param pc.bio.object input data.table, generated from 'build.pubchem.bio' or 'build.taxon.metabolome' functions
#' @param export.file.name valid file path and name.  Extension should be listed as '.tsv'.
#' @return nothing - file written to disk. 
#' @author Corey Broeckling
#' 
#' @export 
#' 

export.sirius <- function(
    pc.bio.object = NULL,
    export.file.name = NULL
) {
  
  out <- data.frame(
    smiles = pc.bio.object$smiles,
    id = paste0(pc.bio.object$cid),
    name = pc.bio.object$name
  )
  utils::write.table(out, file = export.file.name, sep = "\t", col.names=FALSE, row.names = FALSE)
}

# export.sirius(pc.bio.object = pc.bio.sub, export.file.name = "C:/Temp/20250703/pc.bio.sirius.tsv")
