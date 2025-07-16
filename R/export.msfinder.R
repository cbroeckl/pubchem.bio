#' export.msfinder
#'
#' export pubchem.bio pc.bio syle data.table to format suitable for MSFinder input.  
#' @details takes output from 'build.pubchem.bio' or 'build.taxon.metabolome' functions, reformatting, and exporting to input format suitable for MSFinder.  
#' @param pc.bio.object input data.table, generated from 'build.pubchem.bio' or 'build.taxon.metabolome' functions
#' @param export.file.name valid file path and name.  Extension should be listed as '.tsv'.
#' @return nothing - file written to disk. 
#' @author Corey Broeckling
#' 
#' @export 
#' 

export.msfinder <- function(
    pc.bio.object = NULL,
    export.file.name = NULL
) {
  
  out <- data.frame(
    'Title' = pc.bio.object$name,
    'InChIKey' = pc.bio.object$inchikey,
    'Short InChIKey' = pc.bio.object$inchikey.first.block,
    'PubChem CID' = pc.bio.object$cid,
    'Exact mass' = pc.bio.object$monoisotopic.mass,
    'Formula' = pc.bio.object$formula,
    'SMILES' = pc.bio.object$smiles,
    'Database ID' = paste0(pc.bio.object$cid)
  )
  utils::write.table(out, file = export.file.name, sep = "\t", col.names=TRUE, row.names = FALSE)
}

# export.msfinder(pc.bio.object = pc.bio.sub, export.file.name = "C:/Temp/20250703/pc.bio.msfinder.tsv")
