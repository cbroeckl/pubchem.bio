% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build.taxon.metabolome.R
\name{build.taxon.metabolome}
\alias{build.taxon.metabolome}
\title{build.taxon.metabolome}
\usage{
build.taxon.metabolome(
  pc.directory = NULL,
  taxid = c(),
  get.properties = TRUE,
  full.scored = TRUE,
  aggregation.function = max,
  threads = 8,
  db.name = "custom.metabolome",
  rcdk.desc = c("org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.AcidicGroupCountDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.BasicGroupCountDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor"),
  pubchem.bio.object = NULL,
  cid.lca.object = NULL,
  taxid.hierarchy.object = NULL,
  output.directory = NULL
)
}
\arguments{
\item{pc.directory}{directory from which to load pubchem .Rdata files}

\item{taxid}{integer vector of integer NCBI taxonomy IDs.  i.e.  c(9606, 1425170 ) for Homo sapiens and Homo heidelbergensis.}

\item{get.properties}{logical. if TRUE, will return rcdk calculated properties:  XLogP, TPSA, HBondDonorCount and HBondAcceptorCount.}

\item{full.scored}{logincal.  default = FALSE.  When false, only metabolites which map to the taxid(s) are returned.  When TRUE, all metabolites are returned, with scores assigned based on the distance of non-mapped metabolites to the root node.  i.e. specialized metabolites from distantly related species are going to be scored at or near zero, specialized metabolites of mores similar species higher, and more conserved metabolites will score higher than ore specialized.}

\item{aggregation.function}{function. default = max.  can use mean, median, min, etc, or a custom function.  Defines how the aggregate score will be calculated when multiple taxids are used.}

\item{threads}{integer. how many threads to use when calculating rcdk properties.  parallel processing via DoParallel and foreach packages.}

\item{db.name}{character. what do you wish the file name for the saved version of this database to be?  default = 'custom.metabolome', but could be 'taxid.4071' or 'Streptomyces', etc.  Saved as an .Rdata file in the 'pc.directory' location.}

\item{rcdk.desc}{vector. character vector of valid rcdk descriptors.  default = rcdk.desc <- c("org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor", "org.openscience.cdk.qsar.descriptors.molecular.AcidicGroupCountDescriptor", "org.openscience.cdk.qsar.descriptors.molecular.BasicGroupCountDescriptor", "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor"). To see descriptor categories: 'dc <- rcdk::get.desc.categories(); dc' .  To see the descriptors within one category: 'dn <- rcdk::get.desc.names(dc[4]); dn'. Note that the four default parameters are relatively fast to calculate - some descriptors take a very long time to calculate.  you can calculate as many as you wish, but processing time will increase the more descriptors are added.}

\item{pubchem.bio.object}{R data.table, generally produced by build.pubchem.bio; alternatively, define pc.directory}

\item{cid.lca.object}{R data.table, generally produced by build.cid.lca; alternatively, define pc.directory}

\item{taxid.hierarchy.object}{R data.table, generally produced by get.pubchem.ftp; alternatively, define pc.directory}

\item{output.directory}{directory to which the pubchem.bio database is saved.  If NULL, will try to save in pc.directory (if provided), else not saved.}
}
\value{
a data frame containing pubchem CID ('cid'), and lowest common ancestor ('lca') NCBI taxonomy ID integer. will also save to pc.directory as .Rdata file.
}
\description{
utilizes downloaded and properly formatted local pubchem data created by 'get.pubchem.ftp' function to filter a dataset created by 'build.pubchem.bio' function
}
\details{
utilizes downloaded and properly formatted local pubchem data created by 'get.pubchem.ftp' function
}
\examples{
data('cid.lca', package = "pubchem.bio")
data('pubchem.bio', package = "pubchem.bio")
data('taxid.hierarchy', package = "pubchem.bio")
my.taxon.db <- build.taxon.metabolome(pc.bio.out = pubchem.bio, cid.lca.object = cid.lca, 
taxid.hierarchy.object = taxid.hierarchy, get.properties = FALSE, threads = 1)
head(my.taxon.db)

}
\author{
Corey Broeckling
}
