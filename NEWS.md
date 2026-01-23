# pubchem.bio 1.1.0
* add browse.pubchem.bio Rshiny dashboard

# pubchem.bio 1.0.5
* fix export.ComboundDb functionality bugs
* bring CIDs in from PubChemLite as an added data source to expand exogenous compound coverage.

# pubchem.bio 1.0.4
* enable reporting of multiple LCA's in pubchem.bio tables
* add reporting of LCA name(s) in pubchem.bio tables
* add option in build.pubchem.bio to use the neutral child rather than the parent, even if 'use.parent' option is TRUE.
* updated cid.lca.rda data to include lca.level column, used moving forward and reported in pubchem.bio object output.

# pubchem.bio 1.0.3
* added export.CompoundDb function to enable 'conversion to a CompoundDb' sqlite database and compatibility with the RForMassSpectrometry metapackage
* added export.CompoundDb demonstration to vignette

# pubchem.bio 1.0.2
* added build.element.count function
* added filtering demonstration to vignette

# pubchem.bio 1.0.1
* corrected data.table namespace problem
* improve robustness of pathway data download

# pubchem.bio 1.0.0
* Initial CRAN submission.
