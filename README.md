Metabolomics approaches utilize analytical instrumentation to detect and quantify small molecules.  The chemical space covered by metabolomics is potentially massive, with more than 100 million structures in Pubchem.  Efforts to curate custom metabolome databases, such as the [Human Metabolome Database](https://www.hmdb.ca/) are incredibly valuable, enabling annotation efforts to focus on only those chemicals most likely to be found in the samples being analyzed.  Metabolomics, however, is well suited the analysis of metabolomes from _any_ species, and there are very few taxonomically informed databases available.  

Pubchem, which is a massive repository of chemical structures, is also extremely rich in chemical metadata.  Several biologicallly-focused metabolomic databases have deposited their structures data into PubChem, thousands of metabolic pathways are present, and there are rich taxonomic-structure pairs through taxonomy-informed natural products database efforts.  

The pubchem.bio package enables researchers to 

* Download metabolomics-centric subset of PubChem onto their local computer
* Build a metabolomic structured database (data.table) of biological compounds from PubChem
* Define a core biological 'primary' metabolome, comprising metabolites plausibly found in any species
* Develop custom metabolomic structure databases using selected or all available taxonomic data in PubChem, and scoring metabolites based on taxonomic proximity.

As such, this package will help to streamline metabolomic database creation, increase the accuracy of annotation efforts by reducing the chemical search space to the most likely encounted metabolites, in a manner informed by a wealth of biological knowledge. 
