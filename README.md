# C.O.M.I.C.S. A pipeline for the composite identification of selection across multiple genomic scans using Invariant Coordinate Selection in R
(a friendly shiny-app)

Identifying loci that are under selection versus those that are evolving neutrally is a common challenge in evolutionary genetics. Moreover, with the increase in sequence data, genomic studies have begun to incorporate the use of multiple methods to identify candidate loci under selection. Composite methods are usually implemented to transform the data into a multi-dimensional scatter where outliers are identified using a distance metric, the most common being Mahalanobis distance. However, studies have shown that the power of Mahalanobis distance reduces as the number of dimensions or selection tests increases. 
  We implement a new approach: Calling Outlier loci from Multi-dimensional data using Invariant Coordinate Selection, a pipeline for identifying outlier loci from multiple selection scans where the signal is heterogeneous and different selection algorithms capture different signatures of the selective process. The use of ICS provides a better balance between precision and recall than other commonly available methods like the use of Mahalanobis distances.

![Sensitivity_Combined](https://github.com/user-attachments/assets/eb691a57-53e4-4f84-bbda-f8e5f843baed)


# Using COMICS
## Installing Dependencies

COMICS relies on several dependencies:
It is important to note that COMICS can only be applied to R Studios and not the base R console. 
For R Studios visit: https://rstudio.com/products/rstudio/download/


`install.packages(c("shiny", "ggplot2", "tidyr", "bslib", "mvtnorm", "ICS", "moments", "ICSOutlier"), dependencies = TRUE)`

Users may also need to install additional tools for a more efficient ues of COMICS. See:  https://www.rstudio.com/products/rpackages/devtools/

## Installing COMICS

There are two R-scripts that make up COMICS. The first script contains code that designs the user interface for COMICS, otherwise known as the fluid page. -> COMICS.design.R. The second script includes the input of the data and invokes the ICS and figure generation.

For COMICS istallation after dependencies are downloaded:

`install.packages("devtools", dependencies = TRUE)`
 
 `library(devtools)`
 
`devtools::install_github("oeco28/COMICS", build_vignettes = FALSE)`
 
 `library(COMICS)`
 
`run_comics()`

If the above installation does not work, the base COMICS scripts can be downloaded from /scripts from this repository.
   
## Input Files

COMICS incorperates the use of two different input files, both of which are necessary for the complete use of the package. See the example input files for more information about the format of the two files.

input file 1 -> Scan_data.txt -> the main inout file that contains the test statistics for each selection scan for each position on a given chromosome/genome.

input file 2 -> Genome_configuration.txt -> a data table that contains the number of chromosomes and their putativley lengths.

### Additional Questions

If there are any issues that arise, please open a new issue.

# Citing COMICS

The method was succesfully implemented in: Nelson JT, Motamayor JC, Cornejo OE. [Environment and pathogens shape local and regional adaptations to climate change in the chocolate tree, Theobroma cacao L](https://doi.org/10.1111/mec.15754). Mol Ecol. 2021 Feb;30(3):656-669. doi: 10.1111/mec.15754. Epub 2020 Dec 26. PMID: 33247971. 

The publication describing the app will be made available sometime soon: Cornejo, OE and Nelson, JT. COMICS: A pipeline for the composite identification of selection across multiple genomic scans using Invariant Coordinate Selection in R.


