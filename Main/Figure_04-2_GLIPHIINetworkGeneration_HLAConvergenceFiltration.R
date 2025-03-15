library(tidyverse)
library(igraph)
library(graphlayouts)

### NOTE: This script contains the graph-based GLIPHII analysis pipeline that generates three outputs:

### A) 13_IOKIN_GLIPHIINetwork_PostHLAConvergence.rds
###        HLA-filtered version of 09_IOKIN_GLIPHIINetwork_PreHLAConvergence.rds


### B) 14_IOKIN_GLIPHIICommunityFeatures_PostHLAConvergence.csv
###        HLA-filtered version of 10_IOKIN_GLIPHIICommunityFeatures_PreHLAConvergence.csv

### C) 15_IOKIN_GLIPHIISuperNodes_PostHLAConvergence.rds
###        HLA-filtered version of 11_IOKIN_GLIPHIISuperNodes_PreHLAConvergence.rds


### How HLA convergence filtration works:
### A de-orphanized community, where a TRB-CDR3 with an experimentally validated cognate pMHC is present, 
### is classified as HLA-convergent if all patients contributing TRB-CDR3s to the community share the specific 
### HLA allele presenting that pMHC. Nodes that belong to patients without a convergent HLA are removed.
### The shared HLA allele supports the patients' ability to present the same target antigen. 
### For orphan communities (including both IO-KIN-intrinsic and HNSCC-convergent), HLA convergence indicates 
### that all contributing patients share at least one common HLA allele, suggesting the possibility of a 
### shared/homologous pMHC, although the specific target antigen remains unknown. Communities without a minimum
### of one shared HLA allele among all participating patients are dropped.



data_path <- "https://raw.githubusercontent.com/pughlab/IOKIN_TCR-and-ctDNA/refs/heads/main/Data"

### The original GLIPHII network:
trimmedNetwork <- read_rds(file.path(data_path, "09_IOKIN_GLIPHIINetwork_PreHLAConvergence.rds") )

### The original CLIPHII community features dataframe:
GLIPHII_Community_stats <- read_csv( file.path(data_path , "10_IOKIN_GLIPHIICommunityFeatures_PreHLAConvergence.csv" )) 

### ---------------------------------------------------------------------------------------------------------
### Removing nodes and communities with HLA inconsistency:
### ---------------------------------------------------------------------------------------------------------

### Vector of nodes that need to be dropped due to HLA inconsistency:
RemovalNodeList <- read_lines(file.path(data_path , "12_IOKIN_List-of-Nodes-to-Remove_HLAInconsistent.txt") ) %>%
        {\(x) x[x != ""]}()

### Vector of community IDs that need to be dropped due to HLA inconsistency:
RemovalCommunityList <- (GLIPHII_Community_stats %>% 
        filter(HLA_Status == "Divergent"))$Community_id


### Removing the HLA-inconsistent nodes:
trimmedNetwork <- igraph::delete_vertices(trimmedNetwork ,
                                          V(trimmedNetwork)$name [V(trimmedNetwork)$LeidenCommunity %in% RemovalCommunityList])

### Removing the HLA-inconsistent communities:
trimmedNetwork <- igraph::delete_vertices(trimmedNetwork , RemovalNodeList )

### Saving the HLA-filtered GLIPHII network:
saveRDS(trimmedNetwork, "13_IOKIN_GLIPHIINetwork_PostHLAConvergence.rds")
