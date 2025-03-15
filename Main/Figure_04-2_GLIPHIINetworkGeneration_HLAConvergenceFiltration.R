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


### ---------------------------------------------------------------------------------------------------------
### Storing community specificities and their details in a dataframe for downstream analysis:
### ---------------------------------------------------------------------------------------------------------

# Helper function to collapse unique, sorted values into a comma-separated string
collapse_unique <- function(x) {
        paste(sort(unique(x)), collapse = ",")
}

vertex_data <- tibble(
        name = V(trimmedNetwork)$name,
        LeidenCommunity = V(trimmedNetwork)$LeidenCommunity,
        Component = components(trimmedNetwork)$membership,
        Source = V(trimmedNetwork)$Source,
        Patient_id = V(trimmedNetwork)$Patient_id,
        CancerType = V(trimmedNetwork)$CancerType,
        BestResponse = V(trimmedNetwork)$BestResponse,
        HPVStatus = V(trimmedNetwork)$HPVStatus)




GLIPHII_Community_stats <- vertex_data %>%
        group_by(LeidenCommunity) %>%
        summarize(
                Community_size = n(), 
                
                PatientDerived_Nodes = collapse_unique(name[Source == "IOKINBlood"]),
                NumberOf_PatientDerived_Nodes = sum(Source == "IOKINBlood"),
                Number_of_Patients = n_distinct(Patient_id[Source == "IOKINBlood"]),
                Patient_id = collapse_unique(Patient_id[Source == "IOKINBlood"]),
                
                CancerType_Range = collapse_unique(CancerType[Source == "IOKINBlood"]),
                BestResponse_Range = collapse_unique(BestResponse[Source == "IOKINBlood"]),
                HPVStatus_Range = collapse_unique(HPVStatus[Source == "IOKINBlood"]),
                
                Community_specificity = ifelse(
                        any(Source != "IOKINBlood"), 
                        collapse_unique(Source[Source != "IOKINBlood"]),  
                        "IOKIN_intrinsic"  
                ),
                
                KnownExternalTCRs = ifelse(collapse_unique(name[Source != "IOKINBlood" & Source != "HNSCC"]) == "", 
                                           "TBD" , 
                                           collapse_unique(name[Source != "IOKINBlood" & Source != "HNSCC"])),
                
                .groups = "drop") %>%
        rename(Community_id = LeidenCommunity)  


                        
rm (vertex_data)

### ---------------------------------------------------------------------------------------------------------
### Defining colors and abstract annotations:
### ---------------------------------------------------------------------------------------------------------
### Specificity annotation maps and their colors:                        
Annotation_Ref <- tibble(
        ExternalSpecificity = c("IOKIN_intrinsic" ,

                                "HNSCC" ,
                                "HomoSapiens" ,

                                "CEF" ,  "CMV" ,
                                "EBV" ,  "HCV" ,
                                "HPV" ,  "MCPyV" ,
                                "Influenza" , "DENV" ,
                                "HTLV-1" , "YFV" ,

                                "S-pneumoniae" , "M.tuberculosis" ) ,

        Abstract_Annotation = c("IOKIN_intrinsic" ,

                                "HNSCC" ,
                                "Human tumour and auto-antigens" ,

                                "Viral" ,  "Viral" ,
                                "Viral" ,  "Viral" ,
                                "Viral" ,  "Viral" ,
                                "Viral" ,  "Viral" ,
                                "Viral" ,  "Viral" ,

                                "Bacterial" , "Bacterial" ) )


Color_Ref <- tibble(
        Abstract_Annotation = c("IOKIN_intrinsic" ,
                                "HNSCC" ,
                                "Human tumour and auto-antigens" ,
                                "Viral" ,
                                "Bacterial" ,
                                "Cross-species") ,
        Color = c(  "#000000" ,
                    "#CDEDFA" ,
                    "#FF9F01" ,
                    "#AC261B" ,
                    "#BAC70D" ,
                    "#FFDF01"))
                        
### ---------------------------------------------------------------------------------------------------------
### Functions for Annotating Detailed and Abstract Specificities:
### ---------------------------------------------------------------------------------------------------------

# Function to extract Human genes
getGenes <- function(string) {
        fullStrings <- grep("HomoSapiens" ,
                            unlist (strsplit(string , ",")),
                            value = TRUE)
        genes <- paste (unique(sapply(
                strsplit (fullStrings , "__"),
                function(x) x[[3]])),
                collapse = ",")

        return(genes)
}




# Function for Abstract Annotation:
AbstractAnnotate <- function(string) {

        if (string %in% c("IOKIN_intrinsic" )) {
                annotation <- Annotation_Ref$Abstract_Annotation [Annotation_Ref$ExternalSpecificity == string]

        }
        else {

                if (length(unlist (str_split(string , ","))) == 1) {
                        annotation <- Annotation_Ref$Abstract_Annotation [Annotation_Ref$ExternalSpecificity == string]
                }

                else if (length(unlist (str_split(string , ","))) == 2 &
                         grepl("HNSCC" , string) ) {

                        annotation <- Annotation_Ref$Abstract_Annotation [Annotation_Ref$ExternalSpecificity == grep ("HNSCC" ,
                                                                                                                      unlist(str_split(string , ",")) ,
                                                                                                                      value = TRUE ,
                                                                                                                      invert = TRUE)]
                }

                else {
                        annotation <- "Cross-species"
                }

        }



        return(annotation)
}



# Function for Detailed Annotation:
DetailAnnotate <- function(string) {


        if (string %in% c("IOKIN_intrinsic" )) {
                annotation <- string

        }
        else {
                if (length(unlist (str_split(string , ","))) == 1) {
                        annotation <- string
                }
                else {
                        annotation <- paste( grep( "HNSCC",
                                                   unlist (str_split(string , ",")) ,
                                                   value = TRUE ,
                                                   invert = TRUE) , collapse = ",")

                }
        }


        return(annotation)
}

### ---------------------------------------------------------------------------------------------------------
### Adding Abstract and detailed annotation data to the dataframe of community specificities:
### ---------------------------------------------------------------------------------------------------------

GLIPHII_Community_stats <- GLIPHII_Community_stats %>%
        rowwise() %>%
        mutate(Abstract_Annotation = AbstractAnnotate (Community_specificity) ,
               Detailed_Annotation = DetailAnnotate (Community_specificity)) %>%
        mutate(Detailed_Annotation = str_replace(Detailed_Annotation , "HomoSapiens" , getGenes (KnownExternalTCRs))) %>%
        left_join(Color_Ref ,
                  by = "Abstract_Annotation")

                               
### Saving the dataframe:
write_csv(GLIPHII_Community_stats , "14_IOKIN_GLIPHIICommunityFeatures_PostHLAConvergence.csv" )


### ---------------------------------------------------------------------------------------------------------
### Condensing the nodes participating in each community into super-nodes:
### ---------------------------------------------------------------------------------------------------------

GLIPHII_Community_stats$Mock_Community_id <- 1:nrow(GLIPHII_Community_stats)

trimmedNetwork <- set_vertex_attr(trimmedNetwork,
                           name = "Mock_Community_id",
                           value = GLIPHII_Community_stats$Mock_Community_id [match(V(trimmedNetwork)$LeidenCommunity ,
                                                                                    GLIPHII_Community_stats$Community_id) ])

Abstract_Network <- igraph::contract(
        graph = trimmedNetwork ,
        mapping = V(trimmedNetwork)$Mock_Community_id)

                               
### Adding attributes to each super-node from GLIPHII_Community_stats dataframe:

Abstract_Network <- set_vertex_attr(Abstract_Network,
                                    name = "Community_id",
                                    value = GLIPHII_Community_stats$Community_id)


Abstract_Network <- set_vertex_attr(Abstract_Network,
                                    name = "Community_specificity",
                                    value = GLIPHII_Community_stats$Community_specificity)

Abstract_Network <- set_vertex_attr(Abstract_Network,
                                    name = "Detailed_Annotation",
                                    value = GLIPHII_Community_stats$Detailed_Annotation)

Abstract_Network <- set_vertex_attr(Abstract_Network,
                                    name = "Abstract_Annotation",
                                    value = GLIPHII_Community_stats$Abstract_Annotation)

Abstract_Network <- set_vertex_attr(Abstract_Network,
                                    name = "Color",
                                    value = GLIPHII_Community_stats$Color)

Abstract_Network <- set_vertex_attr(Abstract_Network,
                                    name = "HPVStatus",
                                    value = GLIPHII_Community_stats$HPVStatus_Range)

Abstract_Network <- set_vertex_attr(Abstract_Network,
                                    name = "CancerType",
                                    value = GLIPHII_Community_stats$CancerType_Range)

Abstract_Network <- set_vertex_attr(Abstract_Network,
                                    name = "Number_of_Patients",
                                    value = GLIPHII_Community_stats$Number_of_Patients)

Abstract_Network <- set_vertex_attr(Abstract_Network,
                                    name = "BestResponse",
                                    value = GLIPHII_Community_stats$BestResponse_Range)


Abstract_Network <- set_vertex_attr(Abstract_Network,
                                    name = "Patient_id",
                                    value = GLIPHII_Community_stats$Patient_id)



Abstract_Network <- simplify(Abstract_Network)


### ---------------------------------------------------------------------------------------------------------
### Extracting the layout for visualization of the super-node network from pre-HLA convergence super node network 
### for consistency of node coordinates to enable pre , post HLA convergence filtration comparison :
### ---------------------------------------------------------------------------------------------------------


OldAbstract_Network <- read_rds(file.path(data_path , "11_IOKIN_GLIPHIISuperNodes_PreHLAConvergence.rds"))

layout <- as.matrix(data.frame(V(OldAbstract_Network)$LayoutX , V(OldAbstract_Network)$LayoutY , V(OldAbstract_Network)$Community_id) %>%
                            filter(`V.OldAbstract_Network..Community_id` %in% GLIPHII_Community_stats$Community_id) %>%
                            dplyr::select(-`V.OldAbstract_Network..Community_id`))




Abstract_Network <- set_vertex_attr(Abstract_Network,
                                    name = "LayoutX",
                                    value = layout [,1])

Abstract_Network <- set_vertex_attr(Abstract_Network,
                                    name = "LayoutY",
                                    value = layout [,2])

saveRDS(Abstract_Network, "15_IOKIN_GLIPHIISuperNodes_PostHLAConvergence.rds")
    
