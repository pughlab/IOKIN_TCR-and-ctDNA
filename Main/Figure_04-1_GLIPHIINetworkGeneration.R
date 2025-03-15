library(tidyverse)
library(igraph)
library(graphlayouts)

### NOTE: This script contains the graph-based GLIPHII analysis pipeline that generates three outputs:

### A) 09_IOKIN_GLIPHIINetwork_PreHLAConvergence.rds
###        An igraph network representing TRB-CDR3 sequences (nodes) 
###        connected by GLIPH-II identified global or local specificity 
###        signatures (edges, indicated by the 'type' field in the GLIPH-II output). 
###        This network is the result of graph-based filtering, including criteria 
###        such as community size and node connectivity determined by clique formation. 
###        HLA convergence criteria have not been applied.

### B) 10_IOKIN_GLIPHIICommunityFeatures_PreHLAConvergence.csv
###        A dataframe summarizing GLIPHII specificity community features.

### C) 11_IOKIN_GLIPHIISuperNodes_PreHLAConvergence.rds
###        This network is a simplified version of network A, designed for 
###        improved visualization and interpretability.  
###        Individual TRB-CDR3 nodes from A that belong to the same GLIPHII-defined 
###        specificity community are collapsed into a single supernode.
###        Therefore, each supernode represents a unique specificity community, rather than an individual TRB-CDR3 sequence.  
###        Due to the large size of the original network A, this supernode representation 
###        provides a more manageable and visually accessible overview of the overall landscape of specificity communities.


data_path <- "https://raw.githubusercontent.com/pughlab/IOKIN_TCR-and-ctDNA/refs/heads/main/Data"

### ---------------------------------------------------------------------------------------------------------
### Defining the colors for external nodes for visualization purpose:
### ---------------------------------------------------------------------------------------------------------

color_shape_pallette <- data.frame(
        Source = c(
                "HNSCC" , "HomoSapiens" ,
                "EBV" ,  "CMV" ,
                "MCPyV" , "HCV" , "Influenza" ,
                "HPV" ,  "M.tuberculosis" ,
                "CEF"  ,   "S-pneumoniae" ,
                "YFV" , "HTLV" ,
                "IOKINBlood") ,

        colors = c(
                "#CDEDFA" , #"HNSCC" ,
                "#FF9F01" , #"HomoSapiens" ,
                "#AC261B" , #"EBV" ,
                "#266D1F" , #"CMV" ,
                "#BAC70D" , #"MCPyV" ,
                "#FFDF01" , #"HCV" ,
                "#0C99BA" , #"Influenza" ,
                "#1A578F" , #"HPV" ,
                "#FFC86D" , #"M.tuberculosis" ,
                "#C1A132" , #"CEF"  ,
                "#D07F89" , #"S-pneumoniae"  ,
                "#F3A3B5" , #"YFV" ,
                "#6BB8B3" , #"HTLV"
                "#A0001D"   #"IOKINBlood"
        ) )

### ---------------------------------------------------------------------------------------------------------
### Clinical data:
### ---------------------------------------------------------------------------------------------------------

ClinicalData_fname <- "02_IOKIN_ClinicalData.csv"
ClinicalData <- read_csv(file.path(data_path , ClinicalData_fname)) %>%
        dplyr::select(Patient_id , BestResponse , CancerType , 
                      `LA/Metastastic` , `HPV status`) %>%
        rename (`LocalAdvanced-orMetastatic` = `LA/Metastastic`,
                HPVStatus = `HPV status`)

### ---------------------------------------------------------------------------------------------------------
### Reading the GLIPHII output data for network generation: 
### ---------------------------------------------------------------------------------------------------------

GliphOutput_fname <- "07_IOKIN_BuffyCoat_GLIPHIIOutput.csv"
gliph_data <- readr::read_csv(file.path(data_path , GliphOutput_fname ))%>%
        select(
                type , TcRb , V ,
                pattern , Sample ) %>%
        mutate(
                Source = unlist(
                        sapply(
                                strsplit(Sample, split = ":", fixed = TRUE),
                                function(x) x[[1]][1], simplify=FALSE)) ,
                
                Patient_id = ifelse(grepl ("IOKINBlood" , Sample ) ,
                                    str_extract(Sample, "(?<=:).*?(?=_)"),
                                    "TBD"),
                
                Cycle = ifelse(grepl ("IOKINBlood" , Sample ) ,
                               unlist (sapply(strsplit(Sample, split = "_", fixed = TRUE),
                                              function(x) x[2][1],
                                              simplify=FALSE)),
                               "TBD"),
                
                Node = stringr::str_c(TcRb , "_" , Sample)) %>%
        
        left_join(ClinicalData , by = "Patient_id") %>%
        
        mutate_at(vars(BestResponse ,
                       CancerType ,
                       `LocalAdvanced-orMetastatic` ,
                       HPVStatus), ~ ifelse(is.na(.), "TBD", .)) %>%
                  
        left_join (color_shape_pallette , by = "Source")
                  
### ---------------------------------------------------------------------------------------------------------
### Raw network generation: 
### ---------------------------------------------------------------------------------------------------------

# Edges are defined as GLIPH-II identified global and local specificity signatures (denoted as ‘type’ in GLIPHII outputs):

specificity_signatures <- unique(gliph_data$type)
n_count <- 1
for (sig in specificity_signatures) {

        cdr3s <- gliph_data$Node [gliph_data$type == sig]
        if (length(cdr3s) > 1) {

                if (n_count == 1) {

                        EDGES <- data.frame (t (combn(cdr3s , 2 , simplify = TRUE)))

                        n_count <- n_count + 1
                }
                else {

                        edges <- data.frame (t (combn(cdr3s , 2 , simplify = TRUE)))
                        EDGES <- bind_rows(EDGES , edges)
                }

        }

}

network <- graph_from_data_frame(d = EDGES ,
                                 vertices = unique(gliph_data %>%
                                                           group_by(Node) %>%
                                                           mutate(type = paste(sort (unique (type) ) , collapse = ",")) %>%
                                                           unique() %>%
                                                           dplyr::select(Node ,
                                                                         Source ,
                                                                         Patient_id ,
                                                                         Cycle ,
                                                                         CancerType ,
                                                                         BestResponse ,
                                                                         HPVStatus ,
                                                                         `LocalAdvanced-orMetastatic` ,
                                                                         TcRb ,
                                                                         V ,
                                                                         type ,
                                                                         colors)) ,
                                 directed = FALSE)



#Removing multi-edges
network <- simplify(network)

### ---------------------------------------------------------------------------------------------------------
### Network trimming via clique identification: 
### ---------------------------------------------------------------------------------------------------------

### After constructing the networks, filtrations are required to withdraw loosely connected nodes 
### to increase the precision of downstream specificity annotations.
### A clique is a complete subgraph within a larger graph, 
### where every node is directly connected to every other node in that subset. 
### All nodes that do not form a clique with size four are excluded from further analysis.

MaxClique_collection <- max_cliques(network,
                                    min = 4,
                                    max = NULL)


network <- induced_subgraph(network,
                            unlist(MaxClique_collection))

### ---------------------------------------------------------------------------------------------------------
### Community detection for specificity annotation: 
### ---------------------------------------------------------------------------------------------------------
                  
communities_leiden <- igraph::cluster_leiden(network,
                                             objective_function = "CPM",
                                             resolution_parameter = 0.2 ,
                                             n_iterations = 10000,
                                             beta = 0.01)


network <- set_vertex_attr(network,
                           name = "LeidenCommunity",
                           value = communities_leiden$membership)   

### ---------------------------------------------------------------------------------------------------------
### Fine tuning the communities: 
### ---------------------------------------------------------------------------------------------------------

largest_cliques_list <- c()


for (community_id in unique(communities_leiden$membership)) {
        subgraph <- igraph::subgraph(
                network ,
                which(V(network)$LeidenCommunity == community_id))

        if (clique_num (subgraph) > 4) {

                largest_cliques_list <- append(largest_cliques_list ,
                                               unique (names(unlist(largest_cliques(subgraph)))) )
        }
        rm (subgraph)
}


trimmedNetwork <- delete_vertices(network,
                                  setdiff(V(network)$name , largest_cliques_list))

### Saving the network:                        
saveRDS(trimmedNetwork, "09_IOKIN_GLIPHIINetwork_PreHLAConvergence.rds" )

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
write_csv(GLIPHII_Community_stats , "10_IOKIN_GLIPHIICommunityFeatures_PreHLAConvergence.csv" )


### ---------------------------------------------------------------------------------------------------------
### Loop through communities and check for HLA convergence:
### ---------------------------------------------------------------------------------------------------------

HLA_data <- read_csv(file.path(data_path , "08_IOKIN_HLAGenotyping.csv"))

AbstractAnnotationOfInterest <- "" # any of c("Bacterial" , "Cross-species" , "HNSCC" , 
                                   #            "Human tumour and auto-antigens" , 
                                   #            "IOKIN_intrinsic" , "Viral")
                               
CommunityCollection <- (GLIPHII_Community_stats %>%
                                filter(Abstract_Annotation == AbstractAnnotationOfInterest ) %>%
                                arrange(desc(Community_size)))$Community_id

for (i in CommunityCollection ) {
        
        print(i)
        
        nt <- igraph::subgraph(
                trimmedNetwork ,
                which(V(trimmedNetwork)$LeidenCommunity ==  i ))
        print(cat(V(nt)$name [order(V(nt)$Source)], sep = "\n"))
        print(sort(unique (V(nt)$Patient_id [V(nt)$Source == "IOKINBlood"])) )
        
        # print(HLA_data %>%
        #               arrange(Patient_id) %>%
        #               filter(Patient_id %in% c(sort ( unique (V(nt)$Patient_id [V(nt)$Source == "IOKINBlood"]))))%>%
        #               filter( grepl("HLA-A|HLA-B|HLA-C" , Locus1) | grepl("HLA-A|HLA-B|HLA-C" , Locus2)))
        
        
        # print(HLA_data %>%
        #               arrange(Patient_id) %>%
        #               filter(Patient_id %in% c(sort ( unique (V(nt)$Patient_id [V(nt)$Source == "IOKINBlood"]))))%>%
        #               filter( grepl("HLA-A|HLA-B|HLA-C|HLA-DRB1|HLA-DQA1|HLA-DQB1|HLA-DPA1|HLA-DPB1" , Locus1) | 
        #                               grepl("HLA-A|HLA-B|HLA-C|HLA-DRB1|HLA-DQA1|HLA-DQB1|HLA-DPA1|HLA-DPB1"  , Locus2)) %>%
        #               arrange(Locus1 , Locus2 , Patient_id),
        #       n = 1000)
        
        
        
        print("####################################################")
        print("#################################")
        print("##################")
        
}                               

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
        ###The mapping argument passed to contract function is a numeric vector that specifies the mapping. 
        ### We pass the community ids as the mapping argument. 
        ### However, because we've trimmed some communities, some community ids are missing.
        ### For those withdrawn communities, empty super nodes are generated, which are not favorable. 
        ### Therefore, we need to provide mock community ids to prevent empty supernode generation. 
        ### We add a new column to GLIPHII_Community_stats as Mock_Community_id to do this.

                               
### Adding attributes to each super-node from GLIPHII_Community_stats dataframe:
Abstract_Network <- set_vertex_attr(Abstract_Network,
                                    name = "Community_size",
                                    value = GLIPHII_Community_stats$Community_size)

Abstract_Network <- set_vertex_attr(Abstract_Network,
                                    name = "Component_id",
                                    value = GLIPHII_Community_stats$Component_id)

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
### Super-node network visualization:
### ---------------------------------------------------------------------------------------------------------

### Defining the coordination of each super-node using layouts readily available in both igraph and 
### graphlayouts packages. This step is up to the user and based on user's visualization preferences,
### Any other packages or algorithms can be utilized:
InitialLayout <- layout_with_graphopt(Abstract_Network)


### For extra, minimal manual tweaking, we can use tkid:
tkid <- tkplot(
        Abstract_Network,
        vertex.size =  V (Abstract_Network)$Community_size**0.85 ,
        vertex.color = V (Abstract_Network)$Color ,
        vertex.label = NA ,
        vertex.label.color = "#000000" ,
        vertex.label.cex = 0.5 ,
        vertex.label.degree = 75 ,
        vertex.label.family = "Helvetica" ,
        vertex.shape = "circle" ,
        edge.width = 0.5 ,
        edge.color = "#000000" ,
        rescale = FALSE ,
        ###To get all the available layouts: grep("^layout_", ls("package:igraph"), value=TRUE)[-1]
        layout = norm_coords(InitialLayout ,
                ymin=-0.9,
                ymax=0.9,
                xmin=-0.9,
                xmax=0.9))


### Get the modified layout:
layout <- tk_coords(tkid)


### We can save the layout for later use or add the X and Y coordinations 
### as node attributes to the super-node network:

Abstract_Network <- set_vertex_attr(Abstract_Network,
                                    name = "LayoutX",
                                    value = layout [,1])

Abstract_Network <- set_vertex_attr(Abstract_Network,
                                    name = "LayoutY",
                                    value = layout [,2])


### saveRDS(Abstract_Network,"11_IOKIN_GLIPHIISuperNodes_PreHLAConvergence.rds" )
