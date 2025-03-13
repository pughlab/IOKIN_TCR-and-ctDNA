library(tidyverse)
library(igraph)
library(graphlayouts)


### Defining the colors for external nodes for visualization purpose:
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

data_path <- "https://raw.githubusercontent.com/pughlab/IOKIN_TCR-and-ctDNA/refs/heads/main/Data"

### Clinical data:
ClinicalData_fname <- "02_IOKIN_ClinicalData.csv"
ClinicalData <- read_csv(file.path(data_path , ClinicalData_fname)) %>%
        dplyr::select(Patient_id , BestResponse , CancerType , 
                      `LA/Metastastic` , `HPV status`) %>%
        rename (`LocalAdvanced-orMetastatic` = `LA/Metastastic`,
                HPVStatus = `HPV status`)

### Reading the GLIPHII output data for network generation: 
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
### Raw network generation: ---------------------------------------------------------------------------------
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
### Network trimming via clique identification: -------------------------------------------------------------
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
### Community detection for specificity annotation: ---------------------------------------------------------
                  
communities_leiden <- igraph::cluster_leiden(network,
                                             objective_function = "CPM",
                                             resolution_parameter = 0.2 ,
                                             n_iterations = 10000,
                                             beta = 0.01)


network <- set_vertex_attr(network,
                           name = "LeidenCommunity",
                           value = communities_leiden$membership)   

### ---------------------------------------------------------------------------------------------------------
### Fine tuning the communities: ----------------------------------------------------------------------------

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


