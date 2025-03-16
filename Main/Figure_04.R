library(tidyverse)
library(igraph)

data_path <- "https://raw.githubusercontent.com/pughlab/IOKIN_TCR-and-ctDNA/refs/heads/main/Data"

Abstract_Network <- read_rds(file.path(data_path , "11_IOKIN_GLIPHIISuperNodes_PreHLAConvergence.rds"))
network <- readr::read_rds(file.path(data_path , "13_IOKIN_GLIPHIINetwork_PostHLAConvergence.rds" ))
GLIPHII_Community_stats <- read_csv(file.path(data_path , "14_IOKIN_GLIPHIICommunityFeatures_PostHLAConvergence.csv" ))


### Visualization -------------------------------------------------------------------------------------------------

# Fig 4A

plot(
        Abstract_Network,
        vertex.size = V (Abstract_Network)$Community_size**0.85  ,             # Adjust the size of the super-nodes
        vertex.color = V (Abstract_Network)$Color ,
        vertex.label = NA,
        vertex.label.color = "#000000" ,
        # vertex.label.dist = ifelse(V(Abstract_NetworkPrime)$EdgeBetwCommunity_Size > 25  ,
        #                            1.7 , 1
        # ),
        vertex.label.cex = 0.55 ,
        vertex.label.degree = 75 ,
        vertex.label.family = "Helvetica" ,
        # vertex.frame.color = ifelse(V (Abstract_Network)$Source == "HEALTHYcfDNA" ,
        #                             ifelse (V(Abstract_Network)$Detailed_Annotation == "MDavis" ,
        #                                     "#000000" ,
        #                                     "transparent") ,
        #                             "transparent"), # Add a black border around the super-nodes ,
        vertex.frame.width = 0.1 ,
        vertex.shape = "circle" ,
        edge.width = 0.25 ,
        edge.color = "#000000" ,
        rescale = FALSE ,
        ###To get all the available layouts: grep("^layout_", ls("package:igraph"), value=TRUE)[-1]
        layout = norm_coords(as.matrix(data.frame(V(Abstract_Network)$LayoutX , 
                                                  V(OldAbstract_Network)$LayoutY)) ,
                             ymin=-1,
                             ymax=1,
                             xmin=-1,
                             xmax=1))



# Fig 4B

ggplot(data = GLIPHII_Community_stats %>%
               mutate(Abstract_Annotation = fct_relevel(
                       Abstract_Annotation ,
                       "Cross-species" ,
                       "Human tumour and auto-antigens" ,
                       "HNSCC" ,
                       "Viral" ,
                       "Bacterial" ,
                       "IOKIN_intrinsic" ))) +
        geom_violin(
                aes(
                        x = Abstract_Annotation ,
                        y = Number_of_Patients ) ,
                scale = "count" ,
                # width = 1.75 ,
                # fill = "transparent" ,
                # outlier.shape = NA ,
                linewidth = 0.25
        ) +
        ggbeeswarm::geom_beeswarm(
                aes(
                        x = Abstract_Annotation ,
                        y = Number_of_Patients ,
                        fill = Abstract_Annotation) ,
                method = "center" ,
                cex = 0.5 ,
                shape = 21 ,
                stroke = 0.05 ,
                color = "#000000" ,
                size = 2) +
        scale_fill_manual(
                breaks = c ("Human tumour and auto-antigens" ,
                            "HNSCC" ,
                            "Bacterial" ,
                            "Viral" ,
                            "Cross-species" ,
                            "IOKIN_intrinsic"  ) ,
                labels = c ("Human tumour\n(neo)antigen" ,
                            "HNSCC\nconverging signatures" ,
                            "Bacterial" ,
                            "Viral" ,
                            "Cross-species" ,
                            "IOKIN intrinsic" ) ,
                values = c ("#FF9F01" , #"Human tumour and auto-antigens" ,
                            "#CDEDFA" , #"HNSCC" ,
                            "#BAC70D" , #"Bacterial"
                            "#AC261B" , #"Viral" ,
                            "#FFDF01" , #"Cross-species" ,
                            "#000000"  #"IOKIN_intrinsic" ,
                ),
                name = "Community specificity") +
        guides(fill = guide_legend(override.aes = list(size = 3 ,
                                                       stroke = 0.1) ) )+
        
        
        ylab("Number of participating IO-KIN patients\nin each community") +
        xlab ("GLIPHII-defined specificity") +
        
        scale_x_discrete(breaks = c ("Cross-species" ,
                                     "Human tumour and auto-antigens" ,
                                     "HNSCC" ,
                                     "Viral" ,
                                     "IOKIN_intrinsic" ,
                                     "Bacterial" ) ,
                         labels = c ("Cross-species" ,
                                     "Human\ntumour\n(neo)antigen" ,
                                     "HNSCC\nconverging\nsignatures" ,
                                     "Viral" ,
                                     "IOKIN\nintrinsic" ,
                                     "Bacterial" )) +
        scale_y_continuous(breaks = c(1 , 3 , 5 , 7) ,
                           limits = c(0.9 , 7)) +
        
        
        theme_minimal() +
        theme(panel.grid = element_blank() ,
              panel.spacing.x = unit (1.25 , units = "line") ,
              panel.spacing.y = unit (1.25 , units = "line") ,
              
              axis.line = element_line(color = "#000000" , linewidth = 0.25) ,
              axis.ticks.y = element_line(color = "#000000") ,
              
              axis.title = element_text(size = 9 , color = "#000000") ,
              axis.text.y = element_text(size = 9 , color = "#000000") ,
              axis.text.x = element_blank() ,
              
              legend.title = element_text(size = 9 , color = "#000000") ,
              legend.title.position = "top" ,
              legend.text = element_text(size = 9 , color = "#000000" ) ,
              legend.position = "none",
              legend.justification = "left",
              legend.box.spacing = unit(0.05 , units = "line") ,
              legend.key.height = unit(1.5 , units = "line") ,
              
              aspect.ratio = 0.7 )

# Fig 4C

ggplot(data = GLIPHII_Community_stats %>%
               mutate(Abstract_Annotation = fct_relevel(
                       Abstract_Annotation ,
                       "Cross-species" ,
                       "Human tumour and auto-antigens" ,
                       "HNSCC" ,
                       "Viral" ,
                       "Bacterial" ,
                       "IOKIN_intrinsic" ))) +
        geom_violin(
                aes(
                        x = Abstract_Annotation ,
                        y = Community_size ) ,
                scale = "count" ,
                # width = 1.75 ,
                # fill = "transparent" ,
                # outlier.shape = NA ,
                linewidth = 0.25
        ) +
        ggbeeswarm::geom_beeswarm(
                aes(
                        x = Abstract_Annotation ,
                        y = Community_size ,
                        fill = Abstract_Annotation) ,
                method = "center" ,
                cex = 0.5 ,
                shape = 21 ,
                stroke = 0.05 ,
                color = "#000000" ,
                size = 2) +
        scale_fill_manual(
                breaks = c ("Human tumour and auto-antigens" ,
                            "HNSCC" ,
                            "Bacterial" ,
                            "Viral" ,
                            "Cross-species" ,
                            "IOKIN_intrinsic"  ) ,
                labels = c ("Human tumour\n(neo)antigen" ,
                            "HNSCC\nconverging signatures" ,
                            "Bacterial" ,
                            "Viral" ,
                            "Cross-species" ,
                            "IOKIN intrinsic" ) ,
                values = c ("#FF9F01" , #"Human tumour and auto-antigens" ,
                            "#CDEDFA" , #"HNSCC" ,
                            "#BAC70D" , #"Bacterial"
                            "#AC261B" , #"Viral" ,
                            "#FFDF01" , #"Cross-species" ,
                            "#000000"  #"IOKIN_intrinsic" ,
                ),
                name = "Community specificity") +
        guides(fill = guide_legend(override.aes = list(size = 3 ,
                                                       stroke = 0.1) ) )+
        
        
        ylab("Community size\n(i.e. number of participating TRB-CDR3s)") +
        xlab ("GLIPHII-defined specificity") +
        
        scale_x_discrete(breaks = c ("Cross-species" ,
                                     "Human tumour and auto-antigens" ,
                                     "HNSCC" ,
                                     "Viral" ,
                                     "IOKIN_intrinsic" ,
                                     "Bacterial" ) ,
                         labels = c ("Cross-species" ,
                                     "Human\ntumour\n(neo)antigen" ,
                                     "HNSCC\nconverging\nsignatures" ,
                                     "Viral" ,
                                     "IOKIN\nintrinsic" ,
                                     "Bacterial" )) +
        scale_y_log10(breaks = c(5 , 10 , 20 , 30 , 40 , 50 , 60) ,
                      limits = c(4.95 , 70)) +
        
        
        theme_minimal() +
        theme(panel.grid = element_blank() ,
              panel.spacing.x = unit (1.25 , units = "line") ,
              panel.spacing.y = unit (1.25 , units = "line") ,
              
              axis.line = element_line(color = "#000000" , linewidth = 0.25) ,
              axis.ticks.y = element_line(color = "#000000") ,
              
              axis.title = element_text(size = 9 , color = "#000000") ,
              axis.text.y = element_text(size = 9 , color = "#000000") ,
              axis.text.x = element_blank() ,
              
              legend.title = element_text(size = 9 , color = "#000000") ,
              legend.title.position = "top" ,
              legend.text = element_text(size = 9 , color = "#000000" ) ,
              legend.position = "none",
              legend.justification = "left",
              legend.box.spacing = unit(0.05 , units = "line") ,
              legend.key.height = unit(1.5 , units = "line") ,
              
              aspect.ratio = 0.7 )

### ---------------------------------------------------------------------------------------------------------
### Visualizing the specificity community graphs and TRB CDR3 Sequences' LOGOPlots in Figures 4D-L:
### ---------------------------------------------------------------------------------------------------------
### Example:
ID <- Community ID of interest

nt <- igraph::subgraph(
        network ,
        which(V(network)$LeidenCommunity == ID ))



### Zoomed-in community nodes:
plot.igraph(
        nt,
        vertex.size = ifelse(V (nt)$Source == "IOKINBlood"  , 23 , 15) ,
        vertex.color = ifelse(V (nt)$Source == "IOKINBlood"  , "#000000" , V (nt)$colors ) ,
        vertex.frame.color = "#000000" ,
        #vertex.frame.width = ifelse(V (nt)$Source == "INSPIRE" , 16 , 8) ,
        vertex.frame.width = 0.1 ,
        vertex.shape = "circle",
        
        vertex.label = ifelse(V (nt)$Source == "IOKINBlood"  ,
                              V(nt)$Patient_id ,
                              ifelse(V(nt)$Source == "HNSCC" ,
                                     NA ,
                                     ifelse(V(nt)$Source == "HomoSapiens" ,
                                            sub(".+__(.+)$", "\\1", V(nt)$name) ,
                                            V(nt)$Source))),
        vertex.label.family = "Helvetica" ,
        vertex.label.color = "#000000" ,
        vertex.label.cex = 0.4 ,
        vertex.label.dist = ifelse(V (nt)$Source == "IOKINBlood" ,
                                   11 , 4),
        vertex.label.degree = 0 ,
        
        edge.width = 0.25 ,
        edge.color = "#BABABA" ,
        rescale = FALSE ,
        ###To get all the available layouts: grep("^layout_", ls("package:igraph"), value=TRUE)[-1]
        layout = norm_coords(layout_with_graphopt(nt),
                             ymin=-0.98, ymax=0.98,
                             xmin=-0.98, xmax=0.98) )



### LOGOSeq:

ggplot() +
        geom_logo( V(nt)$TcRb ,
                   method = 'prob' ,
                   seq_type='aa' ,
                   font = "akrobat_regular" ,
                   col_scheme = make_col_scheme(chars = LETTERS,
                                                cols = rep ("#000000" , length(LETTERS))) ) +
        theme_minimal() +
        theme(
                panel.grid = element_blank() ,
                axis.text  = element_blank() ,
                axis.title = element_blank() )
