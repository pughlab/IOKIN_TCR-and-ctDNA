library(tidyverse)

### Input data preparation ----------------------------------------------------------------------------------------
data_path <- "https://raw.githubusercontent.com/pughlab/IOKIN_TCR-and-ctDNA/refs/heads/main/Data"

### Clinical data:
ClinicalData_fname <- "02_IOKIN_ClinicalData.csv"
ClinicalData <- read_csv(file.path(data_path , ClinicalData_fname)) %>%
        dplyr::select(Patient_id , BestResponse)

### Sample Inventory and Technologies:
SampleInventory_fname <- "01_IOKIN_SampleInventory-and-Technologies.csv"
SampleInventory <- read_csv(file.path(data_path , SampleInventory_fname))

### Visualization -------------------------------------------------------------------------------------------------
###Patient order:
PatientOrder <- (SampleInventory %>%
                         group_by(Patient_id) %>%
                         summarise(N = n()) %>%
                         ungroup() %>%
                         arrange(N , desc (Patient_id)))$Patient_id


### Figure 1B

ggplot(data = SampleInventory %>%
               left_join(ClinicalData , by = "Patient_id") %>%
               mutate(Technology = fct_relevel(Technology ,
                                               "WES" , "RNA-seq" , "Bespoke ctDNA" , "CapTCR-seq") ,
                      Patient_id = fct_relevel(Patient_id , PatientOrder) ,
                      Timepoint = fct_relevel(Timepoint ,
                                              "Archival" , "Baseline" ,
                                              "Day 2" , "Day 3" , "Day 8" , "Day 15" ,
                                              "Day 22" , "Day 29"),
                      BestResponse = fct_relevel(BestResponse ,
                                                 "CR" , "PR" , "SD" , "PD"))) +
        geom_tile(
                aes(
                        x = Timepoint ,
                        y = Patient_id ,
                        fill = `Sample type`),
                color = "#000000" ,
                linewidth = 0.05 ,
                width = 0.85 , height = 0.75 )+
        
        scale_fill_manual(breaks = c("Tumour" , "Plasma ctDNA" , "Buffy Coat"),
                          labels = c("Tumour" , "Plasma ctDNA" , "Buffy Coat") ,
                          values = c("#5E7814" , "#FCA505" , "#AA2600")) +
        
        
        facet_grid(BestResponse ~ Technology ,
                   scales = "free" ,
                   space = "free" ,
                   labeller = as_labeller(c(`CR` = "CR" ,
                                            `PR` = "PR" ,
                                            `SD` = "SD" ,
                                            `PD` = "PD" ,
                                            `WES` = "WES" ,
                                            `RNA-seq` = "RNA-\nseq" ,
                                            `Bespoke ctDNA` = "Bespoke ctDNA" ,
                                            `CapTCR-seq` = "CapTCR-seq"))) +
        theme_minimal() +
        theme(
                panel.grid = element_blank() ,
                panel.spacing = unit(0.5 , units = "line"),
                axis.ticks.x = element_line(linewidth = 0.25 , color = "#000000"),
                axis.line.x = element_line(color = "#000000" , linewidth = 0.25) ,
                
                
                axis.title = element_blank() ,
                
                axis.text.y = element_text(size = 9 , color = "#000000") ,
                axis.text.x = element_text(size = 9 ,  color = "#000000" , angle = 30 , hjust = 1) ,
                
                strip.text.x = element_text(size = 9 , color = "#000000") ,
                strip.text.y = element_text(size = 9 , color = "#000000" , angle = 0) ,
                
                legend.position = "bottom" ,
                legend.text = element_text(size = 9 , color = "#000000") ,
                legend.title = element_text(size = 9 , color = "#000000"))
