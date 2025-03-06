library(tidyverse)
library(scales)

### Input data preparation ----------------------------------------------------------------------------------------
data_path <- "https://raw.githubusercontent.com/pughlab/IOKIN_TCR-and-ctDNA/refs/heads/main/Data"


### Clinical data:
ClinicalData_fname <- "02_IOKIN_ClinicalData.csv"
ClinicalData <- read_csv(file.path(data_path , ClinicalData_fname)) %>%
        filter(BulkTCRAvailability == "Available") %>%
        mutate(Patient_Annotation = paste(Patient_id , ";" , BestResponse , "\n" , 
                                          "PFS " , PFS , "\n" , 
                                          "OS " , OS , 
                                          sep = "")) 

### Diversity:
diversity_fname <- "05_IOKIN_BuffyCoat_DiversityIndices_DownSampled.csv"
diversity <- read_csv(file.path(data_path , diversity_fname)) %>%
        dplyr::select(Patient_id ,
                      Cycle , Cycle_in_days ,
                      Order_q , Diversity)

### ctDNA Levels:
ctDNAData_fname <- "04_IOKIN_ctDNALevels.csv"
ctDNA <- read_csv (file.path(data_path , ctDNAData_fname)) %>%
        filter(Patient_id %in% unique(ClinicalData$Patient_id)) %>%
        pivot_longer(cols = starts_with("T") ,
                     names_to = "Cycle" ,
                     values_to = "ctDNALevel") 





PatientOrder <- (ClinicalData %>%
        arrange(match(ctDNAShift , c("Decrease" , "Increase")) ,
                desc(OS)) )$Patient_Annotation



VisTable <- full_join(diversity ,
                      ctDNA ,
                      by = c("Patient_id" , "Cycle")) %>%
        left_join(ClinicalData ,
                  by = "Patient_id") %>%
        mutate(Patient_Annotation = fct_relevel(Patient_Annotation , PatientOrder)) %>%
        mutate(ctDNAShift = case_when(
                ctDNAShift == "Decrease" ~ "\u2193 \u0394 ctDNA",
                ctDNAShift == "Increase" ~ "\u2191 \u0394 ctDNA")) %>%
        mutate(ctDNAShift = fct_relevel(ctDNAShift , "\u2193 \u0394 ctDNA", "\u2191 \u0394 ctDNA")) 

### Visualization -------------------------------------------------------------------------------------------------

# Fig 3A


ggplot(data = VisTable ,
       aes(x = Cycle_in_days )) +
        #---------------------------------------------------------
        geom_hline(yintercept = -1 ,
                   color = "#000000" ,
                   linetype = "dashed" ,
                   linewidth = 0.25) +
        #---------------------------------------------------------
        geom_line(
                aes(y = log10 (ctDNALevel + 0.1) ,
                    group = Patient_id) ,
                color = "#000000",
                linewidth = 0.25) +
        
        geom_point(
                aes(y = log10 (ctDNALevel + 0.1)) ,
                shape = 19 ,
                color = "#000000",
                size = 1.0 ,
                stroke = 0) +
        #---------------------------------------------------------
        geom_line(data = VisTable %>% filter(Order_q %in% c(0 , 1)) ,
                aes(x = Cycle_in_days ,
                    y = Diversity/250 ,
                    group = Order_q ,
                    linetype = as.factor(Order_q)) ,
                color = "#71980B",
                linewidth = 0.25 ) +

        
        geom_point(data = VisTable %>% filter(Order_q %in% c(0 , 1)) ,
                aes(x = Cycle_in_days ,
                    y = Diversity/250 ) ,
                shape = 19 ,
                color = "#71980B",
                size = 1.0 ,
                stroke = 0) +
        
        scale_discrete_manual(aesthetics = "linetype" ,
                              breaks = c("0" , "1") ,
                              values = c("dashed" , "solid") ,
                              labels = c("Richness" , "Shannon diversity") ,
                              name = "Diversity metric") +
        #---------------------------------------------------------

        ylab("Mean Tumour Molecule per mL Plasma") +
        xlab ("Time (days)") +
        scale_y_continuous(
                breaks = c(-1 , 0 , 1 , 2 ) ,
                labels = c(0 , 1 , 10 , 100 ) ,
                
                sec.axis = sec_axis(trans = ~ .*1,
                                    name = "Peripheral TRB-CDR3\nrepertoire diversity" ,
                                    breaks = c(30/250 , 230/250 , 430/250 , 630/250),
                                    labels = c(30 , 230 , 430 , 630))) +
        scale_x_continuous(breaks = c (0 , 2 , 3 , 8 , 15 , 22 , 29) ) +
        
        ggh4x::facet_nested_wrap( vars(ctDNAShift , Patient_Annotation) ,
                                  nrow = 2 ,
                                  dir = "h" ) +
        theme_minimal() +
        theme(panel.grid = element_blank() ,
              panel.spacing.x = unit (0.5 , units = "line") ,
              panel.spacing.y = unit (0.1 , units = "line") ,
              
              axis.line = element_line(color = "#000000" , linewidth = 0.25) ,
              axis.line.y.right = element_blank() ,
              
              axis.ticks = element_line(color = "#000000") ,
              axis.ticks.y.right = element_line(color = "#71980B") ,
              
              axis.title = element_text(size = 9 , color = "#000000") ,
              axis.title.y.right = element_text(size = 9 , color = "#71980B") ,
              
              axis.text = element_text(size = 9 , color = "#000000") ,
              axis.text.y.right = element_text(size = 9 , color = "#71980B") ,
              
              strip.text = element_text(size = 9 , color = "#000000" ) ,
              
              legend.title = element_text(size = 9 , color = "#000000") ,
              legend.text = element_text(size = 9 , color = "#000000" ) ,
              legend.position = "top",
              legend.justification = "left",
              legend.box.spacing = unit(0.05 , units = "line") )
