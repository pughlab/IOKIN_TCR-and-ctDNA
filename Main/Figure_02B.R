library(tidyverse)
library(scales)

### Defining the color palette that will be used for each patient -------------------------------------------------
ColorPal_Patient <- tibble(
        Patient_id = c("LIB-25-0012" ,
                       
                       "LIB-25-0018" ,
                       "LIB-25-0017" ,
                       "LIB-25-0009" ,
                       
                       "LIB-25-0013" ,
                       "LIB-25-0004" ,
                       "LIB-25-0002" ,
                       
                       "LIB-25-0014" ,
                       "LIB-25-0007" ,
                       "LIB-25-0008" ,
                       "LIB-25-0010" ,
                       "LIB-25-0011" ,
                       "LIB-25-0001" ,
                       "LIB-25-0005" ,
                       "LIB-25-0003") ,
        Color = c("#08306b" ,
                  
                  "#08519c" ,
                  "#6baed6" ,
                  "#c6dbef" ,
                  
                  "#214C00" ,
                  "#437200" ,
                  "#CADD4C" ,
                  
                  "#FDFD99" ,
                  "#FFD800" ,
                  "#FFB100" ,
                  "#FF8A00" ,
                  "#FE7309" ,
                  "#D0210D" ,
                  "#870110" ,
                  "#440108"))



### Input data preparation ----------------------------------------------------------------------------------------
data_path <- "https://raw.githubusercontent.com/pughlab/IOKIN_TCR-and-ctDNA/refs/heads/main/Data"


### Clinical data:
ClinicalData_fname <- "02_IOKIN_ClinicalData.csv"
ClinicalData <- read_csv(file.path(data_path , ClinicalData_fname)) %>%
        left_join(ColorPal_Patient , by = "Patient_id") %>%
        mutate(Patient_Annotation = paste(Patient_id , "; " , "PFS" , " " , PFS , sep = "")) %>%
        mutate(BestResponse = fct_relevel(BestResponse , "CR" , "PR" , "SD" , "PD"))

### ctDNA Levels:
ctDNAData_fname <- "04_IOKIN_ctDNALevels.csv"
ctDNA <- read_csv (file.path(data_path , ctDNAData_fname)) %>%
        left_join(ClinicalData ,
                  by = "Patient_id")

ctDNALongFormat <- ctDNA %>%
        pivot_longer(cols = starts_with("T") ,
                     names_to = "Timepoint" ,
                     values_to = "ctDNALevel") %>%
        mutate(CycleCode = case_when(
                Timepoint == "T0" ~ 0 ,
                Timepoint == "T1" ~ 2 ,
                Timepoint == "T2" ~ 3 ,
                Timepoint == "T3" ~ 8 ,
                Timepoint == "T4" ~ 15 ,
                Timepoint == "T5" ~ 22 ,
                Timepoint == "T6" ~ 29 ))

### Visualization -------------------------------------------------------------------------------------------------


PatientOrder <- (ClinicalData %>%
                         arrange(match(BestResponse , c("CR" , "PR" , "SD" , "PD")) ,
                                 desc (OS)))$Patient_id



### Figure 2A
ggplot(data = ctDNALongFormat %>%
               mutate(Patient_id = fct_relevel(Patient_id , PatientOrder))) +
        
        
        ### Adding the horizontal line to highlight level 0:
        geom_hline(yintercept = 0 ,
                   linetype = "dashed" ,
                   linewidth = 0.25 ,
                   color = "#000000") +
        #---------------------------------------------------
        geom_line(
                aes(
                        x = CycleCode ,
                        y = ctDNALevel ,
                        group = Patient_id ,
                        color = Patient_id ,
                        linetype = ctDNAShift) ,
                linewidth = 0.4)+
                
        
        geom_point(
                aes(
                        x = CycleCode ,
                        y = ctDNALevel ,
                        fill = Patient_id),
                shape = 21 ,
                size = 1.5 ,
                stroke = 0.15)+
        #---------------------------------------------------

        scale_color_manual(
                values = (ClinicalData %>% arrange(match(Patient_id , PatientOrder)))$Color ,
                breaks = (ClinicalData %>% arrange(match(Patient_id , PatientOrder)))$Patient_id ,
                labels = (ClinicalData %>% arrange(match(Patient_id , PatientOrder)))$Patient_Annotation ,
                name = "Patient ID") +
        scale_fill_manual(
                values = (ClinicalData %>% arrange(match(Patient_id , PatientOrder)))$Color ,
                breaks = (ClinicalData %>% arrange(match(Patient_id , PatientOrder)))$Patient_id ,
                labels = (ClinicalData %>% arrange(match(Patient_id , PatientOrder)))$Patient_Annotation ,
                name = "Patient ID") +
        
        scale_discrete_manual(aesthetics = "linetype" ,
                              breaks = c("Decrease" , "Increase") ,
                              values = c("dashed" , "solid") ,
                              labels = c("\u2193 \u0394 ctDNA", 
                                         "\u2191 \u0394 ctDNA"),
                              name = "ctDNA shift class") +
        
        guides(linetype = guide_legend(title.position = "top", nrow = 1 ) ,
               fill = guide_legend(title.position = "top", direction = "vertical", ncol = 4 ) ,
               color = guide_legend(title.position = "top", direction = "vertical", ncol = 4 )) +
        
        #---------------------------------------------------
        ylab("Mean Tumour Molecule per mL Plasma") +
        xlab ("Time (days)") +
        scale_y_continuous(
                trans = scales::pseudo_log_trans(base = 10),
                breaks = c(0 , 1 , 10 , 100 , 1000 , 10000) ,
                limits = c(-0.05 , 26000) ) +
        scale_x_continuous(breaks = c (0 , 2 , 3 , 8 , 15 , 22 , 29) ) +
        theme_minimal() +
        theme(panel.grid = element_blank() ,
              panel.spacing.x = unit (1.0 , units = "line") ,
              
              axis.line = element_line(color = "#000000" , linewidth = 0.25) ,
              axis.ticks = element_line(color = "#000000") ,
              
              axis.title = element_text(size = 9 , color = "#000000") ,
              axis.text.y = element_text(size = 9 , color = "#000000") ,
              axis.text.x = element_text(size = 9 , color = "#000000") ,
              
              strip.text = element_text(size = 9 , color = "#000000" ) ,
              
              legend.title = element_text(size = 9 , color = "#000000") ,
              legend.text = element_text(size = 9 , color = "#000000" ) ,
              legend.position = "bottom",  
              legend.justification = "left" ,
              legend.box.just = "left" ,
              legend.box = "vertical", 

              aspect.ratio = 1.5) +
        facet_grid(. ~ BestResponse)
