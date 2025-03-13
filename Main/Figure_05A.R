library(tidyverse)
library(scales)
library(ggh4x)
library(ggdist)

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

PatientOrder <- (ClinicalData %>%
        arrange(match(ctDNAShift , c("Decrease" , "Increase")) ,
                desc(OS)) )$Patient_Annotation
