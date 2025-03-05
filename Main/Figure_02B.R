library(tidyverse)
library(scales)

### Input data preparation ----------------------------------------------------------------------------------------
data_path <- "https://raw.githubusercontent.com/pughlab/IOKIN_TCR-and-ctDNA/refs/heads/main/Data"

### Clinical data:
ClinicalData_fname <- "02_IOKIN_ClinicalData.csv"
ClinicalData <- read_csv(file.path(data_path , ClinicalData_fname)) %>%
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
