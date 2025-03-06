library(tidyverse)
library(scales)

### Input data preparation ----------------------------------------------------------------------------------------
data_path <- "https://raw.githubusercontent.com/pughlab/IOKIN_TCR-and-ctDNA/refs/heads/main/Data"


### Clinical data:
ClinicalData_fname <- "02_IOKIN_ClinicalData.csv"
ClinicalData <- read_csv(file.path(data_path , ClinicalData_fname)) %>%
        filter(BulkTCRAvailability == "Available") %>%
        dplyr::select(Patient_id , ctDNAShift)

### Diversity:
diversity_fname <- "05_IOKIN_BuffyCoat_DiversityIndices_DownSampled.csv"
diversity <- read_csv(file.path(data_path , diversity_fname)) %>%
        filter(Order_q %in% c(0 , 1)) %>%
        dplyr::select(Patient_id ,
                      Cycle , Cycle_in_days ,
                      Order_q , Diversity) %>%
        left_join(ClinicalData , by = "Patient_id") %>%
        mutate(Cycle_in_days = case_when(
                Cycle_in_days == 2 ~ 1 ,
                Cycle_in_days == 3 ~ 4 ,
                TRUE ~ Cycle_in_days))

diversitySpezial <- left_join(
        diversity %>%
                filter(Cycle != "T0") ,
        diversity %>%
                filter(Cycle == "T0") %>%
                dplyr::select(Patient_id , Order_q , Diversity) %>%
                pivot_wider(values_from = Diversity ,
                            names_from = Order_q) %>%
                rename(BaselineRichness = `0` ,
                       BaselineShannon = `1`) ,
        by = "Patient_id")

### Visualization -------------------------------------------------------------------------------------------------

THEME <- theme_minimal() +
        theme(panel.grid = element_blank() ,
              panel.spacing.x = unit (1.25 , units = "line") ,
              panel.spacing.y = unit (1.25 , units = "line") ,
              
              axis.line = element_line(color = "#000000" , linewidth = 0.25) ,
              axis.ticks = element_line(color = "#000000") ,
              
              axis.title = element_text(size = 9 , color = "#000000") ,
              axis.text = element_text(size = 9 , color = "#000000") ,
              aspect.ratio = 0.6)

# Fig 3B

ggplot(data = diversitySpezial %>%
               filter(Order_q == 0)) +
        
        geom_hline(yintercept = 0 ,
                   linetype = "dashed" ,
                   linewidth = 0.25 ,
                   color = "#ACABA8") +
        #--------------------------------------------------

        geom_line(data = diversitySpezial %>%
                          filter(Order_q == 0) %>%
                          dplyr::select(Patient_id ,
                                        Cycle , Cycle_in_days ,
                                        BaselineRichness , Diversity) %>%
                          group_by(Cycle_in_days) %>%
                          mutate(MedianRichnessValue = median((Diversity - BaselineRichness) / `BaselineRichness`)) %>%
                          ungroup() %>%
                          dplyr::select(- BaselineRichness , - Diversity , - Patient_id) %>%
                          unique() ,
                  aes(x = Cycle_in_days , y = MedianRichnessValue) ,
                  linewidth = 1 ,
                  lineend = "round" ,
                  color = "#564635" ,
                  alpha = 0.25) +
        #--------------------------------------------------

        ggdist::stat_pointinterval(
                data = diversitySpezial %>%
                        filter(Order_q == 0) %>%
                        filter(ctDNAShift == "Decrease"),
                aes(
                        x = Cycle_in_days - 0.3 ,
                        y = (Diversity - BaselineRichness) / `BaselineRichness`
                ),
                alpha = 1 ,
                # position = position_dodge(width = -0.5) ,
                interval_size_range = c(0.1 , 0.6) ,
                .width = c(.5, 1) ,
                point_size = 0.5 ,
                color = "#5C2F01")  +
                
        ggdist::geom_weave(
                data = diversitySpezial %>%
                        filter(Order_q == 0) %>%
                        filter(ctDNAShift == "Decrease"),
                aes(
                        x = Cycle_in_days - 0.45 ,
                        y = (Diversity - BaselineRichness) / BaselineRichness
                ),
                fill = "#5C2F01" ,
                side = "left" ,
                # layout = "bin",
                overflow = "compress" ,
                slab_linewidth = 0.025 ,
                slab_color = "#5C2F01" ,
                # stackratio = 0.3 ,
                # position = position_dodge(width = - 0.5) ,
                dotsize = 1 ,
                alpha = 0.75 ) +
        
        #--------------------------------------------------

        ggdist::stat_pointinterval(
                data = diversitySpezial %>%
                        filter(Order_q == 0) %>%
                        filter(ctDNAShift == "Increase"),
                aes(
                        x = Cycle_in_days + 0.3 ,
                        y = (Diversity - BaselineRichness) / BaselineRichness
                ),
                alpha = 1 ,
                # position = position_dodge(width = -0.5) ,
                interval_size_range = c(0.1 , 0.6) ,
                .width = c(.5, 1) ,
                point_size = 1 ,
                color = "#E4B200")  +
        
        ggdist::geom_weave(
                data = diversitySpezial %>%
                        filter(Order_q == 0) %>%
                        filter(ctDNAShift == "Increase"),
                aes(
                        x = Cycle_in_days + 0.85 ,
                        y = (Diversity - BaselineRichness) / BaselineRichness
                ),
                fill = "#E4B200" ,
                side = "left" ,
                # layout = "bin",
                overflow = "compress" ,
                slab_linewidth = 0.025 ,
                slab_color = "#000000" ,
                # stackratio = 0.3 ,
                # position = position_dodge(width = - 0.5) ,
                dotsize = 1 ,
                alpha = 0.75 ) +
        
        #-------------------------------------------------------------------------------------------

        xlab("Time (days)") +
        ylab (expression(frac("Day"["n"]*" Richness - baseline Richness",
                              "baseline Richness")))+
        
        scale_x_continuous(breaks = c (1 , 4 , 8 , 15 , 22 , 29),
                           labels = c ("2" , "3" , "8" , "15" , "22" , "29")) +
        scale_y_continuous (breaks = c (-0.5 , 0 , 0.5 , 1 ) ,
                            limits = c(-0.75 , 1.25)) +
        THEME
        
        



# Fig 3C

ggplot(data = diversitySpezial %>%
               filter(Order_q == 1)) +
        
        geom_hline(yintercept = 0 ,
                   linetype = "dashed" ,
                   linewidth = 0.25 ,
                   color = "#ACABA8") +
        #--------------------------------------------------

        geom_line(data = diversitySpezial %>%
                          filter(Order_q == 1) %>%
                          dplyr::select(Patient_id ,
                                        Cycle , Cycle_in_days ,
                                        BaselineShannon , Diversity) %>%
                          group_by(Cycle_in_days) %>%
                          mutate(MedianShannonValue = median((Diversity - BaselineShannon) / `BaselineShannon`)) %>%
                          ungroup() %>%
                          dplyr::select(- BaselineShannon , - Diversity , - Patient_id) %>%
                          unique() ,
                  aes(x = Cycle_in_days , y = MedianShannonValue) ,
                  linewidth = 1 ,
                  lineend = "round" ,
                  color = "#564635" ,
                  alpha = 0.25) +

        #--------------------------------------------------

        ggdist::stat_pointinterval(
                data = diversitySpezial %>%
                        filter(Order_q == 1) %>%
                        filter(ctDNAShift == "Decrease"),
                aes(
                        x = Cycle_in_days - 0.3 ,
                        y = (Diversity - BaselineShannon) / `BaselineShannon`
                ),
                alpha = 1 ,
                # position = position_dodge(width = -0.5) ,
                interval_size_range = c(0.1 , 0.6) ,
                .width = c(.5, 1) ,
                point_size = 0.5 ,
                color = "#5C2F01")  +
        
        ggdist::geom_weave(
                data = diversitySpezial %>%
                        filter(Order_q == 1) %>%
                        filter(ctDNAShift == "Decrease"),
                aes(
                        x = Cycle_in_days - 0.45 ,
                        y = (Diversity - BaselineShannon) / BaselineShannon
                ),
                fill = "#5C2F01" ,
                side = "left" ,
                # layout = "bin",
                overflow = "compress" ,
                slab_linewidth = 0.025 ,
                slab_color = "#5C2F01" ,
                # stackratio = 0.3 ,
                # position = position_dodge(width = - 0.5) ,
                dotsize = 1 ,
                alpha = 0.75 ) +
        
        #--------------------------------------------------

        ggdist::stat_pointinterval(
                data = diversitySpezial %>%
                        filter(Order_q == 1) %>%
                        filter(ctDNAShift == "Increase"),
                aes(
                        x = Cycle_in_days + 0.3 ,
                        y = (Diversity - BaselineShannon) / BaselineShannon
                ),
                alpha = 1 ,
                # position = position_dodge(width = -0.5) ,
                interval_size_range = c(0.1 , 0.6) ,
                .width = c(.5, 1) ,
                point_size = 0.5 ,
                color = "#E4B200")  +
        
        ggdist::geom_weave(
                data = diversitySpezial %>%
                        filter(Order_q == 1) %>%
                        filter(ctDNAShift == "Increase"),
                aes(
                        x = Cycle_in_days + 0.85 ,
                        y = (Diversity - BaselineShannon) / BaselineShannon
                ),
                fill = "#E4B200" ,
                side = "left" ,
                # layout = "bin",
                overflow = "compress" ,
                slab_linewidth = 0.025 ,
                slab_color = "#000000" ,
                # stackratio = 0.3 ,
                # position = position_dodge(width = - 0.5) ,
                dotsize = 1 ,
                alpha = 0.75 ) +

        #--------------------------------------------------

        xlab("Time (days)") +
        ylab (expression(frac("Day"["n"]*" Shannon - baseline Shannon",
                              "baseline Shannon diversity")))+
        
        scale_x_continuous(breaks = c (1 , 4 , 8 , 15 , 22 , 29),
                           labels = c ("2" , "3" , "8" , "15" , "22" , "29")) +
        scale_y_continuous (breaks = c (-0.5 , 0 , 0.5 , 1 ) ,
                            limits = c(-0.75 , 1.25)) +
        
        THEME
