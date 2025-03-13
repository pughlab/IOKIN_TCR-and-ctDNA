# Investigating Early Kinetics in Plasma ctDNA and Peripheral T Cell Receptor Repertoire to Predict Treatment Outcomes to PD-1 Inhibitors
## Abstract
PD-1:PD-L1 axis targeting immune checkpoint blockade (ICB) significantly improves survival in patients with recurrent/metastatic head and neck squamous cell carcinoma (R/M HNSCC). Circulating tumor DNA (ctDNA) and peripheral T cell receptor (TCR) repertoires are emerging as potential predictive biomarkers of response to ICB and understanding their early dynamics during treatment could guide treatment adjustments before progression is evident radiologically. IO-KIN (NCT04606940) is a prospective study of R/M HNSCC patients treated with nivolumab or pembrolizumab. 104 blood samples were collected across 7 time points, spanning baseline and the first 29 days after the initial ICB dose. We examined the early dynamics of ctDNA (N = 15 patients), and the peripheral TCR repertoire (N = 8 patients). While baseline ctDNA levels did not correlate with clinical outcomes, a decrease in ctDNA after day 8 was associated with clinical benefit, prolonged progression-free survival, and a trend toward improved overall survival. TCR repertoire transiently diversified in response to ICB, peaking between days 8-22. GLIPHII (Grouping Lymphocyte Interactions by Paratope Hotspots II) analysis identified signatures recognizing tumor-associated antigens at baseline and emerging in response to ICB as early as three days. We demonstrated that simultaneous monitoring of ctDNA and TCR dynamics during treatment allows for more personalized and timely treatment strategies.

### Key words:
Head and Neck Squamous Cell Carcinoma (HNSCC), Recurrent/Metastatic (R/M), Immune Checkpoint Blockade (ICB), Circulating Tumor DNA (ctDNA), T-Cell Receptor (TCR), Peripheral Blood

### Documents in the data directory:
##### 01_IOKIN_SampleInventory-and-Technologies.csv
Input document to generate Figure 1B Tile plot. <br>

02_IOKIN_ClinicalData.csv <br>
03_IOKIN_BuffyCoat_TRB-CDR3.csv <br>
04_IOKIN_ctDNALevels.csv <br>
05_IOKIN_BuffyCoat_DiversityIndices_DownSampled.csv <br>
06_IOKIN_BuffyCoat_ClonalityIndices_DownSampled.csv <br>

### Session Info:
R version 4.4.1 (2024-06-14) Running under: macOS Big Sur 11.7.10 <br>
RStudio Version 2023.09.1+494 <br>
Package tidyverse version 2.0.0 <br>
Package mgcv version 1.9-1 <br>
Package gratia version 0.10.0 <br>
Package scales version 1.3.0 <br> 
Package survminer version 0.5.0 <br>
Package survival version 3.8-3 <br>
Package igraph version 2.1.4 <br>
Package graphlayouts version 1.2.2

