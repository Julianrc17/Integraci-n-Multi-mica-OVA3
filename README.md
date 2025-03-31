# README: Multi-Omic Analysis of High Grade Serous Ovarian Cancer (HGSOC)

## 1. Description of the Data and Study Objectives

This study analyzes multi-omic data from High Grade Serous Ovarian Cancer (HGSOC) patients, obtained from The Cancer Genome Atlas (TCGA). The primary goal is to explore molecular differences between the immunoreactive and mesenchymal subtypes by integrating multiple layers of omics data.

The study includes the following omics datasets:
- **RNA-Seq (Gene Expression):** Examines gene expression patterns in different HGSOC subtypes.
- **Copy Number Variations (CNV):** Identifies genomic alterations based on changes in DNA copy number.
- **DNA Methylation:** Provides an epigenetic perspective, highlighting potential differences between cancer subtypes.
- **Transcription Factor Expression (TF):** Investigates the activity of transcriptional regulators.

### Study Objectives
1. **Identify a multi-omic signature** that distinguishes immunoreactive from mesenchymal cancer subtypes using DIABLO, an integrative method for multi-omic analysis.
2. **Understand the regulatory mechanisms** of gene expression mediated by CNVs, DNA methylation, and transcription factors. Visualization tools help detect key regulatory patterns.
3. **Use additional tools** such as Cytoscape, PCA, and MultiPower to complement DIABLO and enhance the interpretation of omics interactions.

## 2. Data Processing and Preprocessing

The dataset includes 149 patients for multi-omic analysis, focusing on immunoreactive and mesenchymal subtypes. Given the nature of the study, a vertical integration approach is used, as different omics are measured on the same individuals.

### Boxplot Visualization
Exploratory data analysis confirms minimal variation within each omics type. CNV exhibits greater variability, likely due to genomic instability and tumor heterogeneity.

## 3. Principal Component Analysis (PCA)

PCA is used to explore the main sources of variability and detect potential batch effects. The results indicate that no significant separation exists between the two cancer subtypes based on individual omics datasets.

## 4. Statistical Power Calculation

Statistical power analysis reveals that DNA methylation is the limiting omics dataset, requiring a larger sample size to achieve optimal power for detecting significant differences.

## 5. Multi-Omic Integration with DIABLO

DIABLO is applied to identify a multi-omic signature distinguishing cancer subtypes. Results indicate limited correlation between the omics datasets, with CNV and TF showing the strongest associations. However, the model does not provide a clear separation between subtypes.

## 6. Network Analysis with Cytoscape

To further explore omics interactions, Cytoscape is used for network visualization. The analysis highlights key regulatory interactions, particularly between CNV and TF.

## 7. Conclusions

- No significant molecular differences were found between immunoreactive and mesenchymal subtypes.
- CNV and TF show stronger associations compared to methylation, which has minimal impact.
- DIABLO provides a useful multi-omic signature, but its predictive power is limited.
- Cytoscape highlights key gene regulatory interactions, particularly involving transcription factors like SOX11.

This study suggests that additional omics layers or alternative integration methods may be needed to uncover subtype-specific patterns in HGSOC.

