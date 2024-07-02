# Data analysis scripts for the article "Metabolic-based Polyclonality in a Drosophila Intestinal Tumor Model" by Pierre Delamotte, Mickael Poidevin, Yan Jazczyszyn, Arnaud Le Rouzic, and Jacques Montagne

This repository contains the data and the code (in R) to reproduce two analyses from the article:
* The FACS data analysis (figure 2) and the motility data analysis (figure 7)
* The statistical model for the number of clones (figure 6F)

## Organization of the repository
* /src: R code for the analysis
  * FACS.R: runs the analysis of the FACS data
  * Polyclonality.R: runs the analysis on the number of clones. The file Polyclonality-tools.R contains helper functions. 
* /results: empty directory where the tables and figures will be generated
  * For the FACS analysis:
    * FACS.txt:     Estimates, Standard Errors, raw p-values, and adjusted p-values for all genes, contrasted with the controls
    * Motility.txt: Estimates, Standard Errors, raw p-values, and adjusted p-values for all genes, contrasted with the controls
  * For the Clone analysis
    * poly3.pdf: estimated frequencies of the 4 colors as a function of the heat shock duration and the time before measurement
    * poly4.pdf: estimated number of clines / tumor as a function of the heat shock duration
    * poly6.pdf: data vs predictions for the frequency of each color combination category
    * poly7.pdf: frequency of the number of observations (tumors/non-tumors) as a function of the heat shock duration
* /data:
  * FACSdata.xlsx: Excel file containing the FACS data set. 4 columns: Genotype, number of nuclei, number of GFP-tagged cells, and the Percentage of GFP cells
  * dataMotility.xlsx: same as for FACS
  * FlybowPolyclonality.xlsx: Excel file containing the count categories. The two first columns (Experiment and Name) are tags, Time HS is the heat shock duration, Days post-HS is the time before measurement (less than 21 days: short, more than 21 days; long), four color columns for the non-tumor cells (with approximate counts for the number of clones), Presence of tumor (T/F), and the presence/absence of each color in tumors (if Tumor=T).

## Model information
* FACS and Motility
  Two statistical models were explored. Model1 is a "simple" linear model, using the proportion of GFP cells as a dependent variable; p-values of differences were computed by t-tests. Model3 is a "quasibinomial" glm, considering that the ratio of GFP cells follows an overdispersed binomial distribution. The model was fit by maximum likelihood, and p-values were computed out of t tests. All p-values were adjusted by the Holm-Bonferroni method.
* Clone analysis
  The details of the probabilistic model are provided in the supplementary methods of the manuscript.  
