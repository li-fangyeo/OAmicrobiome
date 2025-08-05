# Orang Asli Microbiome
This study investigates the gut microbiome diversities, predicted pathways and resistomes of three indigenous OA, Jahai, Temiar, and Temuan communities with different lifestyles and degree of urbanisation. We included an urban Malay group as comparison. 

## How to read my codes
Mia.Rmd holds the bulk of the anlaysis, namely alpha, beta, linear models, and figures. 

Then func.Rmd holds the predicted pathway analysis, including prevalence filtering, transforming data to either dichotomised, or inverse-rank normalised, and then using linear models to test for association. 

functionCircular.R has the code to produce the circular figure for pathway analysis. And table.R is just there for ease of producing publication-ready tables.

AMR_scripts folder has all the codes for the resistome analysis, including counting ARG load and diversity, and ARG differential abundance analysis using linear models. 

