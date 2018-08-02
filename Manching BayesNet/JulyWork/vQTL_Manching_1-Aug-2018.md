<!--- pandoc vQTL_Manching_1-Aug-2018.md -f markdown -t html -s -o vQTL_Manching_1-Aug-2018.html --->


# vQTL Analysis on 2012 Manching Plant Height Data.

**Original Data Set:** Manching2012PlantHT.csv, IBM94markerset raw data.csv

## Process

- Clean and format the data to be in the *cross* format needed for scanonevar in vQTL. This results in the ManchingScrubbed.csv file.

- Run the ManchingScrubbed file through the R script vQTL_Manching.R which returns the vQTL results.

vQTL_Manching.R can be found at the following GitHub [link](https://github.com/michael-byrd/Stapleton-Lab/blob/development/Manching%20BayesNet/JulyWork/vQTL_Manching.R).

### 29/July/2018

**Output File:** Manching2012_vQTL_Analysis.csv

- These results use the following scanonevar function which is corrected later. `outv <- scanonevar(cross = dat, mean.formula = Height ~ Low.Water + Low.Nitrogen + Pathogen + mean.QTL.add + mean.QTL.dom, var.formula = ~ var.QTL.add + var.QTL.dom)`

### 31/July/2018

**Output File:** Manching_vQTL_July-31-18.csv

- The scanonevar function was corrected removing the .dom terms and adding that covar.effects return True

`outv <- scanonevar(cross = dat, mean.formula = Height ~ Low.Water + Low.Nitrogen + Pathogen + mean.QTL.add, var.formula = ~ var.QTL.add, return.covar.effects = TRUE)`
