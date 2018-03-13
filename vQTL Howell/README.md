# vQTL Analysis for the Howell Lab
## Overview
This is to provide an in depth record of how this analysis was conducted for those interested in this project.
  
## Initial Data
Our initial Dataset can be was the file [Heat stress.csv](Heat%20stress.csv). 
It contained three columns of data consisting of taxa identification numbers as well as Spliced and Unspliced bZIP60 information.
These taxa are representative of corn lines with very precise genotypes. So our first step was translating these into corresponding genetic information.

### TASSEL5
To do this, Dr. Ann Stapleton directed me to [this website](http://cbsuss05.tc.cornell.edu/hdf5new/query.asp "Genetic Data") to get the genetic data for our specific taxa.
Given that the taxa list to enter was in a format that was more extended than that which was provided in [Heat stress.csv](Heat%20stress.csv), I had to use the [Taxa List](http://cbsuss05.tc.cornell.edu/hdf5new/taxa.asp) found on the same website.
After using a quick [R Script](Full_Zs.R "Converter") to to convert our short names to extended names, I got the genetic data emailed to me.
Once the file was downloaded, I used a program called **TASSEL5** and applied the following filters:
* Many of the Taxa had duplicate rows, so whenever there were two entries that shared a row from the [Heatstress.csv](Heat%20stress.csv), I removed the second one.
* As far as sites go, they needed to have the following characteristics:
	* A site minimum count of 132 (total number of taxa after the first filter)
	* A minimum allele frequency of .1 to eliminate monomorphisms (the .1 number was suggested by Dr. Stapleton)
	* I also filtered the Minor SNP States, as they will not be very useful for prediction

The resulting file is [here](Howell_scrubbed_Z_to_SNPs.txt) and is used to stitch together our preliminary data set.

### Data Cleaning
We begin by replacing the Z values in the [Heat stress](Heat%20stress.csv) file with their corresponding genetic data.
This is done by referring to the [Genetic Data](Howell_scrubbed_Z_to_SNPs.txt) from the previous step.
It is also important to note that the site locations are assumed to be in basepairs. The location was taken from the name in TASSEL.
The code used for this is found [here](HowellvQTL.R) in the first half of the script. This output an object called [Howell-Cross-Object.csv](Howell-Cross-Object.csv).
From here we used the [Howell Sample.R](Howell%20Sample.R) file to clean the data three times. 
1. For some reason TASSEL did not get rid of some minor allele states, so we set them to "N"
2. The vQTL analysis requires the alleles to be in either an "A" or "B" state, so we mapped the A,T,C, and G as follows:  
	* A -> A  
	* T -> B  
	* C -> A  
	* G -> B  
3. There were still many columns that were not useful for prediction. For example, a column consisting of only "A"s and "N"s. So I removed all the columns that did not have three unique entries.

Finally [this object](Howell-Cross-ObjectC3.csv "Howell-Cross-ObjectC3") was ready for vQTL analysis on the Unspliced_bZIP60 and Spliced_bZIP60.
Additionally, I made a file to predict the ratio between the two (Unspliced/Spliced) labeled [Howell-Cross-Object-Ratio.csv](Howell-Cross-Object-Ratio.csv "Ratio Object").

### vQTL Analysis
Finally, after uploading the necessary files to my Stampede2 environment. I ran the following 3 scripts to generate our output for [Spliced](Spliced.R), [Unspliced](Unspliced.R),
and [Ratio](Ratio.R). To submit these to Stampede2, I wrote some quick shell scripts and submitted them with the the sbatch function. 

### Output
After the analysis ran, I took the output csv files
off of Stampede2 and the [Spliced](HowellvQTL_Spliced_LOD,Pvals,EffectSizes.csv),[Unspliced](HowellvQTL_Unspliced_LOD,Pvals,EffectSizes.csv), and [Ratio](HowellvQTL_Ratio_LOD,Pvals,EffectSizes.csv)
can be found here as well.