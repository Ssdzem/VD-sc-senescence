# VD-sc-senescence
Code for reproduction of the experiment in "Single cell analysis dissects the effects of Vitamin D on genetic senescent signatures across murine tissues"


There are 3 main codes: load_input.R, germinal.R and pub_figs.R. 

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

First, load_input.R arranges the directory structure. Then, it downloads the input data from the metadata directory. The metadata directory contains 3 csv: 

-scores.csv contains the senescence gene sets download URL's

-input.csv contains the URL's for downloading the count matrix, barcodes and features of each VD experiment

-ref.csv contains the annotation download URL's for each tissue

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Second, germinal.R is a compilation of all the codes (from preprocessing to cellular communication) to get all the necesary R objects in order to reproduce the results (i.e., the publication figures)

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Finally, pub_figs.R is the code necesary to make all the figures of the article's result. 

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

The rest of the codes in the workflow directory were used for the individual processing of each step of the methods. For example, seurat_analysis_VD.R is the code for the preprocessing, QC and annotation steps of the methods. These codes contain in depth annotation for further understanding of the rationale followed.
