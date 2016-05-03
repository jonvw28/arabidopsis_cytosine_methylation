# Towards a Predictive Model of DNA methylation in *Arabidopsis thaliana*

 Repository for code used in preparing the report for my masters project

*Jonathan Williams*

*03/05/2016*

The code included in this repository is the code used in creating my project report for the above investigation. Included in the sub-directory figures is the code necessary to prepare all the figures in the report. The data in the tables of the report can likewise be recreated by running the scripts in the sub-directory tables. For those who are interested, the code that forms the skeleton of the model presented in the research report is also available in the sub-directory models (this is liable to be updated should I feel like completing some more work on these).

## Preparatory Code

Before being able to run the routines presented here you will need a set of data files which are too large to share here. In addition, some of these are unpublished experimental files prepared by Dr. Marco Catoni of the Sainsbury Laboratory at the University of Cambridge. As such, should you wish to run these routines, you can contact me at jonvw28@gmail.com and I can get permission to share the necessary data with you.

Once you have the data files you will need to place them, without editing their names, into the sub-directory Data. 
To then turn these raw data files into the final processed data used in the project (as per the early sections of the results section) you will then need to run the script data_preprocessing.R This will use the raw data files to complete and save a finalised data file. In particular this will.:

* Coerce the Granges object into a data frame
* Add columns indicating whether the tiles exceed the CpG, CHG and CHH thersholds for methylation in these contexts
* Add a relative methylation column as per the report
* Add columns indicating which side of the definitions for ETELs and RTELs each tile fall
* Add columns indicating if the tile overlaps an annotated gene, as well as for many other feaures from the TAIR 10 annotation
* Add columns indicating if the tile overlaps an annotated Transposable element or tranposon from the TAIR 10 annotation
* Add columns for 3 replicates of expression data collected for WT plants and calculate an average of these
* Add columns for 3 replicates of expression data collected for *met*1-3 mutants and calculate an average of these
* Add a column for relative expression between the above
* Add a column for calculated mappability scores as per the methodology laid out in the report
* Add a column for the simplified annotation for each tile as per the methodology described in the report.

The final data frame that is produced by this will then be saved to the top directory. This data frame can likewise be requested subject to permission.

In order to run this data processing requires some of the functions I have written for this project. All of these functions are included in the script functions_read_in_first.R. All scripts in this repository should source this script directly, but if not, you can manually add all of the functions to your global environment by simply running this script. Not all of these functions will necessarily be needed as some are now defunct for the final direction the project took. Some functions are also just lazy shortcuts. They should all have comments explaining what they do. Some of these may need updating in due course so there may well be revisions to the descriptions, but the functional code will be locked in as at the time of submission of this project.

## Figures

Where figures in the report were produced from the data prepared here, the R script for each figure has been included in the sub-directory figures. In the most part, the scripts will output the figure to the default graphic device, though a simple edit can output it to a graphic device of your choice. The only exception is for figure 7 where a large number of plots will automatically be outputted to the working directory. From these, the 5 quantile CpG count and mappability plots are the figures used in the final report

## Tables

Where tables of data have been included in the report, scripts to generate this data have been included in the sub-directory tables. These should all, in general, generate more data than needed for the table, and should print the result given in the table automatically via the final function call. These can of course be modified towards the end to display any result of interest to you. 

The main exception to this are tables 1 & 2 where the scripts output linear models. To recover the data presented in the report simply pass these to the function summary() in R and you will find all the data presented in the tables.

## Models

In addition to the table outputs, the scripts to run the models put forward in the report are included in the subdirectory models. These are very similar if not identical to the table scripts. This sub-directory is intended for experimenting with for those who are interested. Most of the scripts will have lines near the beginning were parameters such as number of classes are explicitly stated. These should be readily editable with the resultant model able to run. Should you find an issue occurs please do get in touch with me at jonvw28@gmail.com and I will investigate the cause of the error. This is the only section of this repository that is likely to change in time, though I may archive the state of this as it is now and create an updated sub-directory in due course.

## Acknowledgements

I would like to thank Dr Marco Catoni of the Sainsbury laboratory for giving me the idea this project is basedf upon and for his support through out.

I would also like to thank Hajk Drost for his constant advice and input in helping me to find direction for this project. 

Lastly I would like to thank Professor Jerzy Paskowski and the rest of the Paszkowski group for their useful feedback and support, and for being a very helpful practice audience for describing this research to biologists.

## Reference Information

The code presented here was prepared in R studio using R version 3.2.3 in a Windows architecture, with 64-bit operating system. The following packages and version were used:

* dplyr		0.4.3
* fitdistrplus	1.0-6
* GenomicRanges	1.22.3
* ggplot2	2.0.0
* plyr		1.8.3
* RColorBrewer	1.1-2
* rtracklayer	1.30.1
* stringr	1.0.0