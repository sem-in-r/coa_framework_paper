# COA Framework Paper Code and Data

This folder contains code and data for Danks, Ray, and Shmueli (in press) paper title: "The Composite Overfit Analysis Framework: Assessing the Out-of-sample Generalizability of Construct-based Models Using Predictive Deviance, Deviance Trees, and Unstable Paths". For up-to-date version of this folder, please see its Github repository: https://github.com/sem-in-r/coa_framework_paper

The code in this folder should be used to *reproduce* the tables and figures of the manuscript and appendices. It utilize a library for the COA framework described in the paper. The code for COA is found in a separate folder. The COA library is generic and can be *reused* with other models and data. To use the most up-to-date SEMCOA framework code in your own project, please refer to its Github repository: https://github.com/sem-in-r/semcoa

## Files

- `COA_Framework_MS_Paper.R`: Main script that produces the result tables and figures (see Code section below)
- `LICENSE`: License covering use of code in this folder
- `README.md`: This folder/file description file
- `coa_framework_paper.Rproj`: RStudio project file for RStudio users
- `data_dictionary.md`: Data dictionary describing columns in main data file
- `returnlist2_11082021.rda`: Cached results of certain step(s) reused in to run the main code faster (see comments in main code file)
- `utaut2_216.csv`: Main data file of survey responses (see description in data dictionary file)

## Code

There is only a single code file (see file list) that has all the code to produce tables and figures of results in manuscript and appendices. Please see comments in the code file to see which code produces which results.

Note that results are produced by both a simulation and from actual survey data. Simulation code is described in comments in the code file. Survey data is described below in the Data section.

## Data

The empirical demonstration in the main manuscript uses a self-conducted survey. The description of survey collection procedure, including eligibilty requirements, are in the main manuscript. A main data file (see file list) contains all the survey responses. Please see data dictionary file for further description of the columns in the data file.
