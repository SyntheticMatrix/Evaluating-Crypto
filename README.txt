#########################
Instruction of the supplementary files:

#######################
There are four mod.files to replicate our baseline results, where
Cryptocurrency_Full_Anonymity_baseline.mod generates the Impulse Responses of traditional Money demand shock, cryptocurrency demand shock, and Exchange rate depreciation shock under full anonymity cryptocurrency case.

Cryptocurrency_Partial_Anonymitiy_baseline.mod generates the Impulse Responses under partial anonymity case.

Cryptocurrency_Full_Anonymity_Optimal_Monetary_Policy.mod shows how the taylor rules achieve their Targets under full anonymity case. One has to mannually switch the taylor rule between TITO and AIAO, by manipulating the % sign before the certain welfare loss function.

Cryptocurrency_Partial_Anonymity_Optimal_Monetary_Policy.mod shows how the Taylor rules achieve their Targets under partial anonymity case. One has to do the same as Above to switch between TITO and AIAO policies.

NOTICE THAT we run the mod files using Dynare 6.1, so we suggest the same version, otherwise different results might appear due to the updated algorithm. ALSO, the file Name may appear to be too long, so one can mannually rename the mod file as what serves his/her convenience.

ccestrrnew.m file is the cooked data by applying Pfeifer (2014) in which the data should be detrended to fit Bayesian estimation.

#########################
There are also three m.files, which will replicate the Impulse Response figures in the manuscript, which are the same as the mod.files yield, but with better visualization.

#########################
Appendix.tex, cc.bib, and all the figures.png are the latex supported files. Appendix should be placed directly following the manuscript, cc.bib should be put into the same root folder as the manuscript,  in order to compile the latex file with completion.

#########################
We hope These serve well. Any Question or error encountering, please contact our corresponding author.