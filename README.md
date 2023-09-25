## README ##

Biased Statistical Error estimates under resampling with replacement.

Thanks for reading this!

Here's a short overview over the files in the OSF repository.

MATLAB is a folder containing matlab code used for exploratory tests of Marek et al.'s (2022) code. This work would not have been possible without the availability of their code, for which we are grateful.
The project was since moved to the R environment for open-source reproducibility and figures have been included for reproducibility purposes as well.

Please contact charlesdgburns@gmail.com for any issues with reproducibility.

 ## MATLAB ##

Before running 'TestingStatErrors.m', please download the following and add it to your matlab directory:

'abcd_edgewise_correlation_iterative_reliability_single_factor.m' and
'abcd_statisticalerrors.m' 

from the Marek et al. (2022) BWAS code repository:

https://gitlab.com/DosenbachGreene/bwas


 ## R ##

Analyses have been coded in separate R markdown files. 

1) "subject_level.Rmd"" (Figures 1, 2 and Supplementary figures)

2) "correlation_level.Rmd" (Figure 3)

3) "Fig2SrcSim.Rmd" (Supplementary Figure 1)

We recommend navigating through the Chunks, as these have been separated for different tasks, such as generating data, resampling with replacement, arranging plots, et.c.

The workflow in R was enabled thanks to the following packages:

Kassambara, A. (2020). ggpubr:'ggplot2'Based Publication Ready Plots. R package version 0.4. 0. Computer software]. https://cran-r-project. org/web/packages/ggpubr/indes. html.

Ripley, B., Venables, B., Bates, D. M., Hornik, K., Gebhardt, A., Firth, D., & Ripley, M. B. (2013). Package ‘mass’. Cran r, 538, 113-120.

Wickham, H., Averick, M., Bryan, J., Chang, W., McGowan, L. D. A., François, R., ... & Yutani, H. (2019). Welcome to the Tidyverse. Journal of open source software, 4(43), 1686.
