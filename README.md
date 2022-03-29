# llperm

LikeLihood based Permutation of Regression Residuals (PRR) test.

The basic approach is described in the paper (https://svn.r-project.org/Rjournal/trunk/html/_site/archive/2010-1/RJournal_2010-1_Werft+Benner.pdf):
> Werft, Wiebke, and Axel Benner. "glmperm: A Permutation of Regressor Residuals Test for Inference in Generalized Linear Models." R J. 2.1 (2010): 39.

The idea is to calculate a p-value for a variable of interest based on a permutation test. This approach is more robust than using likelihood based inference with a regression model.

We extend the original implementation in three ways:
1. A categorical variable of interest with multiple levels
2. Can be applied to any log-likelihood based model, including over-dispersed and zero-inflated models.
3. Fit both count/zero components in zero-inflated models, then test either or both components separately.

# Installation

Simple:
`devtools::install_git("https://gitl01-int-p.rivm.nl/viljanem/llperm")`

If you cannot install the package, it is possible to just download the .zip and source all files:
`sapply(file.path("llperm", "R", list.files(file.path("llperm", "R"), pattern="*.R")),source,.GlobalEnv)`

# Usage

See the 'vignettes' folder. We have two vignettes:
1. 'introduction.Rmd' shows how to fit the models to a real data set (VEGA subset), and shows that glmperm can fix standard models.
2. 'introduction2.Rmd' shows how to simulate Zero-Inflated Negative Binomial Data, and demonstrates that the standard models do well.

Experiments for the paper draft can be found in 'llperm_experiments.Rmd'.

# License

GPL-3. The package was inspired by the GPL licensed packages glmperm, pscl, ZIBBSeqDiscovery and glm.fit in base R.