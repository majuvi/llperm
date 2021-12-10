# llperm

LikeLihood based Permutation of Regression Residuals (PRR) test.

The basic approach is described in the paper (https://svn.r-project.org/Rjournal/trunk/html/_site/archive/2010-1/RJournal_2010-1_Werft+Benner.pdf):
> Werft, Wiebke, and Axel Benner. "glmperm: A Permutation of Regressor Residuals Test for Inference in Generalized Linear Models." R J. 2.1 (2010): 39.

The idea is to calculate a p-value for a variable of interest based on a permutation test. This approach is more robust than using likelihood based inference with a regression model.

We extend the original implementation in three ways:
1. Can be applied to any log-likelihood based model
2. A categorical variable of interest with multiple levels
3. Fit two models for count/zero component in zero-inflated models.

# Installation

`devtools::install_git("https://gitl01-int-p.rivm.nl/viljanem/llperm")`

# Usage

See the 'vignettes' folder.