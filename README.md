# PoissonRegression
Poisson regression to test correlates of substitution rates

Poireg.R contains an R function to estimate a regression model for substitution rates.

The inputs are:

1) formula: the formula of the regression model. For example, formula=expression(S~AMT+habitat) regresses synonymous substitution on annual mean temperature and habitat. 
2) data: the data matrix, where each row is a species pair and each column is a species trait with the same names as in the fomula. For example, S is used as synonymous substitution in the formula, then column named S1 is for the synonymous substitution in one species in a species pair and column named S2 is for the other species. Similarly, column named AMT1 is the annual mean temperature of one species and column named AMT2 is the annual mean temperature of the other species. column named habitat1 is the habitat of one species and column named habitat2 is the habitat of the other species. In the formula, habtat is a categorical variable, so data$habitat1 and data$habitat2 are factors with the same order of levels.
3) estimate_se: whether to estimate standard error of each coefficient.
4) init: the initial guess of the parameters in the regression model.

The outputs are:

1) mean: the ML estimates of each coefficient for each variable in the same order as in the formula.
2) se: the standard error of each coefficient, if estimate_se=T.
3) CI: the confidence interval of the coefficient, if estimate_se=T.
4) likratio: the likelihood ratio of the model versus the null model with no predictors.
5) pval: the p value of the likelihood test.
6) Rsquare: the pseudo R squared value of the regression model.
