---
title: 'xrnet: Hierarchical Regularized Regression to Incorporate External Data'
authors:
- affiliation: 1
  name: Garrett M Weaver
  orcid: 0000-0002-9918-8386
- affiliation: 1
  name: Juan Pablo Lewinger
date: "31 July 2019"
output:
  html_document:
    df_print: paged
bibliography: paper.bib
tags:
- regularized regression
- lasso regression
- ridge regression
- elastic net regression
- hierarchical regression
affiliations:
- index: 1
  name: Department of Preventive Medicine, University of Southern California
---

# Summary

Regularized regression is a popular method for both variable selection and prediction in high-dimensional data. A number of R packages have been developed to fit regularized regression models in a generalized modelling framework with the option to use a variety of penalty types, including *glmnet* [@friedman2010], *biglasso* [@zeng2017], and *ncvreg* [@breheny2011]. The form of the constraints imposed by different penalties lead to different behavior, but all share one common attribute of shrinking estimates toward zero as regularization is increased. 

The *xrnet* R [@R] package extends the regularized regression framework to enable the inclusion of external data that may be informative for the mean effect of predictors on an outcome, an extension we have termed 'hierarchical regularized regression'. Incorporation of external data can improve predictive performance by allowing predictor effects to shrink towards potentially informative values other than zero. An additional regularization term on the external data enables variable selection and ensures there is little to no penalty when the external data is truly non informative.

Along with this extension, *xrnet* can also fit standard regularized regression models and integrates popular features from the R packages *glmnet* and *biglasso*. Below is a comparison of features that are available in *xrnet*, *glmnet*, *biglasso*:

| Feature | xrnet | glmnet | biglasso |
|---------|-------|--------|----------|
| Matrix types supported | Dense, Sparse, Memory-mapped | Dense, Sparse | Memory-mapped |
| Outcome types supported | Gaussian, Binomial | Gaussian, Multiresponse Gaussian, Binomial, Poisson, Cox | Gaussian, Binomial |
| Feature-specific penalty scaling | yes | yes | yes |
| User controls feature standardization | yes | yes | no |
| User controls inclusion of intercept | yes | yes | no |
| Box (upper/lower constrains on) estimates | yes | yes | no |
| Enhanced feature screening | no | no | yes |
| Integration of external data | yes | no | no |
| Feature-specific penalty types | yes | no | no |

To maintain computational efficiency, the core functionality was developed in C++ with usage of the Eigen linear algebra library [@eigenweb] for the primary data structures. Overall, this R package aims to provide a set of functions to fit and tune hierarchical regularized regression models and attempts to unify some of the best features from currently available R packages for regularized regression into a single easy to use interface.

# Funding and Support

This work is supported by the National Institute of Health (NIH) Grant P01CA196569.

# References
