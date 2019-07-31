---
title: 'xrnet: R Package for Hierarchical Regularized Regression to Incorporate External Data'
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

Regularized regression is a popular method for both variable selection and prediction in high-dimensional data. A number of R packages have been developed to fit these models on various outcomes with a variety of penalty types, including **glmnet**, **biglasso**, and **ncvreg**. The **xrnet** R package extends the regularized regression framework to enable the inclusion of external data that may be informative for the effect of predictors on an outcome, an extension we have termed 'hierarchical regularized regression'. Incorporation of external data can improve predictive performance by allowing predictor effects to shrink towards potentially informative values other than zero. An additional regularization term on the external data enables variable selection and ensures there is little to no penalty when the external data is truly noninformative.

Along with the extesion described above, popular features from the R packages **glmnet** and **biglasso**, both common and unique, were integrated into **xrnet**. 

Features that are not found in both **glment** and **biglasso** or novel in this package, include:

* Integration of external data in a hierarchical regularization framework

* Support for standard R matrices, sparse matrices, and memory-mapped matrices

* Feature-specific penalty types and penalty scaling

* Reduced memory usage by using 'on the fly' variable standardization

* User-specified standardization of features

* Inclusion of unpenalized features as a separate matrix or in the primary feature matrix 

To maintain computational efficiency, the core functionality was developed in C++ with usage of the Eigen linear algebra library for the primary data structures. Overall, this R package aims to provide a set of functions to fit and tune a hierarchical regularized regression models and attempts to unify some of the best features from currently available R packages for regularized regression into a single easy to use interface.

# Funding and Support

This work is supported by the National Institute of Health (NIH) Grant P01CA196569.

# References
