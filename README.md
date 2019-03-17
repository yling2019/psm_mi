# Overview

This repository contains all code required to reproduce analyses in:

A. Ling, M. Mathur, K. Kapphahn, M. Montez-Rath, M. Desai. How to apply multiple imputation in propensity score matching with partially observed confounders? Under review, 2019+.

Address any inquiries to yling [AT] stanford [DOT] edu.

Step1.R
treatment effect estimate using true regression model using full data
- treatment effect estimate using PSM using full data
- treatment effect estimate using complete case analysis (CC), complete variable analysis (CVA), mean imputation, missing indicator

Step2.R
- various MI strategies
 
helperFunction.R
- data generation functions
- missing data generation functions
- implementation of various MI strategies

mice.r, internal.r, and is.r are modifications of source codes in MICE R package for implementation of MI-regPassive
