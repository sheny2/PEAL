## PEAL: Privacy-preserving Efficient Aggregation of Longitudinal Data
This repository contains an R implementation of the PEAL (Privacy-preserving Efficient Aggregation of Longitudinal data) algorithm. PEAL is a novel, one-shot distributed algorithm designed to fit three-level linear mixed-effects models on longitudinal data without sharing individual patient data (IPD).

The key features of PEAL are:

* Privacy-Preserving: Operates exclusively on site-level summary statistics, ensuring patient-level data remains local. 

* Communication-Efficient: Requires only a single round of communication from participating sites. 

* Lossless: Achieves results that are statistically identical to a traditional pooled analysis where all IPD is centralized. 

This tutorial will guide you through simulating a multi-site dataset and using the PEAL engine to fit a model, then comparing the results to a standard linear mixed model fit with `lme4`.


### 

You'll need the PEAL engine script. Make sure the path to the engine file is correct.

```
three_lvl_dat <- readRDS("three_lvl_dat.rds")
head(three_lvl_dat)
```
