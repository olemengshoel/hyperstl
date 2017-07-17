HyperSTL
---

Accompanying code for paper "Optimizing the Decomposition of Time Series using Evolutionary Algorithms: Soil Moisture Analytics". Contains HyperSTL function to smooth the signal while preserving peaks and valleys.

# HyperSTL

## Dependencies

Luca Scrucca (2013). GA: A Package for Genetic Algorithms in R.
  Journal of Statistical Software, 53(4), 1-37. URL
  http://www.jstatsoft.org/v53/i04/.

Mullen, K.M, Ardia, D., Gil, D., Windover, D., Cline, J. (2011). DEoptim: An R Package for
    Global Optimization by Differential Evolution. Journal of Statistical Software, 40(6), 1-26. URL
    http://www.jstatsoft.org/v40/i06/.

Sylvain Gubian, Yang Xiang, Brian Suomela, Julia Hoeng, PMP SA (2016). GenSA: R Functions for Generalized Simulated Annealing
    https://cran.r-project.org/web/packages/GenSA/GenSA.pdf

## Usage

```
y.optimized.stl = HyperSTL(data, "value", 
                                algorithm = 'sga', 
                                max.feval = 1e2,
                                weights = c(5, 5, 8), 
                                data.freq = num.points.in.day, 
                                subsample.rate)

Params:
    #' @param **data** data.frame containing time series: soil moisture and rainfall. One column with name "rainfall" is required. 
    #' @param column name of data.frame-column containing soil moisture time series
    #' @param algorithm Optimization algorithm. Possible values are 
    #' 'sa', 'de', 'jade', 'sga', 'bga', 'irace'. Defaults to 'sa'.
    #' @param max.feval Max number of function evaluations. Defults to 2500.
    #' @param weights Weights (w1, w2, w3) for decomposed time series using STL. 
    #' Expects a vector containing three values. e.g. c(10,10,10) i.e. Equal weights for trend, seasonality and remainder.
    #' @param minparams Minimum params for search
    #' @param maxparams Maximum params for search
    #' @param data.freq Number of data points per season, e.g. season = 24 hrs, 
    #'                  data points are 2 minutes apart => data.freq = 24*60/2 = 720

Returns:
    STL object (object containing trend, seasonality, and remainder).
```

## Setup

* Clone this repository to a desired folder.
`$ cd <folder_of_your_choice>`
`$ git clone https://github.com/aniruddha55/HyperSTL.git`

## Run the code in R Studio

1. Set working directory to the folder where the repo has been cloned.
`> setwd("<folder_of_your_choice>/HyperSTL")`

2. Ensure working directory has been properly set.
`> getwd()`

3. Install necessary packages needed for HyperSTL code.
`> install.packages("GA")`
`> install.packages("DEoptim")`
`> install.packages("GenSA")`
`> install.packages("irace")`

4. Run given sample code that uses HyperSTL.
`> source('sampleRun.R')`
