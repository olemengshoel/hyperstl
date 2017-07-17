library('GenSA')
library('GA')
library('DEoptim')
library('irace')

#' @title HyperSTL
#' @description Optimization algorithm for STL
#' @param data data.frame containing time series: soil moisture and rainfall. One column with name "rainfall" is required. 
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
#'
#' @return Optimized version of time series.
#' @export
#'
#' @examples
HyperSTL <- function(data,
                     column,
                     algorithm = 'sa',
                     max.feval = 2500,
                     weights = c(20, 10, 30),
                     minparams = c(7, 7, 72),
                     maxparams = c(99, 99, 200),
                     data.freq = 24) {
  
  MinParams <- minparams
  MaxParams <- maxparams
  weights = weights / sum(weights)
  
  if (! 'rainfall' %in% colnames(data)) {
     stop('\n No column with name rainfall')
  }
  poly.smooth <- function(data, column) {
    y <- data[[column]]    
    intervals <- findInterval(data$rainfall, c(0, 1))
    index <- 1
    while (index < length(intervals)) {
      startIndex <- index
      
      while (index < length(intervals)
             && intervals[index] != 2) {
        index <- index + 1
      }
      
      ysplit <- y[startIndex:index]
      chunk <- data.frame(y = ysplit, x = 1:length(ysplit))
      
      model <- lm(y ~ stats:::poly(x, 1, raw = TRUE),
                  data = chunk)
      
      y[startIndex:index] <- fitted(model)
      index <- index + 1
      
    }
    y
  }
  
  y <- data[[column]]
  
  #  The "frequency" is the number of observations per season.
  #  Season = 24 hours
  #  For hourly data: frequency = 24
  y.ts <- ts(y, frequency = data.freq)
  yhat <- poly.smooth(data, column)

  getStl <- function(params) {
    y.stl <- stl(
      y.ts,
      s.window = params[1],
      t.window = params[2],
      l.window = params[3],
      robust = TRUE
    )
    y.stl
  }
  
  evalFunc <- function(params) {
    -evalFuncToMinimize(params)
  }
  
  evalFuncToMinimize <- function(params) {
    y.stl <- getStl(params)
    r <- as.numeric(y.stl$time.series[, "remainder"])
    y.trend <- as.numeric(y.stl$time.series[, "trend"])
    mse <- sqrt(mean((y.trend - yhat) ^ 2))
    
    objective <- sum(c(sd(r), abs(max(r) - min(r)), mse) * weights)
    if (is.nan(objective) ||
        is.na(objective) || is.null(objective)) {
      objective = 99999999
    }
    objective
  }
  
  evalFuncForBinary <- function(binary.params) {
    decoded.params =  decode.params.GA(binary.params)
    # penalty for out of range input
    OUT_OF_RANGE_PENALTY = 1e2
    
    
    params = decoded.params
    # # decoding min is set to MinParams
    # ind.small = which(decoded.params < MinParams)
    # params[ind.small] = MinParams[ind.small]
    ind.Big = which(decoded.params > MaxParams)
    params[ind.Big] = MaxParams[ind.Big]
    
    y.stl <- getStl(params)
    r <- as.numeric(y.stl$time.series[, "remainder"])
    y.trend <- as.numeric(y.stl$time.series[, "trend"])
    mse <- sqrt(mean((y.trend - yhat) ^ 2))
    objective <- sum(c(sd(r), abs(max(r) - min(r)), mse) * weights)
    
    objective <- objective
    + OUT_OF_RANGE_PENALTY * norm(params - decoded.params, type = "2") ^
      2
    
    - objective
  }
  
  switch(
    algorithm,
    sa = {
      # Simulated Annealing
      out <-
        GenSA(
          lower = MinParams,
          upper = MaxParams,
          fn = evalFuncToMinimize,
          control = list(max.call = max.feval, verbose = TRUE)
        )
      print(out$par)
      print(cat("Evalfunc(GenSA) = ", evalFuncToMinimize(out$par)))
      solution = out$par
      solution.obj = evalFuncToMinimize(solution)
    },
    # ------------------------------------------------------------------
    de = {
      # DE/rand/1/bin
      itermax = max(1, floor(max.feval / 50))
      outDEoptim <-
        DEoptim(
          evalFuncToMinimize,
          MinParams,
          MaxParams,
          DEoptim.control(
            NP = 50,
            itermax = itermax,
            F = 0.8
          )
        )
      solution = outDEoptim$optim$bestmem
      solution.obj = evalFuncToMinimize(solution)
    },
    # ------------------------------------------------------------------
    jade = {
      # JADE
      itermax = max(1, floor(max.feval / 50))
      outDEoptim <-
        DEoptim(
          evalFuncToMinimize,
          MinParams,
          MaxParams,
          DEoptim.control(
            NP = 50,
            itermax = itermax,
            F = 0.8,
            CR = 0.5,
            strategy = 6,
            c = 0.4
          )
        )
      solution = outDEoptim$optim$bestmem
      solution.obj = evalFuncToMinimize(solution)
      # print(evalFuncToMinimize(solution))
    },
    # ------------------------------------------------------------------
    sga = {
      # Real coded GA
      GAModel <- ga(
        type = "real-valued",
        fitness = evalFunc,
        monitor = NULL,
        min = MinParams,
        max = MaxParams,
        maxfitness = max.feval
      )
      solution = GAModel@solution
      solution.obj = -evalFunc(GAModel@solution)
      
    },
    # ------------------------------------------------------------------
    bga = {
      # Binary-GA
      GAModel <- ga(
        type = "binary",
        fitness = evalFuncForBinary,
        nBits = 21,
        maxfitness = max.feval,
        monitor = NULL
      )
      solution.binary = GAModel@solution
      solution = decode.params.GA(solution.binary)
      ind.Big = which(solution > MaxParams)
      solution[ind.Big] = MaxParams[ind.Big]
      
      solution.obj = -1 * evalFuncForBinary(solution.binary)      
      print(solution.obj)      
    },
    # ------------------------------------------------------------------
    irace = {
      # iRace
      hook.run = function(experiment, config = list()) {
        candidate <- experiment$candidate
        x = c(candidate[["x1"]], candidate[["x2"]], candidate[["x3"]])
        y = evalFuncToMinimize(x)
        y
      }
      # parameters.table <- 'x1 \"\" i \n (7, 99) x2 \"\" i (7, 99) \n x3 \"\" i (72, 200)'
      parameters.table <-
        'x1 "" i (7, 99)\n x2 "" i (7, 99)\n x3 "" i (72, 200)'
      parameters <- readParameters(text = parameters.table)
      print('Read parameters')
      irace.weights <- rnorm(200, mean = 0.9, sd = 0.02)
      result <- irace(
        tunerConfig = list(
          hookRun = hook.run,
          instances = irace.weights[1:100],
          maxExperiments = 4000,
          # nbIterations = 1,
          logFile = ""
        ),
        parameters = parameters
      )
      
      solution = c(result[1, ]$x1, result[1, ]$x2, result[1, ]$x3)
      solution.obj = evalFuncToMinimize(solution)
      print(evalFuncToMinimize(solution))
    },
    # ------------------------------------------------------------------
    {
      print('No algorithm specified')
      solution = c(0, 0, 0)
      solution.obj = 12345.12345
    }
  )
  stl.result <- getStl(solution)  
  stl.result
}
