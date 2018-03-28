resultsSeeds <- setClass(
  'resultsSeeds',
  slots = c(
    stateEstimates = "matrix",
    stateUnscertainLower = "matrix",
    stateUnscertainUpper = "matrix",
    hiddenInputEstimates = "matrix",
    hiddenInputUncertainLower = "matrix",
    hiddenInputUncertainUpper = "matrix",
    outputEstimates = "matrix",
    Data = "matrix",
    DataError = "matrix"
  ),
  prototype = c(
    stateEstimates = matrix(0),
    stateUnscertainLower = matrix(0),
    stateUnscertainUpper = matrix(0),
    hiddenInputEstimates = matrix(0),
    hiddenInputUncertainLower = matrix(0),
    hiddenInputUncertainUpper = matrix(0),
    outputEstimates = matrix(0),
    Data = matrix(0),
    DataError = matrix(0)
  )
)

