##
## Optimization functions
##

## Properly reads a summary table into R
read.summary <- function(filename) {
  data <- read.table(filename,colClasses="character")
}

## Got some more functions to write here...


##
## Utility functions
##

## Creates a list of probabilities for evaluation later
probability.env <- function(probabilities,probability.names=sprintf("P%d",seq(1,length(probabilities)))) {
  data <- as.list(probabilities)
  names(data) <- probability.names
  return (data)
}

## Calculates probability table values, binding the provided variables
simplify.summary <- function(probability.table,rule.probabilities.env) {
  output <- as.data.frame(matrix(0,nrow(probability.table),ncol(probability.table)))
  names(output) <- names(probability.table)
  row.names(output) <- row.names(probability.table)
  for (n in seq(1,ncol(probability.table))) {
    for (m in seq(1,nrow(probability.table))) {
      output[m,n] <- eval(parse(text=probability.table[m,n]),envir=rule.probabilities.env)
    }
  }
  return (output)
}

