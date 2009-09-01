##
## Optimization functions
##

##
## Use this function for optmization.
##
## prob.env - a list of probabilities (with corresponding names) from probability.env()
## probability.table - a data frame of probability strings from read.probabilities()
## observed.table - a data frame of observed probabilities
probability.optmization <- function(prob.env,probability.table,observed.table) {
  result <- simplify.summary(probability.table,prob.env)
  return (sqrt(mean((result - observered.table)^2)))
}

##
## Properly reads a summary table into R
##
## filename - a generation summary file (eg. generation_001_summary.txt)
read.probabilities <- function(filename) {
  data <- read.table(filename,colClasses="character")
}

##
## Creates a list of probabilities for evaluation later
##
## probabilities - a vector of probability values (eg. c(0.1,0.4,...))
## probability.names - a vector of character strings (eg. c("P1","P2",...))
probability.env <- function(probabilities,probability.names=sprintf("P%d",seq(1,length(probabilities)))) {
  data <- as.list(probabilities)
  names(data) <- probability.names
  return (data)
}

##
## Calculates probability table values, binding the provided variables
##
## probability.table - a data frame of probability strings from read.probabilities()
## prob.env - a list of probabilities (with corresponding names) from probability.env()
simplify.summary <- function(probability.table,prob.env) {
  output <- as.data.frame(matrix(0,nrow(probability.table),ncol(probability.table)))
  names(output) <- names(probability.table)
  row.names(output) <- row.names(probability.table)
  for (n in seq(1,ncol(probability.table))) {
    for (m in seq(1,nrow(probability.table))) {
      output[m,n] <- eval(parse(text=probability.table[m,n]),envir=prob.env)
    }
  }
  return (output)
}

