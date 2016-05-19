##
## 
##                This source code is part of
## 
##                          Cell-Diff
## 
##       Exact computation of stochastic stem cell models
## 
##                        VERSION 1.0.0
## Written by Joshua L. Phillips.
## Copyright (c) 2010-2016, Joshua L. Phillips.
## Check out http://www.cs.mtsu.edu/~jphillips/software.html for more
## information.
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
## 
## If you want to redistribute modifications, please consider that
## derived work must not be called official Cell-Diff. Details are found
## in the README & LICENSE files - if they are missing, get the
## official version at github.com/douradopalmares/mdsctk/.
## 
## To help us fund Cell-Diff development, we humbly ask that you cite
## the poster on the package - you can find them in the top README file.
## 
## For more info, check our website at
## http://www.cs.mtsu.edu/~jphillips/software.html
## 
##

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
  return (sqrt(mean((result - observed.table)^2)))
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

