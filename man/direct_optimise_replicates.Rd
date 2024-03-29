% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/replicate_util.R
\name{direct_optimise_replicates}
\alias{direct_optimise_replicates}
\title{Direct maximum likelihood optimizer for a list of frequency count tables.}
\usage{
direct_optimise_replicates(fct_list, penalty = NULL, lambda = NULL,
  multiplier = 10, c_seq_len = 100, starts = NULL,
  alpha_min = 1e-10, delta_min = 1e-10, alpha_max = 1e+05,
  delta_max = 1e+05, search_scheme = "grid")
}
\arguments{
\item{fct_list}{a list of frequency count tables, see
\code{\link{make_frequency_count_table}} for more information on this data
structure.}

\item{penalty}{method of penalizing/regularizing the maximum likelihood solution.}

\item{lambda}{parameter for the h1 penalization method.}

\item{multiplier}{Richness values from the observed richness to \code{multiplier}*(observed richness) will be tested (default: 10).}

\item{c_seq_len}{Number of candidate C values we will maximize over.}

\item{starts}{A dataframe or matrix of optimization starts.  See details.}

\item{alpha_min}{box constraints for optimization over \code{alpha} and \code{delta}.}

\item{delta_min}{box constraints for optimization over \code{alpha} and \code{delta}.}

\item{alpha_max}{box constraints for optimization over \code{alpha} and \code{delta}.}

\item{delta_max}{box constraints for optimization over \code{alpha} and \code{delta}.}

\item{search_scheme}{method for iterating over candidate C values.}
}
\description{
Finds the MLE for a list of replicate frequency count tables
  drawn from the same population using the full likelihood.
}
\details{
\code{\link{direct_optimise_replicates()}} searches over candidate C
  values in a grid spanning from the maximum observed richness (c) to
  \code{multiplier} *c.  At each C value we search over \code{alpha},
  \code{delta} using a gradient search (\code{optim} with method =
  "L-BFGS-B").  We use multiple starts (given in the matrix \code{starts}).
  in optimizing as we have no guaruntee of convexity in this problem.
}
\examples{
list_of_fct <- rre::nb_fct_simulation(1000, 0.1, 0.1, 2) # two replicates
direct_optimise_replicates(list_of_fct)

}
