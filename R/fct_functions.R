#' Creates a frequency count table
#'
#' \code{make_frequency_count_table()} takes a vector of species abundances and
#' returns a frequency count table.
#'
#' A frequency count table is a data frame listing the frequency ('freq') and
#' the number of species observed that many times in the sample ('numTaxa').
#'
#' In order to allow this function to be used easily with simulated data, we
#' omit any zero observations from our frequency count table.
#'
#' If \code{cautious = T}, in the event that all of the species abundances
#' listed are zero, the function returns a data frame with the correct labels
#' for a frequency count table but zero rows.
#'
#' @param abundances A numeric vector of absolute species abundances.
#' @param cautious logical.  If TRUE (default) then a vector of all zeroes will
#'   return a data frame.  If FALSE then an error will be thrown.
#' @return returns a data.frame.
#' @examples
#' taxa <- rnbinom(1000,0.1,0.1) # simulated data
#' make_frequency_count_table(taxa)
#'
#' \dontrun{
#' # The following will produce an error:
#' make_frequency_count_table(rep(0,100),cautious= F)
#' }
#' @export
make_frequency_count_table <- function(abundances, cautious = T) {
  if(cautious & all(abundances == 0)) {
    x <- data.frame(freq = numeric(0),
                    numTaxa = numeric(0))
  } else {
    x <- as.data.frame(table(abundances[abundances != 0]))
    x[, 1] <- as.numeric(as.character(x[, 1]))
    x[, 2] <- as.numeric(as.character(x[, 2]))
    names(x) <- c("freq","numTaxa") # check this
  }
  return(x)
}

#' @title Get formatted frequency count table data
#'
#' @description Formats the frequency count table data in the format
#'   required by the liklihood function
#'
#' @details I may hide this from the user in the final version of this package.  It is primarily useful for running tests on the likelihood function.
#'
#' @param fct A frequency count table.
#' @return This function returns a list with cc (observed richness), and vectors which are just the columns of the frequency count table.
#' @examples
#' toy_fct <- data.frame(freq = c(1,2,3), numTaxa = c(14,7,2))
#' get_cc_ks_fs(toy_fct)
#' @export
get_cc_ks_fs <- function(fct) {
  ks <- fct[,1]
  fs <- fct[,2]
  cc <- sum(fs)
  return(list(cc,ks,fs))
}

#' @title Negative binomial frequency count table simulation
#'
#' @description Simulates a frequency count table according to the negative
#'   binomial model with the parameters input.
#'
#' @details Assuming a gamma mixed Poisson model, the number of times any
#'   taxa is observed follows a single draw from the negative binomial
#'   distribution.  To simulate for C species, we make C draws.
#'
#' @param C The true species richness
#' @param alpha Shape parameter in the gamma distribution (our assumed
#' mixing distribution)
#' @param delta Rate parameter in the gamma distribution (our assumed
#' @param r Replicates.  If r is 1 then a single frequency count table is returned, otherwise a list of frequency count tables is returned.
#' mixing distribution)
#' @examples
#' nb_fct_simulation(1000, 0.1, 0.001) #output is a frequency count table
#' @return Returns a data frame in the form of a frequency count table
#' @export
nb_fct_simulation <- function(ccc, alpha, delta, r = 1) {
  if (r > 1 & r%%1 == 0) {
    fct_list <- rnbinom(n = ccc*r,
                        size = alpha,
                        prob = (delta/(1+delta))) %>%
      matrix(. ,nrow = ccc, ncol = r) %>%
      split(.,col(.)) %>%
      lapply(make_frequency_count_table)
    return(fct_list)
  } else {
    fct <- rnbinom(n = ccc,
                   size = alpha,
                   prob = delta/(1+delta)) %>%
      make_frequency_count_table
    return(fct)
  }
}

#' @title Creates a vector of observations from a frequency count table
#'
#' @description Takes a frequency count table and outputs the observations for
#'   each speices observed at least once in the data set.
#'
#' @details Species observed zero times cannot be included in the output of this
#'   function because we have established a convention that a FCT does not
#'   contain that information.  This function is the inverse of
#'   \code{\link{make_frequency_count_table()}} if all species are observed
#'   at least once.
#' @param fct A data frame in the form of a frequency count table.
#' @examples
#' toy_fct <- data.frame(freq = c(1,2,3), numTaxa = c(14,7,2))
#' abun <- fct_to_observations(toy_fct) # abundances for observed taxa
#' fct <- make_frequency_count_table(abun) # back to where we started
#' @seealso \code{\link{make_frequency_count_table()}}
#' @export
fct_to_observations <- function(fct) {
  obs <- c()
  for (k in fct[,1]) {
    reps <- fct[which(fct[,1] == k),2]
    obs <- c(obs,rep(k,times = reps))
  }
  obs
}

#' @title Expected value and variance of the zero truncated negative binomial
#'
#' @description Functions for calculating the mean or variance of the zero truncated negative binomial
#'   distribution (two parametrizations available).
#'
#' @details If \code{alpha} and \code{p} are passed, then we follow the
#'   parameterization used in \code{\link[stats]{rnbinom}}.  If \code{alpha} and
#'   \code{delta} are passed, then the distribution follows the parametrization
#'   derived by assuming the number of observations for each species is a
#'   gamma-mixed Poisson.  In this derivation the abundances are
#'   \eqn{Gamma(\alpha, \delta)} in the shape/rate parametrization.
#'
#'   These parametrizations are related by \eqn{p = \delta  / (1+\delta)}.
#'
#'   If \code{p} and \code{delta} are both passed then an error is thrown.
#'
#' @param alpha The \code{size} parameter fo \code{\link[stats]{rnbinom}} (also
#'   the shape parameter in the gamma mixing distribution, see 'Details')
#' @param p The \code{prob} parameter of \code{\link[stats]{rnbinom}}.
#' @param delta An alternative parametrization to using p.  (see 'Details')
#' @examples
#' ztnb_ev(0.1,0.5) # using p
#' identical(ztnb_ev(0.1, delta = 1), ztnb_ev(0.1,0.5)) # equivalent.
#'
#' ztnb_var(100,0.95)
#' @return A length 1 numeric vector.
#' @export
ztnb_ev <- function(alpha, p = NULL, delta = NULL) {
  if (is.null(delta)){
    if (is.null(p)) {
      stop("You must pass either delta or p.")
    }
    return((alpha*(1-p))/(p*(1-p^alpha)))
  } else {
    if (!is.null(p)) {
      stop("You cannot pass both delta and p")
    }
    return(alpha/(delta*(1-(delta/(1+delta))^alpha)))
  }
}

# testthat function for later:
#ztnb_ev(alpha = 0.1)
#ztnb_ev(alpha = 0.1, p = 0.1, delta = 55)
#ztnb_ev(alpha = 0.1, p = 0.6)
#rnbinom(10^6, 0.1, 0.6) %>% '['(. > 0) %>% mean
#ztnb_ev(alpha = 0.1, delta = 0.234)
#rnbinom(10^6, 0.1, (0.234)/(1+0.234)) %>% '['(. > 0) %>% mean

#' @rdname ztnb_ev
#' @export
ztnb_var <- function(alpha, p = NULL, delta = NULL) {
  if (is.null(delta)){
    if (is.null(p)) {
      stop("You must pass either delta or p.")
    }
    term_1 <- (alpha^2*(1-p)^2*p^alpha)/((1-p^alpha)*(p^(alpha+2)-p^2))
    term_2 <- alpha*(1-p)/(p^2*(1-p^alpha))
    return(term_1+term_2)
  } else {
    if (!is.null(p)) {
      stop("You cannot pass both delta and p")
    }
    scalar <- 1/(1-(delta/(1+delta))^alpha)
    term_1 <- alpha^2*delta^(alpha-2)*(1/(1+delta))^alpha/((delta/(1+delta))^alpha-1)
    term_2 <- alpha*(1+delta)/(delta^2)
    return(scalar*(term_1+term_2))
  }
}

# testthat for later:

#ztnb_var(alpha = 0.1)
#ztnb_var(alpha = 0.1, p = 0.1, delta = 55)
#ztnb_var(alpha = 0.1, p = 0.6)
#rnbinom(10^6, 0.1, 0.6) %>% '['(. > 0) %>% var
#ztnb_var(alpha = 0.1, delta = 0.234)
#rnbinom(10^6, 0.1, (0.234)/(1+0.234)) %>% '['(. > 0) %>% var
#' @rdname chi_sq_gof
#' @export
chi_sq_gof <- function(fct, ccc_hat, alpha_hat, delta_hat) {
  ks <- fct[,1]
  kmax <- max(fct[,1])
  f_k <- rep(0, kmax)
  f_k[ks] <- fct[, 2]

  p_eta_hat <- dnbinom(1:kmax, alpha_hat, delta_hat/(1+delta_hat))
  f_k_hat <- ccc_hat*p_eta_hat
  finite_sum <- sum((f_k - f_k_hat)^2/f_k_hat)

  p_eta_tail <- pnbinom(kmax+1, alpha_hat, delta_hat/(1+delta_hat), lower.tail = F)
  f_k_tail <- p_eta_tail*ccc_hat

  return(finite_sum+f_k_tail)
}


# Adds the zero rows into a frequency count table.
#' @rdname ztnb_ev
#' @export
expand_fct <- function(fct) {
  k_vals <- fct[,1]
  k_max <- max(k_vals)
  expanded_fct <- matrix(1:k_max,
                         rep(0, times = k_max))
  expanded_fct[k_vals,2] <- fct[,2]
  return(expanded_fct)

}




