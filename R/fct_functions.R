#' Creates a frequency count table
#'
#' \code{make_frequency_count_table()} takes a vector of species abundances and
#' returns a frequency count table.
#'
#' A frequency count table is a data frame listing the frequencies ('freq') and
#' the number of species observed that many times in the sample ('numTaxa').
#'
#' If \code{cautious = T}, in the event that all of the species abundances
#' listed are zero, the function returns a data frame with the correct labels
#' for a frequency count table but zero rows.  This is a precaution mostly
#' designed for small scale simulations where zero species may be observed - it
#' is unlikely to be relevant in real data.
#'
#' @param abundances A numeric vector of absolute species abundances (counts of
#' how many times the species was observed).
#' @param cautious logical.  If TRUE (default) then a vector of all zeroes will
#'   return a data frame with zero rows.  If FALSE then an error will be thrown.
#' @return returns a data.frame.
#' @examples
#' taxa <- rnbinom(1000,0.1,0.1) # simulated abundances
#' make_frequency_count_table(taxa) # frequency count table
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


#' @rdname chi_sq_gof
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


#' @rdname ztnb_ev
expand_fct <- function(fct) {
  k_vals <- fct[,1]
  k_max <- max(k_vals)
  expanded_fct <- matrix(1:k_max,
                         rep(0, times = k_max))
  expanded_fct[k_vals,2] <- fct[,2]
  return(expanded_fct)

}

#' @keywords internal
format_fct <- function(fct) {
  ks <- fct[,1]
  fs <- fct[,2]
  cc <- sum(fs)
  return(list(cc=cc,ks=ks,fs=fs))
}

#' @keywords internal
make_formatted_lists <- function(list_of_fct) {
  len <- length(list_of_fct)
  cc_list <- list()
  fs_list <- list()
  ks_list <- list()
  for (z in 1:len) {
    temp <- format_fct(list_of_fct[[z]])
    cc_list[[z]] <- temp$cc
    fs_list[[z]] <- temp$fs
    ks_list[[z]] <- temp$ks
  }
  stopifnot(length(cc_list) == length(fs_list),
            length(fs_list) == length(ks_list))
  return(list(cc_list = cc_list, fs_list = fs_list, ks_list = ks_list))
}

#' @title Returns the maximum observed richness for a FCT list.
#'
#' @keywords internal
get_cc_max <- function(list_of_fct) {
  if (is.data.frame(list_of_fct) | is.matrix(list_of_fct)) {
    list_of_fct <- list(list_of_fct)
  }
  lapply(list_of_fct, (function(x) x[,2] %>% sum)) %>% unlist %>% max
}





