#' @title Direct maximum likelihood optimizer for a list of frequency count tables.
#'
#' @description Finds the MLE for a list of replicate frequency count tables
#'   drawn from the same population using the full likelihood.
#'
#' @details \code{\link{direct_optimise_replicates()}} searches over candidate C
#'   values in a grid spanning from the maximum observed richness (c) to
#'   \code{multiplier} *c.  At each C value we search over \code{alpha},
#'   \code{delta} using a gradient search (\code{optim} with method =
#'   "L-BFGS-B").  We use multiple starts (given in the matrix \code{starts}).
#'   in optimizing as we have no guaruntee of convexity in this problem.
#'
#' @inheritParams direct_optimise
#' @param fct_list a list of frequency count tables, see
#'   \code{\link{make_frequency_count_table}} for more information on this data
#'   structure.
#' @import foreach
#' @import magrittr
#' @examples
#' list_of_fct <- rre::nb_fct_simulation(1000, 0.1, 0.1, 2) # two replicates
#' direct_optimise_replicates(list_of_fct)
#'
#' @export
direct_optimise_replicates <- function(fct_list,
                                       penalty = NULL, lambda = NULL,
                                       multiplier = 10, c_seq_len = 100,
                                       starts = NULL,
                                       alpha_min = 1e-10, delta_min = 1e-10,
                                       alpha_max = 1e5, delta_max = 1e5,
                                       # search_scheme is always grid, arg for backward compatibility.
                                       search_scheme = "grid") {
  if (is.null(starts)) {
    starts <- t(as.matrix(c(xalpha = 1e-2, delta = 1e-2))) #arbitrarily chosen
  } else if (is.data.frame(starts)) {
    starts <- as.matrix(starts)
  } else if(!is.data.frame(starts) & !is.matrix(starts)) {
    stop("Starts must be a matrix or dataframe.")
  }

  input_starts <- starts

  dat_list_formatted <- make_formatted_lists(fct_list)
  cc_list <- dat_list_formatted$cc_list
  fs_list <- dat_list_formatted$fs_list
  ks_list <- dat_list_formatted$ks_list

  cc_max <- cc_list %>% unlist %>% max
  #print(paste0("cc_max is ", cc_max))

  c_vector = ceiling(seq(cc_max, multiplier*cc_max,length.out = c_seq_len))

  # As we iterate through ccc values we will use the current maximum likelihood
  # alpha and delta as an "extra" start at each ccc
  extra_alpha <- NA
  extra_delta <- NA
  current_max_loglike <- NA

  # gradient_optimize takes a fixed value of ccc and optimizes
  # over eta = (alpha, delta), starting at each start.
  gradient_optimize <- function(ccc_est) {
    if (!is.na(extra_alpha)) { # first iteration only.
      starts <- rbind(input_starts, c(extra_alpha, extra_delta))
    } else {
      starts <- input_starts
    }

    ccc_df <- foreach(j=1:nrow(starts), .combine='rbind') %do% {
      x <- starts[j,]
      y <- optim(x, replicate_likelihood,
                 lower = list(alpha = alpha_min, delta = delta_min),
                 upper = list(alpha = alpha_max, delta = delta_max), method = "L-BFGS-B",
                 control = list(fnscale = -1, maxit = 100),
                 ccc = ccc_est,
                 cc_list = cc_list,
                 ks_list = ks_list,
                 fs_list = fs_list, penalty = penalty, lambda = lambda)
      c(ccc_est, y$value, y$par[1], y$par[2], x[1], x[2])
    }
    if (is.vector(ccc_df)) { # if there was one start.
      ccc_df <- matrix(ccc_df, nrow=1)
    }

    ml_vec <- ccc_df[which(ccc_df[,2] == max(ccc_df[,2]))[1], ]
    # replace the extra starts if needed.
    if (is.na(current_max_loglike) | ml_vec[2] > current_max_loglike) {
      extra_alpha <<- ml_vec[3]
      extra_delta <<- ml_vec[4]
    }
    return(ccc_df)
  }

  full_starts <- lapply(c_vector, gradient_optimize) %>%
    do.call('rbind', .) %>%
    as.data.frame
  colnames(full_starts) <- c("ccc", "likelihood", "alpha", "delta",
                             "a_start","d_start")
  full <- full_starts %>%
    dplyr::group_by(ccc) %>%
    dplyr::arrange(desc(likelihood)) %>%
    dplyr::filter(1:n() == 1) %>%
    dplyr::ungroup()

  best <- full %>%
    dplyr::arrange(desc(likelihood)) %>%
    dplyr::filter(1:n() == 1) %>%
    dplyr::ungroup()

  return(list("best" = best,
              "full" = full,
              "full_starts" = full_starts))
}

#' @title gamma-Poisson log likelihood for replicates
#'
#' @description The gamma-Poisson likelihood (which has parameters \eqn{C, \eta}) for a sample of frequency count tables drawn from the same population.
#'
#' @details The paramaters \code{cc_list}, \code{ks_list}, \code{fs_list} are a
#'   formatted version of the frequency count table list.  The likelihood
#'   function is called many times as we optimise, so we require them to be
#'   formatted this way to save time in the long run.  This function uses RCPP
#'   to speed up the unavoidable factorial calculations in the likelihood.
#' @param x a length 2 vector, \eqn{eta = (\alpha, \delta)}.
#' @param ccc The value of \eqn{C} to evaluate the likelihood at.
#' @param cc_list A list, each item is a length 1 numeric vector stating the observed richness for each replicate.
#' @param ks_list A list, each list item is a vector for that replicate.  The
#'   vector represents the \eqn{k} values for the nonzero frequencies.
#' @param fs_list A list, each list item is a vector for that replicate.  The
#'   vector represents the frequencies \eqn{f_k} associated with the \eqn{k}
#'   values in \code{ks_list}.
#' @param penalty The penalty function form to use, only "h1" is supported currently.
#' @param lambda The penalization parameter.
#' @import foreach
#' @import magrittr
#' @examples
#' list_of_fct <- rre::nb_fct_simulation(1000, 0.1, 0.1, 2) # two replicates
#' direct_optimise_replicates(list_of_fct)
#'
#' @keywords internal
replicate_likelihood <- function(x, ccc, cc_list,
                                 ks_list, fs_list,
                                 penalty, lambda) {
  alpha <- x[1]
  delta <- x[2]
  l <- Map(loglike_unpenalized_cpp, cc_list, ks_list, fs_list,
           MoreArgs = list(x = x, ccc = ccc)) %>%
    unlist %>% sum
  l_p0 <- log_marginal_negative_binomial(0, x[1], x[2])
  if(is.null(penalty)) {
    return(l) # unpenalized return
  } else if(penalty == "h1") {
    if (is.null(lambda)) {
      stop("The h1 penalty requires a lambda value to be set")
    }
    return(l-lambda*l_p0) #h1 penalty
  } else {
    stop("Invalid penalty parameter passed")
  }
}
