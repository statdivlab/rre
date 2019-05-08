#' @title Maximum likelihood estimate for single frequency count table
#'
#' @description Finds the MLE for a single frequency count table with the option of several search schemes and penalization choices.
#'
#' @details THIS STILL NEEDS DETAILS, WORK IN PROGRESS.
#'
#' @param fct a frequency count table, see \code{\link{make_frequency_count_table}} for more information.
#' @param penalty method of penalizing/regularizing the maximum likelihood solution.
#' @param lambda parameter for the h1 penalization method.
#' @param search_scheme method for iterating over candidate C values.
#' @param multiplier Richness values from the observed richness to \code{multiplier}*(observed richness) will be tested (default: 10).
#' @param c_seq_len Number of candidate C values we will maximize over.
#' @param starts A dataframe or matrix of optimization starts.  See details.
#' @param alpha_min,delta_min,alpha_max,delta_max box constraints for optimization over \code{alpha} and \code{delta}.
#' @param forced_ccc_lower_bound allows the user to force a lower bound to the grid search for C.
#' @import foreach
#' @examples
#' direct_optimise(nb_fct_simulation(1000, 100, 0.95))
#'
#' # Suppose we want something quicker:
#' direct_optimise(nb_fct_simulation(1000, 100, 0.95),
#'     search_scheme = "bisection",
#'     c_seq_len = 20)
#'
#' @export
direct_optimise <- function(fct,
                            penalty = NULL, lambda = NULL,
                            search_scheme = "grid",
                            multiplier = 10, c_seq_len = 100, # should add a 'by' option
                            starts = NULL,
                            alpha_min = 1e-10, delta_min = 1e-10,
                            alpha_max = 1e5, delta_max = 1e5,
                            forced_ccc_lower_bound = NA){
  # if starts is not passed we assume the user wants just one reasonable start:
  if(is.null(starts)){
    starts <- t(as.matrix(c(alpha = 1e-2, delta = 1e-2))) #arbitrarily chosen
  } else if(is.data.frame(starts)){
    starts <- as.matrix(starts)
  } else if(!is.data.frame(starts) & !is.matrix(starts)) {
    stop("Starts must be a matrix or dataframe.")
  }

  input_starts <- starts

  ks <- fct[, 1]
  fs <- fct[, 2]
  cc <- sum(fs)

  outcome <- matrix(NA, nrow = 1, ncol = 4)
  outcome_allstarts <- matrix(NA, nrow = 1, ncol = 6)
  a <- 0 # a is just a counter for the row of outcome_allstarts
  b <- 0 # b is a counter for the row of outcome

  switch(search_scheme,
         "grid" = {
           environment(grid_search) <- environment()
           grid_search()
         },
         "bisection" = {
           environment(bisection_search) <- environment()
           bisection_search()
         },
         stop("The search scheme must be 'grid' or 'bisection'.")
  )

  outcome <- outcome[!is.na(outcome[,1]),]
  outcome_allstarts <- outcome_allstarts[!is.na(outcome_allstarts[,1]),] # trim NA

  colnames(outcome) <- c("ccc", "likelihood", "alpha", "delta")
  colnames(outcome_allstarts) <- c("ccc", "likelihood", "alpha", "delta", "a_start","d_start")
  # I'm leaving outcome as "full" to ensure backward compatibility, but its counter-intuitive now because it's not the fullest output.
  list("best" = outcome[which(outcome[,2] == max(outcome[,2], na.rm = T)), ], "full" = outcome, "full_starts" = outcome_allstarts)
}

#' @title helper function for direct_optimise
#' @description helper function for direct_optimise
#' @details See direct_optimise.
#' @keywords internal
grid_search <- function() {
  outcome <- matrix(NA, nrow = c_seq_len+1, ncol = 4)
  outcome_allstarts <- matrix(NA, nrow = c_seq_len*(nrow(starts)+2), ncol = 6)
  if (is.na(forced_ccc_lower_bound)) {
    c_vector = ceiling(seq(cc, multiplier*cc,length.out = c_seq_len))
  } else {
    c_vector = ceiling(seq(forced_ccc_lower_bound,
                           multiplier*forced_ccc_lower_bound,
                           length.out = c_seq_len))
  }

  for (i in 1:length(c_vector)) {
    b <- b+1

    # Adding the nearby (alpha, delta) optima if one exists
    if (i > 1) {
      likelihood_max_index <- (outcome[,2]==max(outcome[,2],na.rm=T)) %>%
        which %>% rev %>% '['(1)
      nearby_alpha_delta <- outcome[likelihood_max_index,3:4]
      starts <- rbind(input_starts,nearby_alpha_delta)
    }
    this_c_outcome <- matrix(NA, nrow = nrow(starts),ncol = 6)
    estimate <- c_vector[i]

    environment(optimise_fixed_c) <- environment()
    optimise_fixed_c() # brings this_c_outcome back to this environment.

    outcome_allstarts[(a+1):(a+nrow(starts)),] <- this_c_outcome
    outcome[b,] <- this_c_outcome[which(this_c_outcome[,2] == max(this_c_outcome[,2]))[1],1:4]
    a <- a+nrow(starts)
  }

  outcome <<- outcome
  outcome_allstarts <<- outcome_allstarts
}

#' @title helper function for direct_optimise
#' @description helper function for direct_optimise
#' @details See direct_optimise.
#' @keywords internal
bisection_search <- function() {
  outcome <- matrix(NA, nrow = c_seq_len+1, ncol = 4)
  outcome_allstarts <- matrix(NA, nrow = c_seq_len*(nrow(starts)+2), ncol = 6)
  c_init = ceiling(seq(cc,multiplier*cc,length=4))

  for (i in 1:length(c_init)) {
    b <- b+1
    starts <- input_starts # redundant
    estimate <- c_init[i]

    environment(optimise_fixed_c) <- environment()
    optimise_fixed_c() # brings this_c_outcome back to this environment.

    outcome_allstarts[(a+1):(a+nrow(starts)),] <- this_c_outcome
    outcome[b,] <- this_c_outcome[which(this_c_outcome[,2] == max(this_c_outcome[,2]))[1],1:4]
    a <- a+nrow(starts)
  }

  # stop flag is flipped to true when we have converged to our MLE, or if we have
  # tested more points than c_seq_len.
  stop_flag <- F
  # bisecting loop
  while(b < c_seq_len && !stop_flag) {
    outcome <- outcome[order(outcome[,1]),] # order by C
    currentML_ind <- which(outcome[,2] == max(outcome[,2],na.rm=T))
    # The following if statement is mostly for my own curiousity.  It can be deleted later, and the current MLE can just subset the first element to be safe.
    if (length(currentML_ind) > 1) {
      print("There were 2 or more points with the highest likelihood, namely...")
      print(currentML_ind)
      currentML_ind <- currentML_ind[1]
    }
    added_points <- numeric(0)
    if (sum(complete.cases(outcome)) == currentML_ind) { # case: the very last point tested has highest likelihood.
      added_points <- c(round((outcome[currentML_ind,1]+outcome[currentML_ind-1,1])/2))
    } else if(currentML_ind == 1) { # case: first slope is negative
      added_points <- c(round((outcome[currentML_ind,1]+outcome[currentML_ind+1,1])/2))
    } else {
      added_points <- c(round((outcome[currentML_ind,1]+outcome[currentML_ind+1,1])/2),
                        round((outcome[currentML_ind,1]+outcome[currentML_ind-1,1])/2))
    }

    if (any(!(added_points %in% outcome[,1]))) {
      added_points <- added_points[!(added_points %in% outcome[,1])]
    } else {
      stop_flag <- T # nothing left to test.
    }

    if (!stop_flag) {
      for (i in 1:length(added_points)) {
        b <- b+1
        estimate <- added_points[i]

        two_closest_C <- unique(outcome[which(abs(outcome[,1]-estimate) %in% sort(abs(outcome[,1] - estimate))[1:2]),1])
        low_C <- unique(outcome[which(outcome[,1]==min(two_closest_C)),3:4])
        high_C <- unique(outcome[which(outcome[,1]==max(two_closest_C)),3:4])
        starts <- rbind(input_starts,low_C,high_C)

        environment(optimise_fixed_c) <- environment()
        optimise_fixed_c()

        outcome_allstarts[(a+1):(a+nrow(starts)),] <- this_c_outcome
        outcome[b,] <- this_c_outcome[which(this_c_outcome[,2] == max(this_c_outcome[,2]))[1],1:4]
        a <- a+nrow(starts)
      }
    }
  }
  outcome <- outcome[order(outcome[,1]),] # one last time for consistency.

  outcome <<- outcome
  outcome_allstarts <<- outcome_allstarts
}


#' @title helper function called by grid_search and bisection_search
#' @description helper function, see direct_optimise.
#' @details See direct_optimise.
#' @keywords internal
optimise_fixed_c <- function() {
  this_c_outcome <- foreach(j=1:nrow(starts), .combine='rbind') %do% {
    x <- starts[j,]
    y <- optim(x, loglike_single_fixed_ccc,
               lower = list(alpha = alpha_min, delta = delta_min),
               upper = list(alpha = alpha_max, delta = delta_max), method = "L-BFGS-B",
               control = list(fnscale = -1, maxit = 100), ccc = estimate, cc = cc,
               ks = ks, fs = fs, penalty = penalty, lambda = lambda)
    c(estimate, y$value, y$par[1], y$par[2], x[1], x[2])
  }
  # stopCluster(cl)

  if (is.vector(this_c_outcome)) {
    this_c_outcome <- matrix(this_c_outcome, nrow=1)
  }
  this_c_outcome <<- this_c_outcome
}

#' @title Likelihood function written in R for reference
#' @description Likelihood function written in R for reference
#' @details Log likelihood function, completely equivalent to the C++ version which is actually used for computation.
#' @keywords internal
loglike_single_fixed_ccc <- function(x, ccc, cc, ks, fs, penalty, lambda) {
  l <- loglike_unpenalized_cpp(x,ccc,cc,ks,fs)
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
  l
}

#' @title Helper for likelihood.
#' @description Helper for likelihood.
#' @details Helper for likelihood.
#' @keywords internal
log_marginal_negative_binomial <- function(k, alpha, delta) {
  dnbinom(x = k, size = alpha, prob = delta/(delta + 1), log  = T)
}

