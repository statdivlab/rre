#' @title Maximum likelihood estimate for a list of frequency count tables.
#'
#' @description Finds the MLE for a list of replicate frequency count tables
#'   (i.e. frequency count tables from the same ecosystem) with the option of
#'   several search schemes and penalization choices.
#'
#' @details THIS STILL NEEDS MUCHO DETAILS, WORK IN PROGRESS.
#' @inheritParams direct_optimise
#' @param fct_list a list of frequency count tables, see \code{\link{make_frequency_count_table}} for more information on this data structure.
#' @import foreach
#' @import magrittr
#' @examples
#' list_of_fct <- list(nb_fct_simulation(1000, 0.1, 0.1),
#'                     nb_fct_simulation(1000, 0.1, 0.1)) # two replicates
#' direct_optimise_replicates(list_of_fct)
#'
#' @export
direct_optimise_replicates <- function(fct_list,
                                       penalty = NULL, lambda = NULL,
                                       search_scheme = "grid",
                                       multiplier = 10, c_seq_len = 100, # should add a 'by' option
                                       starts = NULL,
                                       alpha_min = 1e-10, delta_min = 1e-10,
                                       alpha_max = 1e5, delta_max = 1e5) {
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

  outcome <- matrix(NA, nrow = c_seq_len+1, ncol = 4)
  outcome_allstarts <- matrix(NA, nrow = c_seq_len*(nrow(starts)+2), ncol = 6)
  a <- 0 # a is a counter for the row of outcome_allstarts
  b <- 0 # b is a counter for the row of outcome

  cc_max <- cc_list %>% unlist %>% max

  switch(search_scheme,
         "grid" = {
           environment(grid_search_replicates) <- environment()
           grid_search_replicates()
         },
         "bisection" = {
           environment(bisection_search_replicates) <- environment()
           bisection_search_replicates()
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

#' @keywords internal
grid_search_replicates <- function() {
  outcome <- matrix(NA, nrow = c_seq_len+1, ncol = 4)
  outcome_allstarts <- matrix(NA, nrow = c_seq_len*(nrow(starts)+2), ncol = 6)
  c_vector = ceiling(seq(cc_max, multiplier*cc_max,length.out = c_seq_len))

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

    environment(optimise_fixed_ccc_replicates) <- environment()
    optimise_fixed_ccc_replicates() # brings this_c_outcome back to this environment.

    outcome_allstarts[(a+1):(a+nrow(starts)),] <- this_c_outcome
    outcome[b,] <- this_c_outcome[which(this_c_outcome[,2] == max(this_c_outcome[,2]))[1],1:4]
    a <- a+nrow(starts)
  }

  outcome <<- outcome
  outcome_allstarts <<- outcome_allstarts
}


#' @keywords internal
bisection_search_replicates <- function() {
  c_init = ceiling(seq(cc_max,multiplier*cc_max,length=4))

  for (i in 1:length(c_init)) {
    b <- b+1
    starts <- input_starts # redundant
    this_c_outcome <- matrix(NA, nrow = nrow(starts),ncol = 6)
    estimate <- c_init[i]

    environment(optimise_fixed_ccc_replicates) <- environment()
    optimise_fixed_ccc_replicates()

    outcome_allstarts[(a+1):(a+nrow(starts)),] <- this_c_outcome
    outcome[b,] <- this_c_outcome[which(this_c_outcome[,2] == max(this_c_outcome[,2]))[1],1:4]
    a <- a+nrow(starts)
  }

  stop_flag <- F
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

        environment(optimise_fixed_ccc_replicates) <- environment()
        optimise_fixed_ccc_replicates()

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

optimise_fixed_ccc_replicates <- function() {
  this_c_outcome <- foreach(j=1:nrow(starts), .combine='rbind') %do% {
    x <- starts[j,]
    y <- optim(x, replicate_likelihood,
               lower = list(alpha = alpha_min, delta = delta_min),
               upper = list(alpha = alpha_max, delta = delta_max), method = "L-BFGS-B",
               control = list(fnscale = -1, maxit = 100),
               ccc = estimate,
               cc_list = cc_list,
               ks_list = ks_list,
               fs_list = fs_list, penalty = penalty, lambda = lambda)
    c(estimate, y$value, y$par[1], y$par[2], x[1], x[2])
  }
  if (is.vector(this_c_outcome)) {
    this_c_outcome <- matrix(this_c_outcome, nrow=1)
  }
  this_c_outcome <<- this_c_outcome
}

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
