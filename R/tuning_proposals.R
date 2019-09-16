#' @title Method 0: Unregularized MLE
#'
#' @description Maximum likelihood estimate without regularization.
#'
#' @details This is used as the comparison point for our tuning proposals
#'   and amounts to a wrapper for \code{\link{direct_optimise_replicates}}. The
#'   output \code{selected_lambda} is always NA for this method, but we format
#'   in this way for consistency with methods 1-4.
#'
#' @param fct_list A list of frequency count tables, assumed to be biological
#'   replicates.
#' @param starts Starting values for \code{alpha} and \code{delta} in the MLE
#'   procedure.
#' @param multiplier The upper bound of the grid of candidate C values, stated in terms of a multiple of the maximum observed richess (c).  For example if c is 50 and multiplier is 10, the method evaluates the likelihood in a C grid from 50 to 500.
#' @param c_seq_len The number of points in the C grid search.
#' @examples
#' \donttest{
#' unregularized_mle(nb_fct_simulation(100, 0.1, 0.1, 2))
#' }
#' @export
unregularized_mle <- function(fct_list,
                              starts = data.frame(alpha = c(1e-2,1e-2),
                                                  delta = c(1e-2,1e-4)),
                              multiplier = 20,
                              c_seq_len = 96,
                              ...) {
  start_time <- Sys.time()
  cc_max <- get_cc_max(fct_list)
  result <- direct_optimise_replicates(fct_list,
                                       penalty = "h1",
                                       lambda = 0,
                                       search_scheme = "grid",
                                       multiplier = multiplier,
                                       c_seq_len = c_seq_len,
                                       starts = starts)
  res <- result[["best"]]
  best <- data.frame(selected_lambda = NA,
                     ccc_hat = res["ccc"],
                     likelihood = res["likelihood"],
                     alpha_hat = res["alpha"],
                     delta_hat = res["delta"],
                     cc_max = cc_max)

  print(str(best))

  result[["best"]] <- best
  result[["time"]] <- Sys.time()-start_time

  # warn the user if the estimate is within 10% of the end of the grid.
  if (((best[1,"ccc"] - cc_max) / (cc_max*multiplier - cc_max)) > 0.9) {
    warning('The estimate is near the end of the grid, consider increasing multiplier')
  }

  return(result)
}

#' @title Method 1: Minimum subset distance
#'
#' @description The regularization parameter \eqn{\lambda} is chosen for its
#'   ability to produce subset estimates with low between-subset variance.
#'
#' @details Method 1 is motivated by the belief that if we resample from the
#'   same population, an ideal \eqn{C} estimator should have low variance.
#'   Exploiting the fact that we have replicate data, the idea is to repeatedly
#'   partition the replicates into two subsets and come up with two estimates.
#'   We select the \eqn{\lambda} which yields the lowest between-subset
#'   variance.  This partitioning is repeated \code{partitions} times to average
#'   out the arbitrary choice of subsets.  See paper or source code for more
#'   detail.
#'
#' @param fct_list A list of frequency count tables, assumed to be replicates.
#' @param lambda_vec The values of the penalty parameter we select from.
#' @param starts Starting values for \code{alpha} and \code{delta} in the MLE
#'   procedure.
#' @param partitions An integer indicating the number of times to partition the
#'   data into two subsets
#' @param multiplier The upper bound of the grid of candidate C values, stated in terms of a multiple of the maximum observed richess (c).  For example if c is 50 and multiplier is 10, the method evaluates the likelihood in a C grid from 50 to 500.
#' @param c_seq_len The number of points in the C grid search.
#' @examples
#' \donttest{
#' minimum_subset_distance(nb_fct_simulation(100, 0.1, 0.1, 2))
#' }
#' @export
minimum_subset_distance <- function(fct_list,
                                    lambda_vec = seq(0, 20, by= 2),
                                    starts = data.frame(alpha = c(1e-2,1e-2),
                                                        delta = c(1e-2,1e-4)),
                                    partitions = 10,
                                    multiplier = 20,
                                    c_seq_len = 96,
                                    ...) {
  start_time <- Sys.time()
  r <- length(fct_list)
  if (r < 2) {
    stop("We cannot use this method with less than 2 frequency count tables.")
  }

  cc_max <- get_cc_max(fct_list)

  result_list <- list() # each result in a list until the end
  for (k in 1:partitions) {
    subset_ind <- sample(1:r, round(r/2))
    subset_one <- fct_list[subset_ind]
    subset_two <- fct_list[-subset_ind]

    # Just a string so we can look at the subsetting later if desired:
    subset_str <- paste0( "{", paste0(subset_ind, collapse = ","), "} {",
                          paste0(base::setdiff(1:r, subset_ind),
                                 collapse = ","), "}")

    print(paste0("Working on the partition ",subset_str," at ", Sys.time())) # for cluster so I can monitor.

    # This step is required due to the way the likelihood function is written:
    formed_lists_one <- make_formatted_lists(subset_one)
    cc_list_one <- formed_lists_one$cc_list
    fs_list_one <- formed_lists_one$fs_list
    ks_list_one <- formed_lists_one$ks_list


    formed_lists_two <- make_formatted_lists(subset_two)
    cc_list_two <- formed_lists_two$cc_list
    fs_list_two <- formed_lists_two$fs_list
    ks_list_two <- formed_lists_two$ks_list

    for (lam in lambda_vec) {
      optimise_with_settings <- function(x) {
        res <- direct_optimise_replicates(x, penalty = "h1",
                                          lambda = lam,
                                          search_scheme = "grid",
                                          multiplier = multiplier,
                                          c_seq_len = c_seq_len,
                                          starts = starts)
        res[["best"]]
      }
      result_one <- optimise_with_settings(subset_one)
      result_two <- optimise_with_settings(subset_two)

      ccc_hat_one <- result_one[["ccc"]]
      ccc_hat_two <- result_two[["ccc"]]
      alpha_hat <- mean(c(result_one[["alpha"]], result_two[["alpha"]]))
      delta_hat <- mean(c(result_one[["delta"]], result_two[["delta"]]))


      distance <- var(c(ccc_hat_one, ccc_hat_two))
      ccc_hat <- mean(c(ccc_hat_one, ccc_hat_two))

      df <- data.table::data.table(partition = k,
                                   part_str = subset_str,
                                   distance = distance,
                                   lambda = lam,
                                   ccc_hat = ccc_hat,
                                   alpha_hat = alpha_hat,
                                   delta_hat = delta_hat)
      result_list <- c(result_list, list(df))
    }
  }
  full_results <- data.table::rbindlist(result_list)

  selected_lambda <- full_results %>%
    dplyr::group_by(lambda) %>%
    dplyr::summarize(mean_distance = mean(distance)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(mean_distance) %>%
    dplyr::filter(1:dplyr::n() == 1) %>%
    dplyr::select("lambda") %>%
    as.numeric

  best <- full_results %>%
    dplyr::filter(lambda == selected_lambda) %>%
    # all of the following mean commands are average over the partitions
    dplyr::summarize(selected_lambda = mean(lambda),
                     ccc_hat = ceiling(mean(ccc_hat)),
                     alpha_hat = mean(alpha_hat),
                     delta_hat = mean(delta_hat),
                     mean_distance = mean(distance))

  # The following case is possible due to separate partition estimates:
  if (best["ccc_hat"] < cc_max) best["ccc_hat"] <- cc_max

  # warn the user if the estimate is near the end of the grid.
  if (((best[1,"ccc_hat"] - cc_max) / (cc_max*multiplier - cc_max)) > 0.9) {
    warning('The estimate is near the end of the grid, consider increasing multiplier')
  }
  # warn the user if selected lambda is near the end of the grid.
  if (length(lambda_vec) > 5) { # if the length is less than 5 we cant do much.
    position <- which(sort(lambda_vec) == best[1,"selected_lambda"])
    prop_of_lambda_vec <- position / length(lambda_vec)
    if (prop_of_lambda_vec > 0.8) {
      warning("Selected lambda is near the edge of grid, consider expanding lambda_vec")
    }
  }

  return(list(best = best, full = full_results, time = (Sys.time()-start_time)))
}


#' @title Method 3: Goodness of fit criterion
#'
#' @description A regularization parameter \eqn{\lambda} is selected using a
#'   goodness of fit metric.
#'
#' @details We generate a \eqn{C} estimate for each \eqn{\lambda} in
#'   \code{lambda_vec}.  Using these estimates we use a \eqn{\chi}-square
#'   goodness of fit statistic to evaluate the fit to the sample.  The
#'   \eqn{\lambda} value with the best fit is \code{selected_lambda}, and the
#'   \eqn{C} estimate associated with that \eqn{\lambda} is \code{ccc_hat}.  See
#'   paper for full details.
#'
#'
#' @param fct_list A list of frequency count tables, assumed to be biological
#'   replicates.
#' @param lambda_vec The values of the penalty parameter we consider in
#'   selecting \eqn{\lambda}.
#' @param starts Starting values for \code{alpha} and \code{delta} in the MLE
#'   procedure.
#' @param gof_method The only option currently supported is "chi_sq".
#' @param multiplier The upper bound of the grid of candidate C values, stated in terms of a multiple of the maximum observed richess (c).  For example if c is 50 and multiplier is 10, the method evaluates the likelihood in a C grid from 50 to 500.
#' @param c_seq_len The number of points in the C grid search.
#' @examples
#' \donttest{
#' gof_criterion(nb_fct_simulation(100, 0.1, 0.1, 2))
#' }
#' @export
gof_criterion <- function(fct_list,
                          lambda_vec = seq(0, 20, by= 2),
                          starts = data.frame(alpha = c(1e-2,1e-2),
                                              delta = c(1e-2,1e-4)),
                          gof_method = "chi_sq",
                          multiplier = 20,
                          c_seq_len = 96,
                          ...) {
  start_time <- Sys.time()
  cc_max <- get_cc_max(fct_list)
  supported_methods <- c("chi_sq")
  if ( !(gof_method %in% supported_methods) ) {
    stop(paste0("The supported methods are ",
                paste(supported_methods, collapse = ", ")))
  }
  r <- length(fct_list)
  result_list <- list() # each result in a list until the end
  for (lam in lambda_vec) {
    print(paste0("Running with lambda = ",lam," at ", Sys.time()))
    optimise_with_settings <- function(x) {
      res <- direct_optimise_replicates(x, penalty = "h1",
                                        lambda = lam,
                                        search_scheme = "grid",
                                        multiplier = multiplier,
                                        c_seq_len = c_seq_len,
                                        starts = starts)
      res[["best"]]
    }
    result <- optimise_with_settings(fct_list)

    ccc_hat <- result[,"ccc"] %>% unlist
    alpha_hat <- result[,"alpha"] %>% unlist
    delta_hat <- result[,"delta"] %>% unlist
    if (gof_method == "chi_sq") {
      good_of_fit <- lapply(fct_list, chi_sq_gof,
                            ccc_hat = ccc_hat,
                            alpha_hat = alpha_hat,
                            delta_hat = delta_hat) %>%
        unlist %>% sum
    } else {
      stop("Unknown gof_method passed.")
    }

    df <- data.table::data.table(lambda = lam,
                                 good_of_fit = good_of_fit,
                                 ccc_hat = ccc_hat,
                                 cc_max = cc_max,
                                 alpha_hat = alpha_hat,
                                 delta_hat = delta_hat)
    result_list <- c(result_list, list(df))
  }
  full_results <- data.table::rbindlist(result_list)

  selected_lambda <- full_results %>%
    dplyr::arrange(good_of_fit) %>%
    dplyr::filter(1:dplyr::n() == 1) %>%
    dplyr::select("lambda") %>%
    as.numeric

  best <- full_results %>%
    dplyr::filter(lambda == selected_lambda) %>%
    dplyr::rename(selected_lambda = lambda)

  # warn the user if the estimate is within 10% of the end of the grid.
  if (((best[1,"ccc_hat"] - cc_max) / (cc_max*multiplier - cc_max)) > 0.9) {
    warning('The estimate is near the end of the grid, consider increasing multiplier')
  }

  # warn the user if selected lambda is near the end of the grid.
  if (length(lambda_vec) > 5) { # if the length is less than 5 we cant do much.
    position <- which(sort(lambda_vec) == best[1,"selected_lambda"])
    prop_of_lambda_vec <- position / length(lambda_vec)
    if (prop_of_lambda_vec > 0.8) {
      warning("Selected lambda is near the edge of grid, consider expanding lambda_vec")
    }
  }


  return(list(best = best, full = full_results, time = (Sys.time()-start_time)))
}

#' @title Method 2 and Method 4
#'
#' @description Method 2 is cross-validation using the likelihood in the
#'   evaluation step.  Method 4 is cross-validation using the goodness of fit
#'   statistic in the evaluation step.
#'
#' @details Methods 2 and 4 have very similar structure we we've included them
#'   both in the same function.  To run each method use: \enumerate{ \item
#'   Method 2:  cv_replicates(..., "neg_unreg_like") \item Method 4:
#'   cv_replicates(..., "gof_chi_sq") } In each method we partition the data
#'   \code{partitions} times into training and evaluation subsets.  An estimate
#'   for each \eqn{\lambda} in \code{lambda_vec} is generated and we evaluate
#'   them using the evaluation subset.  The evaluation step depends on the
#'   method, see paper or source code for details of how these functions work.
#'
#'
#' @param fct_list A list of frequency count tables, assumed to be replicates.
#' @param lambda_vec The values of the penalty parameter we consider in
#'   selecting \eqn{\lambda}.
#' @param starts Starting values for \code{alpha} and \code{delta} in the MLE
#'   procedure.
#' @param partitions An integer indicating the number of times to randomly split
#'   the data into testing and validating subsets.
#' @param eval_function A function which evaluates how well a set of parameters
#'   fit a list of frequency count tables.  To conform to goodness of fit, we
#'   use the negative of the likelihood function so that low scores are better.
#' @param multiplier The upper bound of the grid of candidate C values, stated in terms of a multiple of the maximum observed richess (c).  For example if c is 50 and multiplier is 10, the method evaluates the likelihood in a C grid from 50 to 500.
#' @param c_seq_len The number of points in the C grid search.
#' @examples
#' \donttest{
#' cv_replicates(nb_fct_simulation(100, 0.1, 0.1, 2))
#' }
#' @export
cv_replicates <- function(fct_list,
                          lambda_vec = seq(0, 20, by= 2),
                          starts = data.frame(alpha = c(1e-2,1e-2),
                                              delta = c(1e-2,1e-4)),
                          partitions = 10,
                          eval_function = "gof_chi_sq",
                          multiplier = 20,
                          c_seq_length = 96,
                          ...) {
  start_time <- Sys.time()
  r <- length(fct_list)
  if (r < 2) {
    stop("We cannot use this method with less than 2 frequency count tables.")
  }
  valid_functions <- c("gof_chi_sq", "neg_unreg_like")
  if (!(eval_function %in% valid_functions)) {
    stop(paste0("eval_function must be one of: ",
                paste(valid_functions, collapse = ", ")))
  }

  cc_max <- get_cc_max(fct_list)

  result_list <- list() # each result in a list until the end
  for (k in 1:partitions) {
    subset_ind <- sample(1:r, round(r/2))
    subset_train <- fct_list[subset_ind]
    subset_eval <- fct_list[-subset_ind]

    # Just a string so we can look at the subsetting later if desired:
    subset_str <- paste0( "{", paste0(subset_ind, collapse = ","), "} {",
                          paste0(base::setdiff(1:r, subset_ind),
                                 collapse = ","), "}")

    print(paste0("Working on the partition ",subset_str," at ", Sys.time()))

    # This step is required due to the way the likelihood function is written:
    if (eval_function == "neg_unreg_like") {
      formed_lists <- make_formatted_lists(subset_eval)
      cc_list <- formed_lists$cc_list
      fs_list <- formed_lists$fs_list
      ks_list <- formed_lists$ks_list
    }

    for (lam in lambda_vec) {
      optimise_with_settings <- function(x) {
        res <- direct_optimise_replicates(x, penalty = "h1",
                                          lambda = lam,
                                          search_scheme = "grid",
                                          multiplier = multiplier,
                                          c_seq_len = c_seq_length,
                                          starts = starts)
        res[["best"]]
      }
      training_result <- optimise_with_settings(subset_train)

      ccc_hat <- training_result[["ccc"]]
      if (ccc_hat < cc_max) {
        ccc_hat <- cc_max
        # This is a very sensible ad-hoc adjustment, however
        # it is not absolutely necessary for the gof eval fns.
      }
      alpha_hat <- training_result[["alpha"]]
      delta_hat <- training_result[["delta"]]

      if (eval_function == "gof_chi_sq") {
        fn_val <- lapply(subset_eval, chi_sq_gof,
                         ccc_hat = ccc_hat,
                         alpha_hat = alpha_hat,
                         delta_hat = delta_hat) %>%
          unlist %>% sum
      } else if (eval_function == "neg_unreg_like") {
        fn_val <- rre::replicate_likelihood(x = c(alpha_hat, delta_hat),
                                            ccc = ccc_hat,
                                            cc_list = cc_list,
                                            ks_list = ks_list,
                                            fs_list = fs_list,
                                            penalty = "h1",
                                            lambda = 0)
        fn_val <- fn_val * -1 # so that lower is better.
      } else {
        stop("This error should never happen, but just in case.")
      }

      df <- data.table::data.table(partition = k,
                                   part_str = subset_str,
                                   eval_fn_val = fn_val,
                                   lambda = lam,
                                   ccc_hat = ccc_hat,
                                   alpha_hat = alpha_hat,
                                   delta_hat = delta_hat)
      result_list <- c(result_list, list(df))
    }
  }
  full_results <- data.table::rbindlist(result_list)

  selected_lambda <- full_results %>%
    dplyr::group_by(lambda) %>%
    dplyr::summarize(mean_fn_val = mean(eval_fn_val)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(mean_fn_val) %>%
    dplyr::filter(1:dplyr::n() == 1) %>%
    dplyr::select("lambda") %>%
    as.numeric

  best <- full_results %>%
    dplyr::filter(lambda == selected_lambda) %>%
    # all of the following mean commands are average over the partitions
    dplyr::summarize(selected_lambda = mean(lambda),
                     ccc_hat = ceiling(mean(ccc_hat)),
                     alpha_hat = mean(alpha_hat),
                     delta_hat = mean(delta_hat),
                     mean_eval_fn_val = mean(eval_fn_val))

  # warn the user if the estimate is within 10% of the end of the grid.
  if (((best[1,"ccc_hat"] - cc_max) / (cc_max*multiplier - cc_max)) > 0.9) {
    warning('The estimate is near the end of the grid, consider increasing multiplier')
  }

  # warn the user if selected lambda is near the end of the grid.
  if (length(lambda_vec) > 5) { # if the length is less than 5 we cant do much.
    position <- which(sort(lambda_vec) == best[1,"selected_lambda"])
    prop_of_lambda_vec <- position / length(lambda_vec)
    if (prop_of_lambda_vec > 0.8) {
      warning("Selected lambda is near the edge of grid, consider expanding lambda_vec")
    }
  }


  return(list(best = best, full = full_results, time = (Sys.time()-start_time)))
}

#' @title Optimizing the likelihood with a fixed lambda value
#'
#' @description A wrapper function for direct_optimise_replicates to standarize
#'   output for simulations.
#'
#' @details A wrapper function for direct_optimise_replicates to standarize
#'   output for simulations.
#'
#'
#' @param fct_list A list of frequency count tables, assumed to be biological
#'   replicates.
#' @param starts Starting values for \code{alpha} and \code{delta} in the MLE
#'   procedure.
#' @param lambda The fixed lambda value to run the MLE procedure with.
#' @param multiplier The upper bound of the grid of candidate C values, stated in terms of a multiple of the maximum observed richess (c).  For example if c is 50 and multiplier is 10, the method evaluates the likelihood in a C grid from 50 to 500.
#' @param c_seq_len The number of points in the C grid.
fixed_lambda_mle <- function(fct_list,
                             starts = data.frame(alpha = c(1e-2,1e-2),
                                                 delta = c(1e-2,1e-4)),
                             lambda = 0,
                             multiplier = 20,
                             c_seq_length = 96,
                             ...) {
  start_time <- Sys.time()
  result <- direct_optimise_replicates(fct_list,
                                       penalty = "h1",
                                       lambda = lambda,
                                       search_scheme = "grid",
                                       multiplier = multiplier,
                                       c_seq_len = c_seq_len,
                                       starts = starts)
  res <- result[["best"]]
  best <- data.frame(selected_lambda = lambda,
                     ccc_hat = res["ccc"],
                     likelihood = res["likelihood"],
                     alpha_hat = res["alpha"],
                     delta_hat = res["delta"],
                     cc_max = get_cc_max(fct_list))

  result[["best"]] <- best
  result[["time"]] <- Sys.time()-start_time
  return(result)
}

#' @title Supplemental method: single FCT subset estimates
#'
#' @description Divide the list of FCT into single FCT subsets, then compare
#'   those estimates using variance, cv, or gini coefficient.
#'
#' @details This method was explored after the paper was written and simulations
#'   completed as a suggested vairant to \code{\link{minimum_subset_distance()}}
#'
#' @param fct_list A list of frequency count tables, assumed to be replicates.
#' @param lambda_vec The values of the penalty parameter we will train over.
#' @param starts Starting values for \code{alpha} and \code{delta} in the MLE
#'   procedure.
#' @param metric A string which is "variance", "cv", "index_of_dispersion" or
#'   "gini"
#' @examples
#' \donttest{
#' single_fct_subsets(nb_fct_simulation(100, 0.1, 0.1))
#' }
single_fct_subsets <- function(fct_list,
                               lambda_vec = seq(0, 20, by= 2),
                               starts = data.frame(alpha = c(1e-2,1e-2),
                                                   delta = c(1e-2,1e-4)),
                               metric = "variance",
                               ...) {
  start_time <- Sys.time()
  r <- length(fct_list)
  if (r < 2) {
    stop("We cannot use this method with less than 2 frequency count tables.")
  }

  if (!(metric %in% c("variance", "cv", "index_of_dispersion", "gini"))) {
    stop("Unknown metric.")
  }

  metric_fn <- switch(metric, #replace by your input
                      "variance" = (function(x) var(x)),
                      "cv" = (function(x) (sqrt(var(x)))/mean(x)),
                      "index_of_dispersion" = (function(x) (var(x)/mean(x))),
                      "gini" = (function(x) {
                        n <- length(x)
                        gi <- 0
                        for(i in 1:n) {
                          for(j in 1:n) {
                            gi <- gi + abs(x[i] - x[j])
                          }
                        }
                        gi <- gi/(2*n*sum(x))
                        return(gi)
                      })
  )
  print(paste0("Running method with ",r," FCT using ", metric,
               ", which applies the following function to the C estimates:"))
  print(metric_fn)

  cc_max <- get_cc_max(fct_list)

  result_list <- list() # each result in a list until the end

  for (lam in lambda_vec) {
    print(paste0("Working on lambda = ",lam," at ", Sys.time()))
    optimise_with_settings <- function(x) {
      res <- direct_optimise(x, penalty = "h1",
                             lambda = lam,
                             search_scheme = "grid",
                             multiplier = 20, c_seq_len = 96,
                             starts = starts,
                             forced_ccc_lower_bound = cc_max)
      res[["best"]]
    }

    for (j in 1:r) {
      fct <- fct_list[[j]]
      this_fct_result <- optimise_with_settings(fct)
      this_result_as_df <- data.frame(fct_id = j,
                                      lambda = lam,
                                      ccc_hat = this_fct_result['ccc'],
                                      alpha_hat = this_fct_result['alpha'],
                                      delta_hat = this_fct_result['delta'])
      result_list <- c(result_list, list(this_result_as_df))
    }
  }
  full_results <- data.table::rbindlist(result_list)
  full_results %<>%
    dplyr::group_by(lambda) %>%
    dplyr::summarize(mean_ccc_hat = mean(ccc_hat), # temporary name so we can calulate
                     alpha_hat = mean(alpha_hat),
                     delta_hat = mean(delta_hat),
                     metric = metric_fn(ccc_hat)
    ) %>%
    dplyr::rename(ccc_hat = mean_ccc_hat) # change name back

  selected_lambda <- full_results %>%
    dplyr::arrange(metric) %>%
    dplyr::filter(1:dplyr::n() == 1) %>%
    dplyr::select("lambda") %>%
    as.numeric

  best <- full_results %>%
    dplyr::filter(lambda == selected_lambda) %>%
    # all of the following mean commands are average over the partitions
    dplyr::summarize(selected_lambda = mean(lambda), # mean of a constant
                     ccc_hat = mean(ccc_hat),
                     alpha_hat = mean(alpha_hat),
                     delta_hat = mean(delta_hat),
                     metric = mean(metric))
  return(list(best = best, full = full_results, time = (Sys.time()-start_time)))
}
