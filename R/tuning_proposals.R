#' @title Proposal 0: Unregularized MLE wrapper for simulations
#'
#' @description Maximum likelihood estimate without regularization.
#'
#' @details This is used as the comparison point for our tuning proposals, a
#'   wrapper to conform to simulation output.
#'
#'
#' @param fct_list A list of frequency count tables, assumed to be biological
#'   replicates.
#' @param starts Starting values for \code{alpha} and \code{delta} in the MLE
#'   procedure.
#' @export
unregularized_mle <- function(fct_list,
                              starts = data.frame(alpha = c(1e-2,1e-2),
                                                  delta = c(1e-2,1e-4)),
                              ...) {
  result <- direct_optimise_replicates(fct_list,
                                       penalty = "h1",
                                       lambda = 0,
                                       search_scheme = "grid",
                                       multiplier = 20, c_seq_len = 96,
                                       starts = starts)
  res <- result[["best"]]
  best <- data.frame(selected_lambda = NA,
                     ccc_hat = res["ccc"],
                     likelihood = res["likelihood"],
                     alpha_hat = res["alpha"],
                     delta_hat = res["delta"],
                     cc_max = get_cc_max(fct_list))

  result[["best"]] <- best

  return(result)
}

#' @title Proposal 1: Minimum subset distance
#'
#' @description Minimum variance between C estimates in each subset.
#'
#' @details
#'
#' @param fct_list A list of frequency count tables, assumed to be replicates.
#' @param lambda_vec The values of the penalty parameter we will train over.
#' @param starts Starting values for \code{alpha} and \code{delta} in the MLE procedure.
#' @param partitions An integer indicating the number of times to randomly split the data into testing and validating subsets.
#' @export
minimum_subset_distance <- function(fct_list,
                                    lambda_vec = seq(0, 20, by= 2),
                                    starts = data.frame(alpha = c(1e-2,1e-2),
                                                        delta = c(1e-2,1e-4)),
                                    partitions = 10,
                                    ...) {
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
                                          multiplier = 20, c_seq_len = 96,
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
    dplyr::filter(1:n() == 1) %>%
    dplyr::select("lambda") %>%
    as.numeric

  best <- full_results %>%
    dplyr::filter(lambda == selected_lambda) %>%
    # all of the following mean commands are average over the partitions
    dplyr::summarize(selected_lambda = mean(lambda),
                     ccc_hat = mean(ccc_hat),
                     alpha_hat = mean(alpha_hat),
                     delta_hat = mean(delta_hat),
                     mean_distance = mean(distance))
  return(list(best = best, full = full_results))
}


#' @keywords internal
# I'm not planning to export this, its not really a candidate method at all.
# This is an intentionally naive function that optimizes each replicate separately and then averages the c_hat, alpha_hat, delta_hat, etc.
naive_replicates <- function(fct_list) {
  starts <- expand.grid(alpha = c(1e-1,1e-2),
                        delta = c(1e-1,1e-2,1e-3,1e-4))
  results <- lapply(fct_list, direct_optimise_replicates,
                    starts = starts,
                    multiplier = 30) %>%
    lapply(.,'[[',"best") %>%
    do.call(rbind,.) %>%
    colMeans
  best <- c(c_hat = results["ccc"],
            alpha_hat = results["alpha"],
            delta_hat = results["delta"],
            lambda_chosen = NA,
            splits = NA,
            c_hat_lam_zero = results["ccc"])
  best
}

#' @title Goodness of fit regularization criterion
#'
#' @description Optimizes with a grid of regularization parameters and chooses the best one via goodness of fit metrics.
#'
#' @details Currently only chi^2 is supported, but there are potentially better goodness of fit metrics for count data.
#'
#'
#' @param fct_list A list of frequency count tables, assumed to be biological replicates.
#' @param lambda_vec The values of the penalty parameter we will train over.
#' @param starts Starting values for \code{alpha} and \code{delta} in the MLE procedure.
#' @param gof_method Goodness of fit function, currently chi^2 is the only option.
#' @export
gof_criterion <- function(fct_list,
                          lambda_vec = seq(0, 20, by= 2),
                          starts = data.frame(alpha = c(1e-2,1e-2),
                                              delta = c(1e-2,1e-4)),
                          gof_method = "chi_sq",
                          ...) {
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
                                        multiplier = 20, c_seq_len = 96,
                                        starts = starts)
      res[["best"]]
    }
    result <- optimise_with_settings(fct_list)

    ccc_hat <- result["ccc"]
    alpha_hat <- result["alpha"]
    delta_hat <- result["delta"]
    if (gof_method == "chi_sq") {
      good_of_fit <- lapply(fct_list, chi_sq_gof,
                            ccc_hat = ccc_hat,
                            alpha_hat = alpha_hat,
                            delta_hat = delta_hat) %>%
        unlist %>% sum
    } else {
      stop("The existence of this error indicates a coding bug.")
    }

    cc_max <- get_cc_max(fct_list)

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
    dplyr::filter(1:n() == 1) %>%
    dplyr::select("lambda") %>%
    as.numeric

  best <- full_results %>%
    dplyr::filter(lambda == selected_lambda) %>%
    dplyr::rename(selected_lambda = lambda)
  return(list(best = best, full = full_results))
}

#' @title Proposal 4: Cross validation for data with replicates
#'
#' @description Randomly partitions the data into training and evaluation subsets and selects the best parameter estimates based on various evaluation metrics.
#'
#' @details For each partition the we find the MLE solution, repeating for a grid of regularization values given by \code{lambda_vec}.  We then evaluate how well these parameters fit using the training data and a choice of evaluation functions.  Currently we support goodness of fit (chi-square) and unregularized likelihood as possible evaluation functions.
#'
#'
#' @param fct_list A list of frequency count tables, assumed to be replicates.
#' @param lambda_vec The values of the penalty parameter we will train over.
#' @param starts Starting values for \code{alpha} and \code{delta} in the MLE procedure.
#' @param partitions An integer indicating the number of times to randomly split the data into testing and validating subsets.
#' @param eval_function A function which evaluates how well a set of parameters fit a list of frequency count tables.  To conform to goodness of fit, we will assume that lower scores are better.
#' @export
cv_replicates <- function(fct_list,
                          lambda_vec = seq(0, 20, by= 2),
                          starts = data.frame(alpha = c(1e-2,1e-2),
                                              delta = c(1e-2,1e-4)),
                          partitions = 10,
                          eval_function = "gof_chi_sq",
                          ...) {
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
                                          multiplier = 20, c_seq_len = 96,
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
    dplyr::filter(1:n() == 1) %>%
    dplyr::select("lambda") %>%
    as.numeric

  best <- full_results %>%
    dplyr::filter(lambda == selected_lambda) %>%
    # all of the following mean commands are average over the partitions
    dplyr::summarize(selected_lambda = mean(lambda),
                     ccc_hat = mean(ccc_hat),
                     alpha_hat = mean(alpha_hat),
                     delta_hat = mean(delta_hat),
                     mean_eval_fn_val = mean(eval_fn_val))
  return(list(best = best, full = full_results))
}

#' @title Fixed lambda MLE for supplemental simulations
#'
#' @description
#'
#' @details
#'
#'
#' @param fct_list A list of frequency count tables, assumed to be biological replicates.
#' @param starts Starting values for \code{alpha} and \code{delta} in the MLE procedure.
#' @param lambda The fixed lambda value to run the MLE procedure with.
#' @export
fixed_lambda_mle <- function(fct_list,
                             starts = data.frame(alpha = c(1e-2,1e-2),
                                                 delta = c(1e-2,1e-4)),
                             lambda = 0,
                             ...) {
  result <- direct_optimise_replicates(fct_list,
                                       penalty = "h1",
                                       lambda = lambda,
                                       search_scheme = "grid",
                                       multiplier = 20, c_seq_len = 96,
                                       starts = starts)
  res <- result[["best"]]
  best <- data.frame(selected_lambda = lambda,
                     ccc_hat = res["ccc"],
                     likelihood = res["likelihood"],
                     alpha_hat = res["alpha"],
                     delta_hat = res["delta"],
                     cc_max = get_cc_max(fct_list))

  result[["best"]] <- best

  return(result)
}

#' @title Proposal 5: single FCT subset estimates
#'
#' @description Dividing the list of FCT into single FCT and then comparing
#'   those estimates using variance, cv, or gini (metric)
#'
#' @details
#'
#' @param fct_list A list of frequency count tables, assumed to be replicates.
#' @param lambda_vec The values of the penalty parameter we will train over.
#' @param starts Starting values for \code{alpha} and \code{delta} in the MLE procedure.
#' @param
#' @export
single_fct_subsets <- function(fct_list,
                               lambda_vec = seq(0, 20, by= 2),
                               starts = data.frame(alpha = c(1e-2,1e-2),
                                                   delta = c(1e-2,1e-4)),
                               metric = "variance",
                               ...) {
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
    dplyr::filter(1:n() == 1) %>%
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
  return(list(best = best, full = full_results))
}
