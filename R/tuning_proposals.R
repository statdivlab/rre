#' @title Cross validation using likelihood function
#'
#' @description \code{cross_validate_replicates} is a function which finds the maximum likelihood solution over a set of penalty parameter values, then determines which is best by evaluating the likelihood on the test data.
#'
#' @details Coming soon.  Note to self:  use_pen_likelihood, out_dir and print_details_to_file should be removed before this is released to anyone else.
#'
#' @param fct_list A list of frequency count tables, assumed to be replicates.
#' @param split_method The method for choosing test and training data in each split.
#' @param lambda_vec The values of the penalty parameter we will train over.
#' @param starts Starting values for \code{alpha} and \code{delta} in the MLE procedure.
#' @param use_pen_likelihood Logical for whether we will use the penalized likelihood when evaluating on the test data in each split.
#' @param verbose logical.  If true then we additionally return a data frame with the details from each split.
#' @param print_details_to_file Prints the verbose output to a csv file.  For use during simulations if the split details are of potential interest.
#' @param k The number of splits if \code{split_method = "k_splits"}.
#' @param out_dir The file output for print_details_to_file
#' @export
cross_validate_replicates <- function(fct_list, split_method = "leave_one_out",
                                      lambda_vec = c(0,2^(-5:5),4096),
                                      starts = expand.grid(alpha = c(1e2,1e-2),
                                                           delta = c(1e2,1e-2,1e-3,1e-5)),
                                      use_pen_likelihood = F,
                                      verbose = F,
                                      print_details_to_file = T,
                                      k = 0,
                                      out_dir = paste0(getwd(),"/sim_details"),
                                      ...) {
  if (!(split_method %in% c("leave_one_out","k_splits"))) {
    stop("split_method not recognized")
  }

  r <- length(fct_list)
  if ((split_method == "k_splits" & !(k %in% 2:r))) {
    stop("For the k splits parameter you must pass an integer k which is an integer less than the number of replicates")
  }

  if (print_details_to_file == T) {
    dir.create(out_dir, showWarnings = F)
  }


  splits <- ifelse((split_method=="leave_one_out"),length(fct_list),k)
  cv_dat <- data.frame(split = numeric(splits*length(lambda_vec)),
                       lambda = numeric(splits*length(lambda_vec)),
                       l_unpen = numeric(splits*length(lambda_vec)),
                       l_pen = numeric(splits*length(lambda_vec)),
                       c_hat = numeric(splits*length(lambda_vec)),
                       alpha_hat = numeric(splits*length(lambda_vec)),
                       delta_hat = numeric(splits*length(lambda_vec)))
  df_ind <- 1
  switch(split_method,
         "leave_one_out" = {
           for (split in 1:splits) {
             test_list <- fct_list[split]
             train_list <- fct_list[-split]

             environment(one_split) <- environment()
             one_split()

             cv_dat[df_ind:(df_ind+nrow(result)-1),] <- result
             df_ind = df_ind+nrow(result)
           }
         },
         "k_splits" = {
           cuts = c(round(r/k*(1:(k-1)))+1,r+1)
            for (split in 1:splits) {
             test_indices <- cuts[split]:(cuts[split+1]-1)
             test_list <- fct_list[test_indices]
             train_list <- fct_list[-test_indices]

             environment(one_split) <- environment()
             one_split()

             cv_dat[df_ind:(df_ind+nrow(result)-1),] <- result
             df_ind = df_ind+nrow(result)
           }
         })

  c_hat_lam_zero <- cv_dat[which(cv_dat[,"lambda"] == 0),"c_hat"] %>% mean
  if (use_pen_likelihood) {
    temp <- cv_dat %>% aggregate(formula=l_pen~lambda,FUN=mean, data= .)
    lambda_chosen <- temp[which(temp[,2] == max(temp[,2]))[1],1]
    c_hat <- cv_dat[which(cv_dat[,"lambda"] == lambda_chosen),"c_hat"] %>% mean
    alpha_hat <- cv_dat[which(cv_dat[,"lambda"] == lambda_chosen),"alpha_hat"] %>% mean
    delta_hat <- cv_dat[which(cv_dat[,"lambda"] == lambda_chosen),"delta_hat"] %>% mean
  } else {
    temp <- cv_dat %>% aggregate(formula=l_unpen~lambda,FUN=mean, data= .)
    lambda_chosen <- temp[which(temp[,2] == max(temp[,2]))[1],1]
    c_hat <- cv_dat[which(cv_dat[,"lambda"] == lambda_chosen),"c_hat"] %>% mean
    alpha_hat <- cv_dat[which(cv_dat[,"lambda"] == lambda_chosen),"alpha_hat"] %>% mean
    delta_hat <- cv_dat[which(cv_dat[,"lambda"] == lambda_chosen),"delta_hat"] %>% mean
  }

  cv_best <- c(c_hat = c_hat,
               alpha_hat = alpha_hat,
               delta_hat = delta_hat,
               lambda_chosen = lambda_chosen,
               folds = ifelse((split_method=="leave_one_out"),
                                        length(fct_list),k),
               c_hat_lam_zero = c_hat_lam_zero)

  if (print_details_to_file) {
    write_csv(cv_dat, path = paste0(out_dir, "/",
                                    (Sys.time() %>%
                                        str_replace_all(., " ","_") %>%
                                        str_replace_all(., ":","_")),
                                    ".csv"))
  }

  if (verbose) {
    return(list(cv_best = cv_best, cv_dat = cv_dat))
  } else {
    return(cv_best)
  }
}


#' @title Helper function for \code{cross_validate_replicates}
#' @keywords internal
#  Runs the cross validation procedure for one split. pen_likelihood is a
#T/F flag for whether we should also store the penalized likelihood in addition
#to the unpenalized.  This is a design decision which we hope to resolve after
#more thought. Note: It can now be removed, I just haven't done it yet.
one_split <- function(pen_likelihood = T){

  split_df <- data.frame(lambda = lambda_vec,
                         likelihood_unpenalised = 0)
  if (pen_likelihood) {
    split_df$likelihood_penalised = 0
  }
  split_df$c_hat = 0
  split_df$alpha_hat = 0
  split_df$delta_hat = 0

  test_list_formatted <- make_formatted_lists(test_list)
  cc_list <- test_list_formatted$cc_list
  fs_list <- test_list_formatted$fs_list
  ks_list <- test_list_formatted$ks_list

  for (lambda in lambda_vec) {
    MLE <- direct_optimise_replicates(train_list,
                                      starts = starts,
                                      penalty = "h1",
                                      lambda = lambda,
                                      ...)[["best"]]
    alpha <- MLE["alpha"]
    delta <- MLE["delta"]
    ccc <-MLE["ccc"]
    x <- c(alpha,delta)

    split_df[which(split_df$lambda == lambda),"c_hat"] <- ccc
    split_df[which(split_df$lambda == lambda),"alpha_hat"] <- alpha
    split_df[which(split_df$lambda == lambda),"delta_hat"] <- delta



    # Calculate the likelihood on the test set:
    if (ccc < max(unlist(cc_list))) {
      ccc_old <-ccc
      ccc <- max(unlist(cc_list))
      print(paste0("Ad hoc line that shifts C used with lambda = ",lambda))
    }


    l_unpen <- replicate_likelihood(x,ccc,cc_list,
                                    ks_list,fs_list,
                                    penalty = NULL, lambda = NULL)

    l_pen <- replicate_likelihood(x,ccc,cc_list,
                                  ks_list,fs_list,
                                  penalty = "h1", lambda = lambda)


    split_df[which(split_df$lambda == lambda),"likelihood_unpenalised"] <- l_unpen
    if (pen_likelihood) {
      split_df[which(split_df$lambda == lambda),"likelihood_penalised"] <- l_pen
    }
  }
  split_filler <- rep(split,times = nrow(split_df))
  result <- cbind(split_filler,split_df)
  result <<- result
}


#' @title Tuning penalty parameter by minimum variance criterion
#'
#' @description This function chooses the \code{lambda} value which produces C
#'   estimates with the smallest variance.
#'
#' @details We found by simulation that this consistently chooses the largest
#'   \code{lambda} we test.  Therefore this function is not recommended, but we
#'   keep it for organizing our ideas.
#'
#' @param replicate_list A list of frequency count tables
#' @param lamdba_vec A set of lambda values which we propose testing over.
#' @param rtn_full_list logical.  If F we only return the best lambda and the
#'   corresponding estimate.  If T we return all lambdas and estimates in a
#'   dataframe.
#'
#' @export
min_var_C_hat <- function(replicate_list, lambda_vec, rtn_full_list = F, ...) {
  if (is.list(lambda_vec)) {
    lambda_vec <- unlist(lambda_vec)
  }
  if (!is.vector(lambda_vec)) {
    stop("The values of lambda to be tested must be a list or vector")
  }
  lam_df <- data.frame(lambda = lambda_vec,
                       mean_C_hat = rep(NA,length(lambda_vec)),
                       var_C_hat = rep(NA,length(lambda_vec)))
  for(i in 1:length(lambda_vec)) {
    lam_df[i,2:3] <- mean_var_C_hat(replicate_list,
                                    lambda = lambda_vec[i],
                                    penalty = "h1", ...)
  }

  if (rtn_full_list == T) {
    return(lam_df)
  }  else {
    ind <- which(lam_df[,3] == min(lam_df[,3]))[1]
    return(lam_df[ind,]) # So it returns lambda, mean and variance.
  }
}

#' @keywords internal
mean_var_C_hat <- function(replicate_list, ...) {
  C_hats <- get_C_hat_replicates(replicate_list, ...)
  return(c(mean=mean(C_hats),var=var(C_hats)))
}

#' @keywords internal
get_C_hat_replicates <- function(replicate_list, detail = 0, ...) {
  if (!(detail %in% c(0,1,2))) {
    stop("detail parameter must be an integer in {0,1,2}")
  }
  result_list <- lapply(replicate_list,direct_optimise_bisect,...)
  C_hats <- lapply(result_list,(function(li) '$'(li,best) %>%
                                  '['(.,1))) %>% unlist

  if (detail == 0) { # detail = 0 means just give the bestC_vec
    return(C_hats)
    # If detail == 1, returns the table with the optima at each tested C for each replicate.
  } else if (detail >= 1) {
    # this function just extracts the right data frame and adds a column to indicate
    # which replicate this row is associated with.
    full <- Map((function(li,it) '$'(li,full) %>% cbind(.,rep(it,length(li)))),
                result_list, 1:length(result_list)) %>%
      do.call(rbind,.)
    colnames(full)[5] <- "replicate"
    if (detail == 1) {
      list(C_hats = C_hats, full = full)
    } else if (detail == 2) {# If detail == 2, returns the matrix with all the info from every start (potentially very large
      full_starts <- Map((function(li,it) '$'(li,full_starts) %>% cbind(.,rep(it,length(li)))),
                         result_list, 1:length(result_list)) %>%
        do.call(rbind,.)
      colnames(full_starts)[7] <- "replicate"
      list(C_hats = C_hats, full = full, full_starts = full_starts)
    } else {
      # should never execute:
      stop("There is an error in the coding of the detail parameter")
    }
  }
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




