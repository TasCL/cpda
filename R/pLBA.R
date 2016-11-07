#' Calculate Gaussian Log-likelihood and Probability Density
#'
#' This function takes all relevant inputs from the driver and runs a
#' 2-block, differential evolution MCMC
#'
#' @param MCMC_params Provides all MCMC and KDE parameters along with a
#' few other things
#' @param Subject_data Stores the RT data to be fit.
#' @param Model_specifics Stores function handles and block information
#' for the model
#' @param savename Specifies where you want the MCMC output stored.
#' @param ... other arguments
#' @return the chain parameter information (param_chain). This ends up
#' saved in a data file however. The matrix has a structure
#' \code{param_chain(MCMC_iteration,parameter num,chain num)}
#' @param report the iteration interval of returning a progress report.
#' Default 100
#' @keywords DE_MCMC_2block_fun
#' @export
#' @references Holmes, W. (2015). A practical guide to the Probability Density
#' Approximation (PDA) with improved implementation and error characterization.
#' \emph{Journal of Mathematical Psychology}, \bold{68-69}, 13--24,
#' doi: http://dx.doi.org/10.1016/j.jmp.2015.08.006.
#' @importFrom pracma fliplr interp1 std
DE_MCMC_2block_fun <- function(MCMC_params, Subject_data, Model_specifics,
  report=100) {
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## Translated from William R. Holmes's MATLAB codes
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ## Extract some model specifics
  block_1_ind <- Model_specifics$block_1_ind;
  block_2_ind <- Model_specifics$block_2_ind;

  ## These are DE MCMC parameters the determine proposal step size
  gamma_1 <- 2.38/sqrt(2*length(block_1_ind));
  if (!any(is.na(block_2_ind))) { gamma_2 <- 2.38/sqrt(2*length(block_2_ind)) }

  ## Extract MCMC parameters
  Nstep        <- MCMC_params$Nstep;
  Nchain       <- MCMC_params$Nchain;
  noise_size   <- MCMC_params$noise_size;
  burnin       <- MCMC_params$burnin;
  resample_mod <- MCMC_params$resample_mod;

  ## Initialize data structures and stagnation counter
  nparameter <- length(block_1_ind)+length(block_2_ind)
  init <- initialize_structures(nmc=Nstep, npar=nparameter, nchain=Nchain, b=2.7)

  init$direction[,1] <- NA
  direction <- init$direction[,1]
  init$proposal[,1] <- NA

  epsilon <- init$proposal[,1]
  prior_new <- prior_old <- matrix(numeric(Nchain*nparameter),ncol=Nchain)

  proposal <- init$proposal
  param_old <- init$param_old
  LL_keep <- init$LL_keep

  log_lik_old <- matrix(numeric(Nchain),nrow=1)
  log_lik_new <- matrix(numeric(Nchain),nrow=1)

  ## Initialize LL and prior for the first parameter set.
  for(i in 1:Nchain) {
    LL <- Compute_log_likelihood_FFT(Subject_data=Subject_data, params=param_old[,i],
      MCMC_params=MCMC_params);
    log_lik_old[1,i] <- LL;
    prior_old[1,i] <- SwitchModel_Prior(pVec=param_old[,i]);
  }

  LL_keep[1,] <- log_lik_old

  ## Initialize the acceptance counter rate counters.
  acceptance <- 0;
  accept1 <- 0;
  accept2 <- 0;

  param_chain <- array(numeric(Nstep*nparameter*Nchain), dim=c(Nstep,nparameter,Nchain))
  ## Loop over number of chain iterations.
  for(ii in 2:Nstep) {
    ## This is just something to keep track of simulation progress.
    if((ii %% report)==0) { cat(ii," ") }

    ## loop over all chains.
    for(nn in 1:Nchain) {
      ## Parameter block 1
      ## Generate proposal
      rand_int <- sample(Nchain);
      ind1 <- rand_int[-nn][1];
      ind2 <- rand_int[-nn][2];

      ## Note that we are only updating the first block of parameters in
      ## the proposal here.
      direction[block_1_ind] <- param_old[block_1_ind,ind1]-param_old[block_1_ind,ind2];
      direction[block_2_ind] <- 0;

      gamma <- gamma_1;
      epsilon[block_1_ind] <- noise_size*(runif(length(direction[block_1_ind]))-0.5);
      epsilon[block_2_ind] <- 0;

      proposal[,nn] <- param_old[,nn] + gamma * direction + epsilon;

      ## Compute prior
      prior_new[1,nn] <- SwitchModel_Prior(pVec=proposal[,nn]);

      ## Compute likelihood. If the prior is 0, then just assign negative
      ## infinity to ensure it is rejected.
      if(prior_new[1,nn] != 0) {
        LL <- Compute_log_likelihood_FFT(Subject_data,proposal[,nn],MCMC_params);
        log_lik_new[1,nn] <- LL;
      } else {
        log_lik_new[1,nn] <- -Inf;
      }

      ## Compute log acceptance probability
      log_acceptance <- log_lik_new[1,nn] - log_lik_old[1,nn] + log(prior_new[1,nn]) -
        log(prior_old[1,nn]);

      if(is.na(log_acceptance)) {log_acceptance <- -Inf}
      ## Draw a random number to determine if the chain will be accepted.
      prob <- log(runif(length(log_acceptance)));

      ## Update any chain that is accepted
      if(prob < log_acceptance) {

        param_old[,nn] <- proposal[,nn];
        param_chain[ii,,nn] <- proposal[,nn];

        log_lik_old[1,nn] <- log_lik_new[1,nn];
        prior_old[1,nn] <- prior_new[1,nn];

        accept_temp <- 1;
        accept1 <- accept1+1;
      } else {
        ##Leave any chain where the proposal is not accepted.
        param_chain[ii,,nn] <- param_old[,nn];
        accept_temp <- 0;
      }

      ## Pareameter Block 2 - Do everything the same as block 1 but only update the second
      ## block of parameters.
      rand_int <- sample(Nchain);
      ind1 <- rand_int[-nn][1];
      ind2 <- rand_int[-nn][2];

      direction[block_2_ind] <- param_old[block_2_ind,ind1]-param_old[block_2_ind,ind2];
      direction[block_1_ind] <- 0;

      ## Generate proposal
      gamma <- gamma_2;
      epsilon[block_2_ind] <- noise_size*(runif(length(direction[block_2_ind]))-0.5);
      epsilon[block_1_ind] <- 0;

      proposal[,nn] <- param_old[,nn]+gamma*direction+epsilon;

      ## Compute prior
      prior_new[1,nn] <- SwitchModel_Prior(pVec=proposal[,nn]);

      ## Compute likelihood, SPECIFIC to model
      if(prior_new[1,nn] != 0) {
        LL <- Compute_log_likelihood_FFT(Subject_data, proposal[,nn], MCMC_params);
        log_lik_new[1,nn] <- LL;
      } else {
        log_lik_new[1,nn] <- -Inf
      }

      log_acceptance <- log_lik_new[1,nn]-log_lik_old[1,nn] + log(prior_new[1,nn]) -
        log(prior_old[1,nn]);
      if(is.na(log_acceptance)) {log_acceptance <- -Inf}

      prob <- log(runif(length(log_acceptance)));

      ## Update any chain that is accepted
      if(prob < log_acceptance) {
        param_old[,nn] <- proposal[,nn];
        param_chain[ii,,nn]=proposal[,nn];

        log_lik_old[1,nn] <- log_lik_new[1,nn];
        prior_old[1,nn] <- prior_new[1,nn];
        acceptance <- acceptance+1;
        accept2 <- accept2+1;
      } else {
        ## Leave any chain where the proposal is not accepted and increase
        ## the stagnation counter
        param_chain[ii,,nn] <- param_old[,nn];
        if(accept_temp==1) { acceptance <- acceptance+1}
      }
    }

    ## Store the likelihood in the chain data structure.
    LL_keep[ii,] <- log_lik_old;

    ## Resample likelihood at a user defined frequency to unstick chains.
    if(ii %% resample_mod == 0) {
      for(nn in 1:Nchain) {
        ## Compute likelihood, SPECIFIC to model
        params <- param_chain[ii,,nn];
        LL <-  Compute_log_likelihood_FFT(Subject_data, params, MCMC_params);
        log_lik_old[1,nn] <- LL;
      }
    }

    ## Reset outliers after burnin/2
    if(ii == burnin/2) {
      num_param <- dim(param_old)[1];
      num_chain <- dim(param_old)[2];
      t_delay_std <- sd(param_old[6,]);
      t_delay_mean <- mean(param_old[6,]);
      muw1_std <- std(param_old[3,]);
      muw1_mean <- mean(param_old[3,]);
      A_std <- std(param_old[1,]);
      A_mean <- mean(param_old[1,]);

      for(nn in 1:Nchain) {
        if( (abs(param_old[6,nn] - t_delay_mean) > t_delay_std) |
            (abs(param_old[3,nn] - muw1_mean) > muw1_std) |
            (abs(param_old[1,nn] - A_mean) > A_std ) ) {

          temp_param_mat <- param_old;
          temp_param_mat <- temp_param_mat[,-nn]
          for(ij in 1:num_param) {
            temp_vec <- temp_param_mat[ij,];
            replacement <- mean(temp_vec);
            param_old[ij,nn] <- replacement;
            remove(temp_vec)
          }
          param_chain[ii,,nn] <- param_old[,nn];
        }
      }
      acceptance=0;
      accept1=0;
      accept2=0;
    }
  }
  ## Compute some acceptance rates.
  acceptance_rate  <- acceptance / (Nchain*(Nstep-burnin))
  acceptance_rate1 <- accept1 / (Nchain*(Nstep-burnin))
  acceptance_rate2 <- accept2 / (Nchain*(Nstep-burnin))

  # Store simulation infomartion.
  cat("\n")
  return(param_chain)
}

#' @export
DE_MCMC_2block_fun2 <- function(MCMC_params, Subject_data, Model_specifics,
  report=100) {

  ## Extract some model specifics
  block_1_ind <- Model_specifics$block_1_ind;
  block_2_ind <- Model_specifics$block_2_ind;

  ## These are DE MCMC parameters the determine proposal step size
  gamma_1 <- 2.38/sqrt(2*length(block_1_ind));
  if (!any(is.na(block_2_ind))) { gamma_2 <- 2.38/sqrt(2*length(block_2_ind)) }

  ## Extract MCMC parameters
  Nstep        <- MCMC_params$Nstep;
  Nchain       <- MCMC_params$Nchain;
  noise_size   <- MCMC_params$noise_size;
  burnin       <- MCMC_params$burnin;
  resample_mod <- MCMC_params$resample_mod;

  ## Initialize data structures and stagnation counter
  nparameter <- length(block_1_ind)+length(block_2_ind)
  init <- initializeStructures(nmc=Nstep, npar=nparameter, nchain=Nchain)

  init$direction[,1] <- NA
  direction <- init$direction[,1]
  init$proposal[,1] <- NA

  epsilon <- init$proposal[,1]
  prior_new <- prior_old <- matrix(numeric(Nchain*nparameter),ncol=Nchain)

  proposal <- init$proposal
  param_old <- init$param_old
  LL_keep <- init$LL_keep

  log_lik_old <- matrix(numeric(Nchain),nrow=1)
  log_lik_new <- matrix(numeric(Nchain),nrow=1)

  ## Initialize LL and prior for the first parameter set.
  for(i in 1:Nchain) {
    LL <- logLik_pLBA(Subject_data, param_old[,i], MCMC_params);
    log_lik_old[1,i] <- LL;
    prior_old[1,i] <- SwitchModelPrior(pVec=param_old[,i]);
  }

  LL_keep[1,] <- log_lik_old

  ## Initialize the acceptance counter rate counters.
  acceptance <- 0;
  accept1 <- 0;
  accept2 <- 0;

  param_chain <- array(numeric(Nstep*nparameter*Nchain),
    dim=c(Nstep,nparameter,Nchain))
  ## Loop over number of chain iterations.
  for(ii in 2:Nstep) {
    ## This is just something to keep track of simulation progress.
    if((ii %% report)==0) { cat(ii," ") }

    ## loop over all chains.
    for(nn in 1:Nchain) {
      ## Parameter block 1
      ## Generate proposal
      rand_int <- sample(Nchain);
      ind1 <- rand_int[-nn][1];
      ind2 <- rand_int[-nn][2];

      ## Note that we are only updating the first block of parameters in
      ## the proposal here.
      direction[block_1_ind] <- param_old[block_1_ind,ind1]-
        param_old[block_1_ind,ind2];
      direction[block_2_ind] <- 0;

      gamma <- gamma_1;
      epsilon[block_1_ind] <- noise_size*(runif(length(direction[block_1_ind]))-0.5);
      epsilon[block_2_ind] <- 0;

      proposal[,nn] <- param_old[,nn] + gamma * direction + epsilon;

      ## Compute prior
      prior_new[1,nn] <- SwitchModelPrior(pVec=proposal[,nn]);

      ## Compute likelihood. If the prior is 0, then just assign negative
      ## infinity to ensure it is rejected.
      if(prior_new[1,nn] != 0) {
        LL <- logLik_pLBA(Subject_data, proposal[,nn], MCMC_params);
        log_lik_new[1,nn] <- LL;
      } else {
        log_lik_new[1,nn] <- -Inf;
      }

      ## Compute log acceptance probability
      log_acceptance <- log_lik_new[1,nn] - log_lik_old[1,nn] +
        log(prior_new[1,nn]) - log(prior_old[1,nn]);

      if(is.na(log_acceptance)) {log_acceptance <- -Inf}
      ## Draw a random number to determine if the chain will be accepted.
      prob <- log(runif(length(log_acceptance)));

      ## Update any chain that is accepted
      if(prob < log_acceptance) {

        param_old[,nn] <- proposal[,nn];
        param_chain[ii,,nn] <- proposal[,nn];

        log_lik_old[1,nn] <- log_lik_new[1,nn];
        prior_old[1,nn] <- prior_new[1,nn];

        accept_temp <- 1;
        accept1 <- accept1+1;
      } else {
        ##Leave any chain where the proposal is not accepted.
        param_chain[ii,,nn] <- param_old[,nn];
        accept_temp <- 0;
      }

      ## Pareameter Block 2 - Do everything the same as block 1 but only
      ## update the second block of parameters.
      rand_int <- sample(Nchain);
      ind1 <- rand_int[-nn][1];
      ind2 <- rand_int[-nn][2];

      direction[block_2_ind] <- param_old[block_2_ind,ind1]-param_old[block_2_ind,ind2];
      direction[block_1_ind] <- 0;

      ## Generate proposal
      gamma <- gamma_2;
      epsilon[block_2_ind] <- noise_size*(runif(length(direction[block_2_ind]))-0.5);
      epsilon[block_1_ind] <- 0;

      proposal[,nn] <- param_old[,nn]+gamma*direction+epsilon;

      ## Compute prior
      prior_new[1,nn] <- SwitchModelPrior(pVec=proposal[,nn]);

      ## Compute likelihood, SPECIFIC to model
      if(prior_new[1,nn] != 0) {
        LL <- logLik_pLBA(Subject_data, proposal[,nn], MCMC_params);
        log_lik_new[1,nn] <- LL;
      } else {
        log_lik_new[1,nn] <- -Inf
      }

      log_acceptance <- log_lik_new[1,nn]-log_lik_old[1,nn] +
        log(prior_new[1,nn]) - log(prior_old[1,nn]);
      if(is.na(log_acceptance)) {log_acceptance <- -Inf}

      prob <- log(runif(length(log_acceptance)));

      ## Update any chain that is accepted
      if(prob < log_acceptance) {
        param_old[,nn] <- proposal[,nn];
        param_chain[ii,,nn]=proposal[,nn];

        log_lik_old[1,nn] <- log_lik_new[1,nn];
        prior_old[1,nn] <- prior_new[1,nn];
        acceptance <- acceptance+1;
        accept2 <- accept2+1;
      } else {
        ## Leave any chain where the proposal is not accepted and increase
        ## the stagnation counter
        param_chain[ii,,nn] <- param_old[,nn];
        if(accept_temp==1) { acceptance <- acceptance+1}
      }
    }

    ## Store the likelihood in the chain data structure.
    LL_keep[ii,] <- log_lik_old;

    ## Resample likelihood at a user defined frequency to unstick chains.
    if(ii %% resample_mod == 0) {
      for(nn in 1:Nchain) {
        ## Compute likelihood, SPECIFIC to model
        params <- param_chain[ii,,nn];
        LL <-  logLik_pLBA(Subject_data, params, MCMC_params);
        log_lik_old[1,nn] <- LL;
      }
    }

    ## Reset outliers after burnin/2
    if(ii == burnin/2) {
      num_param <- dim(param_old)[1];
      num_chain <- dim(param_old)[2];
      t_delay_std <- sd(param_old[6,]);
      t_delay_mean <- mean(param_old[6,]);
      muw1_std <- std(param_old[3,]);
      muw1_mean <- mean(param_old[3,]);
      A_std <- std(param_old[1,]);
      A_mean <- mean(param_old[1,]);

      for(nn in 1:Nchain) {
        if( (abs(param_old[6,nn] - t_delay_mean) > t_delay_std) |
            (abs(param_old[3,nn] - muw1_mean) > muw1_std) |
            (abs(param_old[1,nn] - A_mean) > A_std ) ) {

          temp_param_mat <- param_old;
          temp_param_mat <- temp_param_mat[,-nn]
          for(ij in 1:num_param) {
            temp_vec <- temp_param_mat[ij,];
            replacement <- mean(temp_vec);
            param_old[ij,nn] <- replacement;
            remove(temp_vec)
          }
          param_chain[ii,,nn] <- param_old[,nn];
        }
      }
      acceptance=0;
      accept1=0;
      accept2=0;
    }
  }
  ## Compute soem acceptance rates.
  acceptance_rate  <- acceptance / (Nchain*(Nstep-burnin))
  acceptance_rate1 <- accept1 / (Nchain*(Nstep-burnin))
  acceptance_rate2 <- accept2 / (Nchain*(Nstep-burnin))
  cat("\n")
  return(param_chain)
}

#' Initialize a Structure
#'
#' This file initializes the initial conditions and various data structures
#' for the MCMC
#'
#' @param nmc number of MCMC iteration (steps). DMC's nmc. Default 100
#' @param npar number of parameter. Default 2
#' @param nchain number of chains. Default 3
#' @param b LBA's b parameter. Default 2.7
#' @keywords initialize_structures
#' @return a list with 5 elements: param_old, param_chain, proposal, direction,
#' and LL_keep
#' @export
#' @examples
#' nchain <- 3
#' nmc    <- 100
#' npar   <- 7
#' initialize_structures(nmc, npar, nchain)
initialize_structures <- function(nmc=100, npar=7, nchain=3, b=2.7) {
  ## initializeStructures(nmc=20, npar=7, nchain=3)
  param_old <- matrix(numeric(npar*nchain), nrow=npar)
  param_old[1,] <- 0 +  b*rnorm(nchain);  # A
  param_old[2,] <- 0 +  5*runif(nchain);               # muv1
  param_old[3,] <- 0 +  5*runif(nchain);               # muw1
  param_old[4,] <- 0 +  5*runif(nchain);               # muv2
  param_old[5,] <- 0 +  5*runif(nchain);               # muw2
  param_old[6,] <- 0 + .5*runif(nchain);              # t_delay
  param_old[7,] <- 0 + .5*runif(nchain);              # t_nd

  ## This will store the MCMC chain history.
  param_chain <- array(numeric(nmc*npar*nchain),
                       dim = c(nmc, npar, nchain))

  ## Populate first value of the chain with initial guess
  param_chain[1,1,] <- param_old[1,];
  param_chain[1,2,] <- param_old[2,];
  param_chain[1,3,] <- param_old[3,];
  param_chain[1,4,] <- param_old[4,];
  param_chain[1,5,] <- param_old[5,];
  param_chain[1,6,] <- param_old[6,];
  param_chain[1,7,] <- param_old[7,];

  ## This is a temporary structure to hold a proposal.
  proposal <- direction <- matrix(numeric(npar*nchain), nrow=npar)
  LL_keep  <- matrix(numeric(nmc*nchain), nrow=nmc)

  out <- list(param_old,param_chain,proposal,direction,LL_keep)
  names(out) <- c("param_old","param_chain","proposal","direction","LL_keep")
  return(out)
}

#' Switch Model's Prior Distributions
#'
#' This function computes the prior for the given parameter set. The prior
#' for each parameter here is just a uniform distribution for simplicity.
#'
#' @param pVec a parameter vector. Similar to DMC's p.vector
#' @return a prior probability density summed across parameters
#' @export
#' @examples
#' p.vector <- c(1.51, 3.32, 1.51, 2.24, 3.69, 0.31, 0.08)
#' SwitchModel_Prior(pVec=p.vector)
SwitchModel_Prior <- function(pVec) {
  ## b    A       muv1    muw1    muv2    muw2    t_delay t_ND
  ## 2.7  pVec[1] pVec[2] pVec[3] pVec[4] pVec[5] pVec[6] pVec[7]
  b       <- 2.7;
  A       <- pVec[1];
  muv1    <- pVec[2];
  muw1    <- pVec[3];
  muv2    <- pVec[4];
  muw2    <- pVec[5];
  t_delay <- pVec[6];
  t_ND    <- pVec[7];

  if(A < 0 | A > 10) { pA <- 0 } else { pA <- 1/10 }
  if(b < 0 | b > 10) { pb <- 0 } else { pb <- 1/10 }
  if(muv1 < -3 | muv1>7) { pV1 <- 0 } else { pV1 <- 1/10 }
  if(muv2 < -3 | muv2>7) { pV2 <- 0 } else { pV2 <- 1/10 }
  if(muw1 < -3 | muw1 > 7) { pW1 <- 0 } else { pW1 <- 1/10 }
  if(muw2 < -3 | muw2 > 7) { pW2 <- 0 } else { pW2 <- 1/10 }
  if(t_delay < 0 | t_delay > 1) { pT <- 0 } else { pT <- 1/1}
  if(t_ND < 0 | t_ND > 1) { pND <- 0} else { pND <- 1/1}

  prior <- pA*pb*pT*pV1*pV2*pW1*pW2*pND;
  return(prior)
}

#' Compute Log-likelihood Using Fast Fourier Transform
#'
#' This function takes subject data, the parameters for the current chain
#' under consideration, and MCMC parameters and computes the log likelihood.
#'
#' @param Subject_data subject data
#' @param params LBA parameters
#' @param MCMC_params MCMC parameters
#' @return Log-likelihood
#' @export
Compute_log_likelihood_FFT <- function(Subject_data, params, MCMC_params) {
  ## This function takes subject data, the parameters for the current chain
  ## under consideration, and MCMC parameters and computes the log likelihood.

  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## Translated from William R. Holmes's MATLAB codes
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ## Extract the fixed parameter sigma
  sigma <- MCMC_params$sigma_exact;
  ## Extract the KDE parameters
  bandwidth <- MCMC_params$bandwidth;
  Nsample <- MCMC_params$LL_NSAMPLE;

  ## Extract model parameters.
  A      <- params[1];
  muv1   <- params[2];
  muw1   <- params[3];
  muv2   <- params[4];
  muw2   <- params[5];
  t_delay<- params[6];
  t_ND   <- params[7];

  ## Extract data
  T0 <- as.vector(attr(d, "SwitchTime"))
  C1time <- attr(d,"C1time") - t_ND ## Subtract off the non decision time
  C2time <- attr(d,"C2time") - t_ND
  ##  Add the delay to the switch time.
  T0 <- T0 + t_delay;

  ## Do sampling. Run LL_NSAMPLE sets of accumulators forward in time and see which one terminates first.
  x1_0 <- A*runif(Nsample) ##  Uniform, start point
  x2_0 <- A*runif(Nsample) ##  Uniform, start point

  v1 <- muv1+sigma*rnorm(Nsample) ## Normal, rate
  w1 <- muw1+sigma*rnorm(Nsample) ## Normal, rate

  v2 <- muv2+sigma*rnorm(Nsample) ## Normal, rate
  w2 <- muw2+sigma*rnorm(Nsample) ## Normal, rate

  ## Compute the time at which the accumulators under the influence of muv1 and
  ## muv2 will terminate.
  T1_s_first <- (2.7-x1_0) / (v1);
  T2_s_first <- (2.7-x2_0) / (v2);

  ## Handle negative drift rates
  T1_s_first[v1 < 0] <- 1000000;
  T2_s_first[v2 < 0] <- 1000000;

  ## Compute the time at which the accumulators under the subsequent influence
  ## of muw1 and muw1 will terminate.
  T1_s_second <- T0+(2.7-x1_0-v1*T0) / w1;
  T2_s_second <- T0+(2.7-x2_0-v2*T0) / w2;

  ## Handle negative drift rates
  T1_s_second[w1 < 0] <- 1000000;
  T2_s_second[w2 < 0] <- 1000000;

  ## Separate out indices for pre switch termination and post switch
  ## termination
  ind_C1_before <- (T1_s_first < T0) & (T1_s_first < T2_s_first);
  ind_C2_before <- (T2_s_first < T0) & (T2_s_first < T1_s_first);

  ind_C1_after <- (T1_s_first > T0) & (T2_s_first > T0) & (T1_s_second<T2_s_second) &
    (T1_s_second < 10000);
  ind_C2_after <- (T1_s_first > T0) & (T2_s_first > T0) & (T2_s_second<T1_s_second) &
    (T2_s_second < 10000);

  RT_s1_after <- c(T1_s_first[ind_C1_before], T1_s_second[ind_C1_after]);
  RT_s2_after <- c(T2_s_first[ind_C2_before], T2_s_second[ind_C2_after]);

  ## Do FFT Smoothing
  ## Extract the minimum and masimum switch times.
  m <- min(min(C1time),min(C2time))-3*bandwidth;
  M <- max(max(C1time),max(C2time))+3*bandwidth;

  ## Set the histogram grid size and construct the bin centers / edges
  N_grid <- 2^10;
  grid_centers <- seq(m, M, length.out=N_grid);

  dt <- grid_centers[2]-grid_centers[1];
  term1 <- grid_centers - dt/2
  term2 <- grid_centers[length(grid_centers)] + dt/2
  bin_edges <- c(term1,term2);

  ## Initialize the Gaussian filter on the spectral domain.
  freq <- 2*pi * N_grid/(M-m) * 0.5 * seq(0,1, length.out=N_grid/2+1);

  fil0 <- exp(-0.5 * freq^2 * bandwidth ^ 2);

  ## require(pracma)
  fil1 <- fliplr( t(as.matrix(fil0[2:(length(fil0)-1)])) )
  filter <- c(fil0, fil1)

  ## Compute the LL contribution from the option 1 data.
  LL <- 0;
  if(any(!is.na(C1time))) {
    ## a solution to the inconsistency between histc and hist
    if(min(term1) > min(RT_s1_after)) {bin_edges[1] <- min(RT_s1_after)}
    if(max(term2) < max(RT_s1_after)) {bin_edges[length(bin_edges)] <- max(RT_s1_after)}

    ## Bin the data to create the noisy histogram.
    tmp <- hist(RT_s1_after, breaks=bin_edges, plot = FALSE, right = TRUE)
    tmp$counts[1] <- 0
    tmp$counts[length(tmp$counts)] <- 0
    bincount <- c(tmp$counts, 0);

    PDF1_hist <- (1 * bincount) / (dt * Nsample);

    ## Apply the FFT
    PDF1_fft <- fft(PDF1_hist[1:(length(PDF1_hist))-1])

    ## Multiply the data by the filter in the spectral domain.
    PDF1_fft_filt <- filter * PDF1_fft;

    ## Transform back to data space.
    PDF1_smoothed <- fft(PDF1_fft_filt, inverse=TRUE) / length(PDF1_fft_filt);

    ## Interpolate the grid likelihood to the data
    PDF1 <- interp1(grid_centers, as.vector(Re(PDF1_smoothed)), C1time, method="linear");

    ## Reset any zero values to a minimum value.
    PDF1 <- pmax(PDF1, 10^(-5));

    ## Compute the log likelihood
    LL <- LL + sum(log(PDF1));
  }


  ## Compute the LL contribution from the option 2 data.
  if (any(!is.na(C2time))) {
    ## a solution to the inconsistency between histc and hist
    if(min(term1) > min(RT_s2_after)) {bin_edges[1] <- min(RT_s2_after)}
    if(max(term2) < max(RT_s2_after)) {bin_edges[length(bin_edges)] <- max(RT_s2_after)}

    tmp <- hist(RT_s2_after, breaks=bin_edges, plot = FALSE, right = TRUE)
    tmp$counts[1] <- 0
    tmp$counts[length(tmp$counts)] <- 0
    bincount <- c(tmp$counts, 0);

    PDF2_hist <- (1 * bincount) / (dt * Nsample);

    ## Apply the FFT
    PDF2_fft <- fft(PDF2_hist[1:(length(PDF2_hist))-1])

    ## Multiply the data by the filter in the spectral domain.
    PDF2_fft_filt <- filter * PDF2_fft;

    ## Transform back to data space.
    PDF2_smoothed <- fft(PDF2_fft_filt, inverse=TRUE) / length(PDF2_fft_filt);

    ## Interpolate the grid likelihood to the data
    PDF2 <- interp1(grid_centers, as.vector(Re(PDF2_smoothed)), C2time, method="linear");

    ## Reset any zero values to a minimum value.
    PDF2 <- pmax(PDF2, 10^(-20));

    ## Compute the log likelihood
    LL <- LL + sum(log(PDF2));
  }
  rm(filter)
  return(LL)
}


#' Compute Kernel Bandwidth
#'
#' Use Silverman's run of thumb to calculate kernel bandwidth.
#'
#' @param object a data vector s
#' @param ns number of simulation
#' @param ... other arguments
#' @return a scalar value of bandwidth
#' @export
#' @examples
#' datavec <- rnorm(1000, 5, 1)
#' get_bandwidth(datavec, 10000)
get_bandwidth <- function(object, ns, ...) {
  iqr <- IQR(object, type=1)
  out <- .9 * min(iqr, sd(object)) * ns^(-.2)
  return(out)
}
