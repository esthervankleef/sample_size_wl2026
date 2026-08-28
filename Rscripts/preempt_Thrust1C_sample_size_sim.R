##########################################################################
# PRE-EMPT SAMPLE-SIZE / PRECISION SIMULATION
# Household continuous-time Markov model (CTMC)
#
# Date created: 16 August 2026
# Date last updated: 19 August 2026
#
# PURPOSE
# Evaluate candidate PRE-EMPT cohort sizes for Thrusts 1A, 1B and 1C.
#
# THRUST 1A: CTMC PARAMETER PRECISION
# Household pARG carriage is modelled as a continuous-time process:
#
#   pARG-negative -> pARG-positive   (acquisition)
#   pARG-positive -> pARG-negative   (clearance)
#
# Acquisition hazard for a susceptible household member:
#
#   background_rate + n_positive_contacts * p2p_rate
#
# Clearance occurs at rate gamma.
#
# Each simulated study is refitted using the same CTMC structure to
# estimate background acquisition, beta (person-to-person transmission)
# and gamma (clearance). Precision is assessed using log-scale Wald 95%
# CIs, with a design target of +/-30% for each parameter.
#
# THRUST 1B: INFORMATIVE INDEX YIELD
# Thrust 1B is treated here as a cohort-yield calculation, NOT as a power
# calculation for R0-within. An informative index is a participant who is
# pARG-positive at the time of a captured priority-antibiotic exposure.
# These participants are eligible for intensive sampling at days
# 0, 3, 7, 14 and 28.
#
# The simulation estimates the expected index yield and probability of
# obtaining >=40, 60, 80 or 100 informative episodes. These are design
# sensitivity thresholds, not biological requirements for R0-within.
#
# THRUST 1C: R0-BETWEEN PRECISION
# Household R0-between is approximated as:
#
#   R0-between = (hhsize - 1) * [1 - exp(-beta/gamma)]
#
# where 1/gamma is mean carriage duration and
# 1-exp(-beta/gamma) is the probability of transmission to one susceptible
# household contact during a carriage episode.
#
# R0-between and its approximate 95% CI are calculated from each fitted
# CTMC using the delta method on the log scale. The primary precision
# target is +/-50%, with +/-25% considered as a stricter sensitivity.
#
# IMPORTANT
# This is a design-stage simulation, not the final PRE-EMPT analysis.
# Results are conditional on the assumed prevalence, acquisition,
# clearance, household transmission, antibiotic-use and exposure-capture
# parameters and should be interpreted alongside sensitivity analyses.
##########################################################################

rm(list=ls())

set.seed(12)

library(pacman)
pacman::p_load(msm,writexl)

##############################
# STUDY DESIGN
##############################

# Site currently being explored.
site <- "Kenya"

# Number of people per household.
hhsize <- 5

# Total follow-up.
follow_up_days <- 365

# Routine sampling every 3 months.
routine_times <- c(0,90,180,270,365)

# Intensive sampling after a priority-antibiotic exposure.
index_times <- c(0,3,7,14,28)

# Candidate cohort sizes in households.
household_grid <- c(100, 150, 200)

# Number of simulated studies per cohort size.
# Increase further for final analysis if needed.
numsims <- 100

##############################
# THRUST 1B: INDEX-YIELD SENSITIVITY
##############################

# These are not formal biological requirements.
# They show the probability of generating at least 40,60,80 or 100
# pARG-positive participants receiving a priority antibiotic.
index_thresholds <- c(40,60,80,100)

##############################
# PRECISION TARGETS
##############################

# Primary R0-between criterion:
# 95% CI entirely within +/-50% of the point estimate.
primary_R0_precision_target <- 0.50

# Optional stricter R0 sensitivity criterion:
# 95% CI entirely within +/-25% of the point estimate.
secondary_R0_precision_target <- 0.25

# Thrust 1A parameter-level criterion:
# 95% CI entirely within +/-30% of the point estimate.
parameter_precision_target <- 0.30

##############################
# BURKINA-INSPIRED HIGH-TRANSMISSION SCENARIO
##############################

# Design-stage baseline pARG prevalence.
baseline_prevalence <- 0.50

# Overall daily acquisition hazard used as the calibration target.
overall_acquisition_rate <- 0.0135

# Mean carriage duration and corresponding clearance rate.
mean_carriage_days <- 102
clearance_rate <- 1/mean_carriage_days

# Fraction of overall acquisition attributed to background/environmental
# acquisition rather than household person-to-person transmission.
background_fraction <- 0.30

# Calibrated background acquisition rate.
background_rate <- overall_acquisition_rate*background_fraction

# Calibrated household person-to-person transmission rate.
#
# overall acquisition =
# background acquisition +
# household acquisition
#
# household acquisition is approximated as:
# (hhsize-1) * baseline_prevalence * p2p_rate
p2p_rate <- (overall_acquisition_rate-background_rate)/
  ((hhsize-1)*baseline_prevalence)

##############################
# ANTIBIOTIC PARAMETERS
##############################

# Average number of antibiotic courses per person-year.
abx_courses_per_year <- 0.80

# Proportion of antibiotic courses that are priority antibiotics.
priority_abx_fraction <- 0.20

# Probability that a relevant antibiotic exposure is captured.
abx_capture_probability <- 0.90

##############################
# HOUSEHOLD STATE SPACE
##############################

# Each household has 2^hhsize possible pARG-positive/negative states.
nstates <- 2^hhsize

state_bits <- matrix(0,nrow=nstates,ncol=hhsize)

for(s in 1:nstates){
  state_bits[s,] <- as.integer(
    intToBits(s-1)
  )[1:hhsize]
}

# Convert individual household states into the CTMC state number.
state_number <- function(x){
  1+sum(
    x*2^(0:(hhsize-1))
  )
}

# Identify allowable acquisition and clearance transitions.
add_matrix <- matrix(FALSE,nstates,nstates)
remove_matrix <- matrix(FALSE,nstates,nstates)

for(i in 1:nstates){
  for(j in 1:nstates){
    if(i==j) next
    x <- state_bits[i,]
    y <- state_bits[j,]
    if(sum(x!=y)==1 && sum(y-x)==1){
      add_matrix[i,j] <- TRUE
    }
    if(sum(x!=y)==1 && sum(x-y)==1){
      remove_matrix[i,j] <- TRUE
    }
  }
}

##############################
# ANTIBIOTIC EXPOSURE
##############################

# Simulate priority-antibiotic exposure dates for one participant.
make_abx_history <- function(){
  daily_p <- 1-exp(
    -abx_courses_per_year/365
  )
  days <- integer(0)
  for(day in 1:follow_up_days){
    if(
      runif(1)<daily_p &&
      runif(1)<=abx_capture_probability &&
      runif(1)<priority_abx_fraction
    ){
      days <- c(days,day)
    }
  }
  days
}

##############################
# INITIAL HOUSEHOLD STATE
##############################

# Simulate baseline pARG carriage independently for household members.
sample_initial_state <- function(){
  infected <- rbinom(
    hhsize,
    1,
    baseline_prevalence
  )
  state_number(infected)
}

##############################
# TRUE HOUSEHOLD CTMC
##############################

# Acquisition rate:
# background_rate + number of positive household contacts * p2p_rate
#
# Clearance rate:
# gamma
make_Q <- function(bg,beta,gamma){
  Q <- matrix(
    0,
    nstates,
    nstates
  )
  for(i in 1:nstates){
    n_positive <- sum(
      state_bits[i,]
    )
    Q[
      i,
      remove_matrix[i,]
    ] <- gamma
    Q[
      i,
      add_matrix[i,]
    ] <- bg+n_positive*beta
    Q[i,i] <- -sum(Q[i,])
  }
  Q
}

Q_true <- make_Q(
  background_rate,
  p2p_rate,
  clearance_rate
)

##############################
# CACHE TRUE TRANSITION MATRICES
##############################

# MatrixExp() is expensive, so calculate each unique time interval once.
P_cache <- list()

get_P_true <- function(dt){
  key <- as.character(dt)
  if(is.null(P_cache[[key]])){
    P_cache[[key]] <<- MatrixExp(
      Q_true,
      t=dt
    )
  }
  P_cache[[key]]
}

##############################
# SIMULATE ONE HOUSEHOLD
##############################

# Simulate only at routine sampling dates and antibiotic exposure dates,
# rather than every day, to keep the simulation efficient.
simhh <- function(){
  abx_history <- vector(
    "list",
    hhsize
  )
  
  for(i in 1:hhsize){
    abx_history[[i]] <- make_abx_history()
  }
  
  abx_days <- unique(
    unlist(abx_history)
  )
  
  times <- sort(
    unique(
      c(
        routine_times,
        abx_days
      )
    )
  )
  
  states <- integer(
    length(times)
  )
  
  state <- sample_initial_state()
  states[1] <- state
  
  if(length(times)>1){
    for(k in 2:length(times)){
      dt <- times[k]-times[k-1]
      P <- get_P_true(dt)
      state <- sample(
        1:nstates,
        1,
        prob=P[state,]
      )
      states[k] <- state
    }
  }
  
  list(
    times=times,
    states=states,
    abx_history=abx_history
  )
}

##############################
# SIMULATE COHORT
##############################

makefakehhdata <- function(numhh){
  households <- vector(
    "list",
    numhh
  )
  for(h in 1:numhh){
    households[[h]] <- simhh()
  }
  households
}

##############################
# GET STATE AT A SPECIFIC DAY
##############################

get_state <- function(household,day){
  household$states[
    match(
      day,
      household$times
    )
  ]
}

##############################
# THRUST 1B: INFORMATIVE INDICES
##############################

# Identify every pARG-positive priority-antibiotic exposure episode.
# These episodes define the potential intensive-follow-up population.
get_index_episodes <- function(households){
  output <- list()
  row_id <- 1
  
  for(h in seq_along(households)){
    hh <- households[[h]]
    
    for(index in 1:hhsize){
      abx_days <- hh$abx_history[[index]]
      
      if(length(abx_days)==0){
        next
      }
      
      for(abx_day in abx_days){
        index_state <- get_state(
          hh,
          abx_day
        )
        
        if(
          state_bits[
            index_state,
            index
          ]!=1
        ){
          next
        }
        
        output[[row_id]] <- data.frame(
          household=h,
          index=index,
          antibiotic_day=abx_day
        )
        
        row_id <- row_id+1
      }
    }
  }
  
  if(length(output)==0){
    return(data.frame())
  }
  
  do.call(
    rbind,
    output
  )
}

##############################
# THRUST 1C: INFORMATIVE PAIRS
##############################

# An informative pair requires:
# 1. pARG-positive index at antibiotic exposure
# 2. pARG-negative contact at that time
# 3. a subsequent routine sample
#
# Acquisition is observed as an interval-censored event between the
# antibiotic exposure and the next routine sampling time.
get_pairs <- function(households){
  pair_list <- list()
  row_id <- 1
  
  for(h in seq_along(households)){
    hh <- households[[h]]
    
    for(index in 1:hhsize){
      abx_days <- hh$abx_history[[index]]
      
      if(length(abx_days)==0){
        next
      }
      
      for(abx_day in abx_days){
        index_state <- get_state(
          hh,
          abx_day
        )
        
        if(
          state_bits[
            index_state,
            index
          ]!=1
        ){
          next
        }
        
        future_times <- routine_times[
          routine_times>abx_day
        ]
        
        if(length(future_times)==0){
          next
        }
        
        followup_day <- future_times[1]
        
        followup_state <- get_state(
          hh,
          followup_day
        )
        
        for(contact in setdiff(
          1:hhsize,
          index
        )){
          if(
            state_bits[
              index_state,
              contact
            ]!=0
          ){
            next
          }
          
          pair_list[[row_id]] <- data.frame(
            household=h,
            index=index,
            contact=contact,
            antibiotic_day=abx_day,
            followup_day=followup_day,
            acquired=as.integer(
              state_bits[
                followup_state,
                contact
              ]==1
            )
          )
          
          row_id <- row_id+1
        }
      }
    }
  }
  
  if(length(pair_list)==0){
    return(data.frame())
  }
  
  do.call(
    rbind,
    pair_list
  )
}

##############################
# ROUTINE CTMC TRANSITION COUNTS
##############################

# Aggregate identical household state transitions.
# This does not change the likelihood; it simply avoids repeatedly
# evaluating identical transition probabilities.
get_transition_counts <- function(households){
  counts <- list()
  
  for(h in seq_along(households)){
    hh <- households[[h]]
    
    states <- sapply(
      routine_times,
      function(t){
        get_state(
          hh,
          t
        )
      }
    )
    
    for(k in 2:length(routine_times)){
      dt <- routine_times[k]-routine_times[k-1]
      
      key <- paste(
        dt,
        states[k-1],
        states[k],
        sep="_"
      )
      
      if(is.null(counts[[key]])){
        counts[[key]] <- 1
      }else{
        counts[[key]] <- counts[[key]]+1
      }
    }
  }
  
  keys <- names(counts)
  
  out <- matrix(
    0,
    nrow=length(keys),
    ncol=4
  )
  
  for(i in seq_along(keys)){
    z <- as.numeric(
      strsplit(
        keys[i],
        "_",
        fixed=TRUE
      )[[1]]
    )
    
    out[i,] <- c(
      z[1],
      z[2],
      z[3],
      counts[[keys[i]]]
    )
  }
  
  data.frame(
    dt=out[,1],
    from=out[,2],
    to=out[,3],
    n=out[,4]
  )
}

##############################
# FITTED HOUSEHOLD CTMC
##############################

# Same model as the data-generating CTMC, but the rates are estimated
# from each simulated study.
make_Q_fit <- function(bg,beta,gamma){
  Q <- matrix(
    0,
    nstates,
    nstates
  )
  
  for(i in 1:nstates){
    n_positive <- sum(
      state_bits[i,]
    )
    
    Q[
      i,
      remove_matrix[i,]
    ] <- gamma
    
    Q[
      i,
      add_matrix[i,]
    ] <- bg+n_positive*beta
    
    Q[i,i] <- -sum(
      Q[i,]
    )
  }
  
  Q
}

##############################
# FAST MARKOV LIKELIHOOD
##############################

# Routine sampling has intervals of 90,90,90 and 95 days.
# Calculate a transition matrix for each unique interval length.
markov_nll <- function(par,data){
  bg <- exp(
    par[1]
  )
  
  beta <- exp(
    par[2]
  )
  
  gamma <- exp(
    par[3]
  )
  
  Q <- make_Q_fit(
    bg,
    beta,
    gamma
  )
  
  interval_times <- unique(
    data$dt
  )
  
  P_list <- lapply(
    interval_times,
    function(x){
      MatrixExp(
        Q,
        t=x
      )
    }
  )
  
  names(P_list) <- as.character(
    interval_times
  )
  
  ll <- 0
  
  for(i in 1:nrow(data)){
    P <- P_list[[as.character(data$dt[i])]]
    
    p <- P[
      data$from[i],
      data$to[i]
    ]
    
    ll <- ll+
      data$n[i]*
      log(p+1e-12)
  }
  
  -ll
}

##############################
# FIT CTMC
##############################

fit_markov <- function(households){
  data <- get_transition_counts(
    households
  )
  
  fit <- try(
    optim(
      par=log(
        c(
          background_rate,
          p2p_rate,
          clearance_rate
        )
      ),
      fn=markov_nll,
      data=data,
      method="BFGS",
      hessian=TRUE,
      control=list(
        maxit=100
      )
    ),
    silent=TRUE
  )
  
  if(
    inherits(
      fit,
      "try-error"
    )
  ){
    return(NULL)
  }
  
  if(
    fit$convergence!=0
  ){
    return(NULL)
  }
  
  V <- try(
    solve(
      fit$hessian
    ),
    silent=TRUE
  )
  
  if(
    inherits(
      V,
      "try-error"
    )
  ){
    return(NULL)
  }
  
  # Do not impose parameter bounds.
  # Reject only clearly unusable covariance matrices.
  if(
    any(!is.finite(V)) ||
    any(diag(V)<=0)
  ){
    return(NULL)
  }
  
  list(
    estimates=exp(fit$par),
    vcov=V
  )
}

##############################
# R0-BETWEEN
##############################

# Current design-stage household R0 approximation:
# (number of household contacts) *
# probability of transmission to one contact during carriage.
calculate_R0 <- function(beta,gamma){
  duration <- 1/gamma
  p_contact <- 1-exp(
    -beta*duration
  )
  (hhsize-1)*p_contact
}

##############################
# R0 CONFIDENCE INTERVAL
##############################

# Calculate an approximate 95% CI for R0 using the delta method
# on the log scale.
#
# R0 = (hhsize-1)*(1-exp(-beta/gamma))
calculate_R0_ci <- function(fit){
  beta <- fit$estimates[2]
  gamma <- fit$estimates[3]
  
  R0 <- calculate_R0(
    beta,
    gamma
  )
  
  x <- beta/gamma
  
  # Derivative of log(R0) with respect to log(beta) and log(gamma).
  g <- if(
    abs(x)<1e-8
  ){
    1
  }else{
    x*exp(-x)/(1-exp(-x))
  }
  
  var_log_R0 <- g^2*(
    fit$vcov[2,2]+
      fit$vcov[3,3]-
      2*fit$vcov[2,3]
  )
  
  se_log_R0 <- sqrt(
    max(
      var_log_R0,
      0
    )
  )
  
  c(
    R0=R0,
    lower=exp(
      log(R0)-1.96*se_log_R0
    ),
    upper=exp(
      log(R0)+1.96*se_log_R0
    )
  )
}

##############################
# THRUST 1A: PARAMETER PRECISION
##############################

# Calculate log-scale Wald 95% CIs for the three fitted CTMC rates.
#
# The rates themselves need not be normally distributed.
# The approximation is that their log estimates are approximately normal.
#
# We then apply the +/-30% criterion directly:
# lower >= 70% of estimate AND upper <= 130% of estimate.
calculate_param_precision <- function(fit){
  se_log <- sqrt(
    diag(fit$vcov)
  )
  
  est <- fit$estimates
  
  lower <- exp(
    log(est)-1.96*se_log
  )
  
  upper <- exp(
    log(est)+1.96*se_log
  )
  
  within_30 <- (
    lower>=0.70*est &
      upper<=1.30*est
  )
  
  names(lower) <- c(
    "bg",
    "beta",
    "gamma"
  )
  
  names(upper) <- c(
    "bg",
    "beta",
    "gamma"
  )
  
  names(within_30) <- c(
    "bg",
    "beta",
    "gamma"
  )
  
  list(
    lower=lower,
    upper=upper,
    within_30=within_30
  )
}

##############################
# ONE SIMULATED STUDY
##############################

run_one_study <- function(n_households){
  
  households <- makefakehhdata(
    n_households
  )
  
  ##############################
  # ANTIBIOTIC EXPOSURE BURDEN
  ##############################
  
  # Number of individuals with at least one captured priority-antibiotic
  # exposure.
  abx_exposed_individuals <- 0
  
  # Total number of captured priority-antibiotic exposure episodes.
  abx_exposure_episodes <- 0
  
  for(hh in households){
    for(i in 1:hhsize){
      abx_days <- hh$abx_history[[i]]
      
      if(length(abx_days)>0){
        abx_exposed_individuals <-
          abx_exposed_individuals+1
        
        abx_exposure_episodes <-
          abx_exposure_episodes+
          length(abx_days)
      }
    }
  }
  
  ##############################
  # THRUST 1B
  ##############################
  
  index_episodes <- get_index_episodes(
    households
  )
  
  n_index_episodes <- nrow(
    index_episodes
  )
  
  # Number of unique pARG-positive individuals who receive at least one
  # priority antibiotic.
  if(n_index_episodes>0){
    n_index_individuals <- length(
      unique(
        paste(
          index_episodes$household,
          index_episodes$index,
          sep="_"
        )
      )
    )
  }else{
    n_index_individuals <- 0
  }
  
  # Five intensive samples per index individual:
  # days 0,3,7,14,28.
  intensive_samples <-
    n_index_individuals*
    length(index_times)
  
  ##############################
  # THRUST 1C
  ##############################
  
  pairs <- get_pairs(
    households
  )
  
  n_pairs <- nrow(
    pairs
  )
  
  n_events <- ifelse(
    n_pairs>0,
    sum(
      pairs$acquired
    ),
    0
  )
  
  # Two samples represented per eligible contact:
  # at index antibiotic exposure + next routine sample.
  contact_samples <- n_pairs*2
  
  ##############################
  # FIT HOUSEHOLD CTMC
  ##############################
  
  fit <- fit_markov(
    households
  )
  
  if(
    is.null(fit)
  ){
    return(
      data.frame(
        households=n_households,
        participants=n_households*hhsize,
        abx_exposed_individuals=
          abx_exposed_individuals,
        abx_exposure_episodes=
          abx_exposure_episodes,
        T1B_index_individuals=
          n_index_individuals,
        T1B_index_episodes=
          n_index_episodes,
        T1B_intensive_samples=
          intensive_samples,
        T1C_eligible_pairs=
          n_pairs,
        T1C_transmission_events=
          n_events,
        T1C_contact_samples=
          contact_samples,
        T1C_R0=NA,
        T1C_R0_lower=NA,
        T1C_R0_upper=NA,
        T1C_precision_50=NA,
        T1C_precision_25=NA,
        T1A_precision_bg=NA,
        T1A_precision_beta=NA,
        T1A_precision_gamma=NA,
        T1A_all_precision_30=NA
      )
    )
  }
  
  ##############################
  # R0 PRECISION
  ##############################
  
  R0_ci <- calculate_R0_ci(
    fit
  )
  
  R0_precision <- (
    R0_ci["upper"]-
      R0_ci["lower"]
  )/
    (
      2*R0_ci["R0"]
    )
  
  ##############################
  # PARAMETER PRECISION
  ##############################
  
  param_precision <-
    calculate_param_precision(
      fit
    )
  
  data.frame(
    households=n_households,
    participants=n_households*hhsize,
    abx_exposed_individuals=
      abx_exposed_individuals,
    abx_exposure_episodes=
      abx_exposure_episodes,
    T1B_index_individuals=
      n_index_individuals,
    T1B_index_episodes=
      n_index_episodes,
    T1B_intensive_samples=
      intensive_samples,
    T1C_eligible_pairs=
      n_pairs,
    T1C_transmission_events=
      n_events,
    T1C_contact_samples=
      contact_samples,
    T1C_R0=
      R0_ci["R0"],
    T1C_R0_lower=
      R0_ci["lower"],
    T1C_R0_upper=
      R0_ci["upper"],
    T1C_precision_50=
      R0_precision<=primary_R0_precision_target,
    T1C_precision_25=
      R0_precision<=secondary_R0_precision_target,
    T1A_precision_bg=
      as.numeric(
        param_precision$within_30["bg"]
      ),
    T1A_precision_beta=
      as.numeric(
        param_precision$within_30["beta"]
      ),
    T1A_precision_gamma=
      as.numeric(
        param_precision$within_30["gamma"]
      ),
    T1A_all_precision_30=
      all(
        param_precision$within_30
      )
  )
}

##############################
# SAMPLE SIZE SIMULATIONS
##############################

results <- data.frame()

for(hh in household_grid){
  
  cat(
    "\nHouseholds:",
    hh,
    "\n"
  )
  
  sim_results <- vector(
    "list",
    numsims
  )
  
  for(sim in 1:numsims){
    print(sim)
    if(
      sim%%10==0
    ){
      cat(
        "simulation",
        sim,
        "of",
        numsims,
        "\n"
      )
    }
    
    sim_results[[sim]] <-
      run_one_study(
        hh
      )
  }
  
  d <- do.call(
    rbind,
    sim_results
  )
  
  ##############################
  # THRUST 1B YIELD
  ##############################
  
  # Probability of generating at least 40,60,80 or 100
  # pARG-positive antibiotic-exposed index episodes.
  index_threshold_results <- sapply(
    index_thresholds,
    function(x){
      mean(
        d$T1B_index_episodes>=x,
        na.rm=TRUE
      )
    }
  )
  
  names(
    index_threshold_results
  ) <-
    paste0(
      "T1B_probability_indices_",
      index_thresholds
    )
  
  ##############################
  # SAMPLE-SIZE SUMMARY
  ##############################
  
  results <- rbind(
    results,
    data.frame(
      households=hh,
      participants=hh*hhsize,
      
      # Antibiotic exposure burden.
      mean_abx_exposed_individuals=
        mean(
          d$abx_exposed_individuals,
          na.rm=TRUE
        ),
      
      mean_abx_exposure_episodes=
        mean(
          d$abx_exposure_episodes,
          na.rm=TRUE
        ),
      
      # Thrust 1B.
      mean_T1B_index_individuals=
        mean(
          d$T1B_index_individuals,
          na.rm=TRUE
        ),
      
      mean_T1B_index_episodes=
        mean(
          d$T1B_index_episodes,
          na.rm=TRUE
        ),
      
      T1B_probability_indices_40=
        index_threshold_results[
          "T1B_probability_indices_40"
        ],
      
      T1B_probability_indices_60=
        index_threshold_results[
          "T1B_probability_indices_60"
        ],
      
      T1B_probability_indices_80=
        index_threshold_results[
          "T1B_probability_indices_80"
        ],
      
      T1B_probability_indices_100=
        index_threshold_results[
          "T1B_probability_indices_100"
        ],
      
      # Intensive within-host sample burden.
      mean_T1B_intensive_samples=
        mean(
          d$T1B_intensive_samples,
          na.rm=TRUE
        ),
      
      # Thrust 1C.
      mean_T1C_eligible_pairs=
        mean(
          d$T1C_eligible_pairs,
          na.rm=TRUE
        ),
      
      mean_T1C_transmission_events=
        mean(
          d$T1C_transmission_events,
          na.rm=TRUE
        ),
      
      mean_T1C_contact_samples=
        mean(
          d$T1C_contact_samples,
          na.rm=TRUE
        ),
      
      mean_T1C_R0=
        mean(
          d$T1C_R0,
          na.rm=TRUE
        ),
      
      mean_T1C_R0_lower=
        mean(
          d$T1C_R0_lower,
          na.rm=TRUE
        ),
      
      mean_T1C_R0_upper=
        mean(
          d$T1C_R0_upper,
          na.rm=TRUE
        ),
      
      # Primary Thrust 1C criterion.
      T1C_probability_precision_50=
        mean(
          d$T1C_precision_50,
          na.rm=TRUE
        ),
      
      # Optional stricter sensitivity.
      T1C_probability_precision_25=
        mean(
          d$T1C_precision_25,
          na.rm=TRUE
        ),
      
      # Thrust 1A: parameter-level +/-30% precision.
      T1A_probability_precision_bg=
        mean(
          d$T1A_precision_bg,
          na.rm=TRUE
        ),
      
      T1A_probability_precision_beta=
        mean(
          d$T1A_precision_beta,
          na.rm=TRUE
        ),
      
      T1A_probability_precision_gamma=
        mean(
          d$T1A_precision_gamma,
          na.rm=TRUE
        ),
      
      # Probability that ALL THREE parameters meet +/-30%.
      T1A_probability_all_parameters_30=
        mean(
          d$T1A_all_precision_30,
          na.rm=TRUE
        ),
      
      # Number of successful model fits.
      successful_fits=
        sum(
          !is.na(
            d$T1C_R0
          )
        ),
      
      fit_success_rate=
        mean(
          !is.na(
            d$T1C_R0
          )
        ),
      
      # Routine sampling burden.
      routine_samples_per_person=
        length(
          routine_times
        ),
      
      routine_samples_per_household=
        length(
          routine_times
        )*hhsize,
      
      routine_samples_total=
        hh*hhsize*
        length(
          routine_times
        ),
      
      # Approximate additional samples from intensive index follow-up
      # and contact sampling.
      mean_total_additional_samples=
        mean(
          d$T1B_intensive_samples+
            d$T1C_contact_samples,
          na.rm=TRUE
        )
    )
  )
}

##############################
# DISPLAY RESULTS
##############################

print(
  results
)

##############################
# SAVE ASSUMPTIONS AND RESULTS
##############################

# The Excel file stores:
# Assumptions = exact parameters used
# Results = simulation results for each candidate cohort size.
assumptions <- data.frame(
  parameter=c(
    "site",
    "household_size",
    "follow_up_days",
    "routine_sampling_days",
    "index_sampling_days",
    "routine_samples_per_person",
    "intensive_samples_per_index",
    "baseline_prevalence",
    "overall_acquisition_rate",
    "background_fraction",
    "background_rate",
    "p2p_rate",
    "mean_carriage_days",
    "clearance_rate",
    "abx_courses_per_year",
    "priority_abx_fraction",
    "abx_capture_probability",
    "primary_R0_precision_target",
    "secondary_R0_precision_target",
    "parameter_precision_target"
  ),
  value=c(
    site,
    hhsize,
    follow_up_days,
    paste(
      routine_times,
      collapse=", "
    ),
    paste(
      index_times,
      collapse=", "
    ),
    length(
      routine_times
    ),
    length(
      index_times
    ),
    baseline_prevalence,
    overall_acquisition_rate,
    background_fraction,
    background_rate,
    p2p_rate,
    mean_carriage_days,
    clearance_rate,
    abx_courses_per_year,
    priority_abx_fraction,
    abx_capture_probability,
    primary_R0_precision_target,
    secondary_R0_precision_target,
    parameter_precision_target
  )
)

if(
  !dir.exists(
    "./Output"
  )
){
  dir.create(
    "./Output",
    recursive=TRUE
  )
}

write_xlsx(
  list(
    Assumptions=assumptions,
    Results=results
  ),
  "./Output/PREEMPT_sample_size_results_hh7.xlsx"
)
