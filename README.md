# Sample size simulation scripts

This folder contains R scripts used for design-stage sample size and feasibility calculations for PRE-EMPT AMR (Wellcome LEAP Application - 2608B9002/DC7).

## `preempt_Thrust1B_sample_size_sim.R`

Simulation for the **UK Thrust 1B controlled antibiotic challenge**.

The script evaluates how many participants need to be screened to obtain sufficient pARG-positive challenge participants and detectable plasmid-host acquisition events.

It varies assumptions including:

* baseline pARG prevalence;
* challenge participation;
* acquisition-event rate; and
* probability of detecting an acquisition event.

The number of detected events is used as a design-stage proxy for the amount of information available for estimating within-gut plasmid transmission (`R0-within`). This is a sample-size approximation and not the final mechanistic analysis.

Results are exported to `Output/Thrust1b/UK_Thrust1B_simple_R0within_sample_size.xlsx`.

## `preempt_Thrust1C_sample_size_sim.R`

Simulation for the **PRE-EMPT household cohort**.

The script simulates longitudinal pARG carriage within households using a continuous-time Markov model (CTMC), allowing for:

* background acquisition;
* person-to-person household transmission; and
* clearance of carriage.

For each candidate cohort size, simulated routine household observations are used to re-estimate the model parameters and assess their precision.

The script evaluates:

* precision of the household reproduction number (`R0`);
* precision of acquisition, transmission (`beta`), and clearance (`gamma`) parameters;
* expected numbers of informative antibiotic-exposed index participants; and
* expected index-contact pairs and observed acquisitions.

## Running

From the repository root:

```bash
Rscript Rscripts/preempt_Thrust1B_sample_size_sim.R
Rscript Rscripts/preempt_Thrust1C_sample_size_sim.R
```

The simulations are intended for **study design and sensitivity analysis**. Results depend on the epidemiological, biological, recruitment, and measurement assumptions specified near the beginning of each script.
