# CLOVERS Study: Data Analysis

- `create_longitudinal_outcome.R`: code to create daily states from the CLOVERS data

- `oddsratio28days.R`: code to estimate the odds ratio at the end of follow-up for the association between treatment and mortality under the 12 estimation procedures considered in the paper (Figure 3).

-  `riskdifference.R`: code to estimate the risk differences at the end of follow-up for the association between treatment and mortality (Figure 3). The code additionally includes estimation via survival analyses using both a Cox proportional hazard model and a multistate model.

-  `plotclovers.R`: code to plot the CLOVERS outcome data (Figure 1) and the results of each estimation procedure (Figure 3). To plot the results from each estimation procedure, one would need to 1) run `oddsratio28days.R` and `riskdifference.R`, and 2) take the saved results.

To be able to reproduce the results one would need the original data from the CLOVERS study and run the code in `oddsratio28days.R` first, followed by `oddsratio28days.R` and `riskdifference.R`.
