# seabird-causal Code

## Recommended code workflow and script overview

1.  [*A_kernels_extract.R*](A_kernels_extract.R) - from light logger data, estimate the position of birds during the breeding and non-breeding seasons. Using these positions, 'core areas' are estimated for each population (the 50th percentile contours) which are taken forward to estimate stressor values.

2.  [*B_download_environmentals.R*](B_download_environmentals.R) - download air temperature, and sea surface temperature (SST) data from ERA5 (<https://doi.org/10.24381/cds.4991cf48>).

3.  [*C_sst_extract.R*](C_sst_extract.R) - using ERA5 data (<https://doi.org/10.24381/cds.4991cf48>), extracts mean and variation in sea surface temperature for a given 'core area' polygon over the date range of bird data.

4.  [*D_airtemp_extract.R*](D_airtemp_extract.R) - using ERA5 data (<https://doi.org/10.24381/cds.4991cf48>), extracts mean air temperature at the colony coordinates over the date range of bird data.

5.  [*E_gls_classification.R*](E_gls_classification.R) - from light and immersion logger data, estimate the proportion of time birds spend performing certain behaviours (in flight, foraging, resting, at colony).

6.  [*F_calc_energetics.R*](F_calc_energetics.R) - following Frankish et al. (REF), combine behaviourial classfications with SST and air temperature to estimate daily energetic expenditure.

7.  [*G_TRIM.R*](G_TRIM.R) - estimate population abundance using TRIM (TRends and Indices for Monitoring data - <https://www.cbs.nl/en-gb/society/nature-and-environment/indices-and-trends--trim-->) following Knape 2016 (<https://doi.org/10.1111/1365-2664.12720>).

8.  [*H_clean_data.R*](H_clean_data.R) - data loading and cleaning. Involves changing of variable names to be consistent across research sites.

9.  [*I_population_RD_stability_behaviour.R*](I_population_RD_stability_behaviour.R) - from the [cleaned](G_clean_data.R) data, fit behaviourial reaction norms via generalised additive models (GAM) and estimate response diversity from GAM derived first derivatives.

10. [*J_population_RD_stability_energy.R*](J_population_RD_stability_energy.R) - from the [cleaned](G_clean_data.R) data, fit energetic expenditure reaction norms via generalised additive models (GAM) and estimate response diversity from GAM derived first derivatives.

11. [*K_figure2.R*](K_figure2.R) - generate figure2.png.

12. [*L_stability_rd_spp.R*](L_stability_rd_spp.R) - correlate the stability of each species' colony with its response diversity estimate from the GAMs fit by [*H_population_RD_stability_behaviour.R*](H_population_RD_stability_behaviour.R) and [*I_population_RD_stability_energy.R*](I_population_RD_stability_energy.R).

13. [*M_figure3.R*](M_figure3.R) - generate figure3.png.

## Supporting functions

The [utilities folder](https://github.com/duncanobrien/seabird-RD/tree/master/R/utilities) contains supporting functions used in the primary workflow scripts.
