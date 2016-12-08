# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/) 
and this project adheres to [Semantic Versioning](http://semver.org/).

## 0.0.7 - 2016-12-08
### Added
- A new activity domain computation method is added, "compute_wake_sleep_domains_3.m". This method, as opposed to the previous ones, does not rely on contor lines being connected and anything like that. It searches the heat map from bottom to top and splits into two clusters using MatLab's built-in k-means function. If clusters are distant enough, method stops and identifies activity domains.

### Changed
- Minor changes in "generate_random_graph.m".

## 0.0.6 - 2016-12-07
### Fixed
- The function name in "update_SW_network.m" is corrected.

### Changed
- Event probabilities are changed for the two-state neuron for simulation time improvement.
- addpath ../../power-law-estimator removed from "scale_free_graph_generation.m". It was taking too much time when called too many times in a loop.
- Minor changes in "generate_random_graph.m".

## 0.0.5 - 2016-12-06
### Fixed
- Fixed zero degree removal in "scale_free_graph_generation.m".
- Fixed missing option in "generate_random_graph.m" for SF-Chung-Lu option.

### Changed
- Mean degree verification now comes before graph connectedness verification to improve simulation time.
- Counter is added to "generate_random_graph.m" for 'ER' and 'SF-Chung-Lu' options.
- Minor changes in the main file "SW_network.m". New input options are added.

## 0.0.4 - 2016-12-05
### Changed
- Minor changes in Chung-Lu.

### Added
- A new option is added to "generate_inhibitory_links". Previously, power-law exponent determination was weird. Now, it is more systematic.

## 0.0.3 - 2016-12-01
### Changed
- Minor changes

### Added
- "plot_SE_WE_heat_map.m" is added. All related parts of "SE_network.m" are moved there.
- "plot_sample_SE_WE_trajectories.m" is added. All related parts of "SW_network.m" are moved there.

## 0.0.2 - 2016-11-30
### Changed
- Minor changes in "generate_random_graph.m"

### Added
- "do_SW_drift_diffusion_analysis.m" is added. All related parts of "SW_network.m" are moved there.
- "fit_PL_and_EXP_to_SW_bouts.m" is added. All parts of "SW_network.m" related to bounded power-law and exponential distribution fitting to SW bout durations are moved there. 
- "interval_to_bout.m" was in the RBM folder. It is now copied into SW-MatLab and renamed to "compute_bout_durations.m".
- "scale_free_graph_generation.m", "determine_alpha.m" and "chung_lu.m" are moved into a new folder titled "generate_network_configuration/SF_Chung_Lu".
- "generate_random_graph.m" and "generate_inhibitory_links.m" are moved into a new folder titled "generate_network_configuration".
- "compute_wake_sleep_domains_NEW.m", "compute_wake_sleep_domains.m", "compute_duration_matrix.m", "calculate_ints_tpb.m" are moved into a new folder titled "compute_activity_domains"
- "compute_activity_domains.m" is added. All related parts of "SW_network.m are moved there.
- "record_sim.m" is added. All parts of "SW_network.m" that are related to recording simulation output in a single structure are moved there.
- "generate_network_configuration.m" is added. All graph generation related parts of "SW_network.m" are now moved there.

## 0.0.1 - 2016-11-30
### Added
- This CHANGELOG file will be updated as my MatLab code for the Sleep-Wake project is updated.

